"""Pair continuous-heading irregular passive yaw with pinned MATLAB output."""

import os
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest
from scipy.interpolate import CubicSpline
from scipy.signal import fftconvolve

from wecsim import PMWave, WEC, WorldPoint
from wecsim.bodyClass import BodyClass
from wecsim.irregularWave import pm_equal_energy_components
from wecsim.passiveYaw import (
    HeldPassiveYawExcitation, SampledPassiveYawExcitation,
)


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="paired MATLAB irregular passive-yaw output not provided",
)


def _max_error(actual, expected, limit, label):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = np.max(np.abs(actual - expected))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


def test_published_irregular_passive_yaw_with_source_force(tmp_path):
    """Replay the published forcing to isolate wave and motion integration."""
    source = Path(REFERENCE)
    hydro = (Path(APPLICATIONS)
             / "_Common_Input_Files/OSWEC/hydroData/oswec.h5").resolve()
    prefix = "OSWEC_PASSIVE_YAW_IRR"
    case = "PassiveYawRegression"
    components_csv = np.loadtxt(source / f"{prefix}_components.csv", delimiter=",")
    flap = np.loadtxt(source / f"{prefix}_{case}_body1.csv", delimiter=",")
    base = np.loadtxt(source / f"{prefix}_{case}_body2.csv", delimiter=",")
    pto = np.loadtxt(source / f"{prefix}_{case}_pto1.csv", delimiter=",")
    wave = np.loadtxt(source / f"{prefix}_wave.csv", delimiter=",")
    assert components_csv.shape == (500, 4)
    assert flap.shape == base.shape == pto.shape == (25001, 25)
    assert wave.shape == (25001, 2)

    phases = tmp_path / "phases.csv"
    np.savetxt(phases, components_csv[:, 3, None], delimiter=",")
    components = pm_equal_energy_components(
        hydro, significant_height=2.5, peak_period=8,
        directions=[10], spreading=[1], phase=components_csv[:, 3, None],
    )
    for actual, expected, label in (
        (components.omega, components_csv[:, 0], "frequency"),
        (components.spectral_amplitude, components_csv[:, 1], "spectrum"),
        (components.d_omega, components_csv[:, 2], "frequency width"),
    ):
        _max_error(actual, expected, 1e-12, label)

    body = BodyClass(str(hydro))
    body.bodyNumber = 1
    body.bodyTotal = 2
    body.readH5file()
    sampled = SampledPassiveYawExcitation.from_hydro_data(
        body.hydroData, components, dt=0.01, end_time=250,
        ramp_time=100, rho=1000, g=9.81,
    )
    _max_error(sampled.elevation, wave[:, 1], 1e-11, "wave elevation")

    class LoggedSourceForce:
        elevation = sampled.elevation

        def force(self, at_time, yaw):
            return flap[round(at_time / 0.01), 19:25]

    wec = WEC("OSWEC published irregular passive yaw")
    moving = wec.body("flap", hydro, mass=12700,
                      inertia=(1.85e6,) * 3, passive_yaw=True)
    wec.body("base", hydro, mass=999, inertia=(999,) * 3)
    yaw = wec.coordinate(
        "yaw", moving.move("yaw", pivot=WorldPoint(0, 0, -8.9)),
    )
    wec.rotational_pto("hinge", yaw, damping=120000)
    # This is a dynamics diagnostic. The published MATLAB run holds its
    # heading coefficients for 1 degree; the Python default stays continuous.
    with patch.object(SampledPassiveYawExcitation, "from_hydro_data",
                      return_value=LoggedSourceForce()):
        result = wec.run(PMWave(2.5, 8, direction=10, phase_file=phases),
                         dt=0.01, end_time=250, ramp_time=100,
                         radiation_memory=40)
    _max_error(result.time, flap[:, 0], 1e-10, "time")
    _max_error(result.wave_elevation, wave[:, 1], 1e-11, "wave elevation")
    _max_error(result.bodies["flap"].position[:, 5], flap[:, 6],
               0.003, "source-forced flap yaw")
    _max_error(result.bodies["flap"].velocity[:, 5], flap[:, 12],
               0.0002, "source-forced yaw speed")
    _max_error(result.bodies["base"].position, base[:, 1:7],
               1e-10, "fixed base position")
    _max_error(result.bodies["base"].velocity, base[:, 7:13],
               1e-10, "fixed base velocity")
    hinge = result.ptos["hinge"]
    _max_error(hinge.force, pto[:, 17], 25, "source-forced PTO torque")
    _max_error(-hinge.absorbed_power, pto[:, 23], 4,
               "source-forced source-signed PTO power")
    _max_error(pto[:, 17], -120000 * pto[:, 11], 1e-5,
               "source PTO damping law")


def test_published_irregular_passive_yaw_heading_threshold(tmp_path):
    """Check the published 1-degree hold, force balance, and native motion."""
    source = Path(REFERENCE)
    hydro = (Path(APPLICATIONS)
             / "_Common_Input_Files/OSWEC/hydroData/oswec.h5").resolve()
    prefix = "OSWEC_PASSIVE_YAW_IRR"
    flap = np.loadtxt(source / f"{prefix}_PassiveYawRegression_body1.csv",
                      delimiter=",")
    pto = np.loadtxt(source / f"{prefix}_PassiveYawRegression_pto1.csv",
                     delimiter=",")
    forces = np.loadtxt(source / f"{prefix}_body1_forces.csv", delimiter=",")
    mass = np.loadtxt(source / f"{prefix}_body1_mass.csv", delimiter=",")
    components_csv = np.loadtxt(source / f"{prefix}_components.csv", delimiter=",")
    phases = tmp_path / "phases.csv"
    np.savetxt(phases, components_csv[:, 3, None], delimiter=",")
    components = pm_equal_energy_components(
        hydro, significant_height=2.5, peak_period=8,
        directions=[10], spreading=[1], phase=components_csv[:, 3, None],
    )
    body = BodyClass(str(hydro))
    body.bodyNumber = 1
    body.bodyTotal = 2
    body.readH5file()
    assert forces.shape == flap.shape == (25001, 25)
    assert mass.shape == (4,)
    _max_error(forces[:, 0], flap[:, 0], 1e-10, "force sample time")
    _max_error(mass, [12700, 1.85e6, 1.85e6, 1.85e6], 1e-9,
               "published flap mass and inertia")
    irf = body.hydroData["hydro_coeffs"]["radiation_damping"][
        "impulse_response_fun"]
    kernel = CubicSpline(
        np.asarray(irf["t"]).ravel(),
        np.asarray(irf["K"])[5, 5] * 1000,
    )(np.arange(4001) * 0.01)
    speed = flap[:, 12]
    radiation = 0.01 * fftconvolve(speed, kernel)[:len(speed)]
    radiation -= 0.005 * kernel[0] * speed
    radiation[4000:] -= 0.005 * kernel[-1] * speed[:-4000]
    _max_error(forces[:, 6], radiation, 1e-6,
               "source yaw radiation from saved speed and HDF5 IRF")
    added_mass = body.hydroData["hydro_coeffs"]["added_mass"][
        "inf_freq"][5, 5] * 1000
    _max_error(forces[:, 12], added_mass * forces[:, 24], 1e-6,
               "source yaw added mass from saved acceleration")
    _max_error(forces[:, 18], np.zeros(len(forces)), 1e-8,
               "published yaw hydrostatic restoring")
    _max_error(flap[:, 18], flap[:, 24] - forces[:, 6]
               - forces[:, 12] - forces[:, 18], 1e-6,
               "source yaw hydrodynamic force balance")
    _max_error(flap[:, 18] + pto[:, 17], mass[3] * forces[:, 24],
               1e-6, "source yaw rigid inertia balance")
    model = SampledPassiveYawExcitation.from_hydro_data(
        body.hydroData, components, dt=.01, end_time=250,
        ramp_time=100, rho=1000, g=9.81,
    )
    held = HeldPassiveYawExcitation(model, threshold=1)
    for at_time, angle in zip(flap[:, 0], flap[:, 6]):
        held.commit(at_time, np.array([angle]))
    _max_error(np.asarray(held.force_history), flap[:, 19:25], 1e-6,
               "published sampled heading force on MATLAB yaw")

    wec = WEC("OSWEC sampled irregular passive yaw")
    moving = wec.body(
        "flap", hydro, mass=12700, inertia=(1.85e6,) * 3,
        passive_yaw=True, passive_yaw_threshold=1,
    )
    wec.body("base", hydro, mass=999, inertia=(999,) * 3)
    yaw = wec.coordinate(
        "yaw", moving.move("yaw", pivot=WorldPoint(0, 0, -8.9)),
    )
    wec.rotational_pto("hinge", yaw, damping=120000)
    result = wec.run(
        PMWave(2.5, 8, direction=10, phase_file=phases),
        dt=.01, end_time=250, ramp_time=100, radiation_memory=40,
    )
    _max_error(result.time, flap[:, 0], 1e-10, "time")
    _max_error(result.bodies["flap"].position[:5001, 5], flap[:5001, 6],
               1e-5, "motion before first threshold crossing")
    _max_error(result.bodies["flap"].position[:, 5], flap[:, 6], .21,
               "sampled-heading trajectory sensitivity envelope")
    _max_error(result.bodies["flap"].velocity[:, 5], flap[:, 12], .03,
               "sampled-heading speed sensitivity envelope")
    _max_error(result.ptos["hinge"].force, -120000 * result.ptos["hinge"].velocity,
               1e-8, "native PTO damping law")
    assert pto.shape == (25001, 25)


def test_continuous_heading_irregular_passive_yaw(tmp_path):
    """Check the same PM realization and 250 s flap/PTO trajectory."""
    source = Path(REFERENCE)
    hydro = (Path(APPLICATIONS)
             / "_Common_Input_Files/OSWEC/hydroData/oswec.h5").resolve()
    prefix = "OSWEC_PASSIVE_YAW_IRR_CONT"
    case = "PassiveYawRegressionContinuous"
    components_csv = np.loadtxt(source / f"{prefix}_components.csv", delimiter=",")
    flap = np.loadtxt(source / f"{prefix}_{case}_body1.csv", delimiter=",")
    base = np.loadtxt(source / f"{prefix}_{case}_body2.csv", delimiter=",")
    pto = np.loadtxt(source / f"{prefix}_{case}_pto1.csv", delimiter=",")
    wave = np.loadtxt(source / f"{prefix}_wave.csv", delimiter=",")
    assert components_csv.shape == (500, 4)
    assert flap.shape == base.shape == pto.shape == (25001, 25)
    assert wave.shape == (25001, 2)

    phases = tmp_path / "phases.csv"
    np.savetxt(phases, components_csv[:, 3, None], delimiter=",")
    components = pm_equal_energy_components(
        hydro, significant_height=2.5, peak_period=8,
        directions=[10], spreading=[1], phase=components_csv[:, 3, None],
    )
    for actual, expected, label in (
        (components.omega, components_csv[:, 0], "frequency"),
        (components.spectral_amplitude, components_csv[:, 1], "spectrum"),
        (components.d_omega, components_csv[:, 2], "frequency width"),
    ):
        _max_error(actual, expected, 1e-12, label)

    body = BodyClass(str(hydro))
    body.bodyNumber = 1
    body.bodyTotal = 2
    body.readH5file()
    force_model = SampledPassiveYawExcitation.from_hydro_data(
        body.hydroData, components, dt=0.01, end_time=250,
        ramp_time=100, rho=1000, g=9.81,
    )
    source_yaw_force = np.stack([
        force_model.force(t, angle)
        for t, angle in zip(flap[:, 0], flap[:, 6])
    ])
    _max_error(source_yaw_force, flap[:, 19:25], 1e-6,
               "six-component excitation on source yaw")

    wec = WEC("OSWEC irregular passive yaw")
    moving = wec.body("flap", hydro, mass=12700,
                      inertia=(1.85e6,) * 3, passive_yaw=True)
    wec.body("base", hydro, mass=999, inertia=(999,) * 3)
    yaw = wec.coordinate(
        "yaw", moving.move("yaw", pivot=WorldPoint(0, 0, -8.9)),
    )
    wec.rotational_pto("hinge", yaw, damping=120000)
    result = wec.run(PMWave(2.5, 8, direction=10, phase_file=phases),
                     dt=0.01, end_time=250, ramp_time=100,
                     radiation_memory=40)
    _max_error(result.time, flap[:, 0], 1e-10, "time")
    _max_error(result.wave_elevation, wave[:, 1], 1e-12, "wave elevation")
    _max_error(result.bodies["flap"].position[:, 5], flap[:, 6],
               1e-4, "flap yaw")
    _max_error(result.bodies["flap"].velocity[:, 5], flap[:, 12],
               2e-5, "flap yaw speed")
    _max_error(result.bodies["flap"].position[:, :5], flap[:, 1:6],
               1e-10, "stationary flap positions")
    _max_error(result.bodies["flap"].velocity[:, :5], flap[:, 7:12],
               1e-10, "stationary flap velocities")
    _max_error(result.bodies["base"].position, base[:, 1:7],
               1e-10, "fixed base position")
    _max_error(result.bodies["base"].velocity, base[:, 7:13],
               1e-10, "fixed base velocity")
    hinge = result.ptos["hinge"]
    _max_error(hinge.stroke, pto[:, 5], 1e-4, "PTO angle")
    _max_error(hinge.velocity, pto[:, 11], 2e-5, "PTO angular speed")
    _max_error(hinge.force, pto[:, 17], 3, "PTO torque")
    _max_error(-hinge.absorbed_power, pto[:, 23], 0.05,
               "PTO source-signed power")
    _max_error(pto[:, 17], -120000 * pto[:, 11], 1e-5,
               "source PTO damping law")
