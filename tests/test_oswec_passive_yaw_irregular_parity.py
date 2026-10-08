"""Pair continuous-heading irregular passive yaw with pinned MATLAB output."""

import os
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from wecsim import PMWave, WEC, WorldPoint
from wecsim.bodyClass import BodyClass
from wecsim.irregularWave import pm_equal_energy_components
from wecsim.passiveYaw import SampledPassiveYawExcitation


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
