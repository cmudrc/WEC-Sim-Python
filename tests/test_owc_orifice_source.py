"""Pair the published OWC orifice block with its saved MATLAB signal path."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat
from scipy.signal import fftconvolve

from wecsim import OrificePTO
from wecsim.bodyClass import BodyClass
from wecsim.irregularWave import (
    pm_equal_energy_components, synthesize_irregular_response,
)


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="fresh MATLAB OWC orifice output absent",
)


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


def _load(name):
    return np.loadtxt(Path(REFERENCE) / f"OWC_ORIFICE_{name}.csv",
                      delimiter=",", ndmin=2)


def test_published_orifice_force_on_matlab_piston_motion():
    reference = Path(REFERENCE)
    logged = _load("orifice")
    parameters = _load("parameters").ravel()
    flex = loadmat(reference / "OWC_ORIFICE_Flex_out.mat",
                   simplify_cells=True)["Flex_out"]
    time = np.asarray(flex["time"])
    values = np.asarray(flex["signals"]["values"])
    assert logged.ndim == 2 and logged.shape[1] == 6
    assert values.ndim == 2 and values.shape[1] == 10
    assert len(time) == len(values) and np.all(np.diff(time) > 0)
    assert np.all(np.diff(logged[:, 0]) > 0)
    model = OrificePTO(*parameters)
    # The orifice block logs every ode23t internal evaluation (418,828 rows),
    # while Flex_out is sampled at 0.005 s (26,001 rows). The logged flow gives
    # the exact internal piston speed. Its interpolated value agrees with the
    # independent flexible-mode speed on the output grid.
    piston_speed = logged[:, 4] / model.piston_area
    sampled_speed = np.interp(time, logged[:, 0], piston_speed)
    _max_error(sampled_speed, values[:, 1], 1e-10,
               "orifice input versus flexible-mode speed")
    result = model.evaluate(piston_speed)

    # MATLAB logs [force N, Mach flag, pressure kPa, flow m^3/s, power kW].
    _max_error(result.force, logged[:, 1], 1e-8, "orifice reaction")
    np.testing.assert_array_equal(result.compressibility_flag, logged[:, 2])
    _max_error(result.pressure_drop / 1000, logged[:, 3], 1e-10,
               "orifice pressure")
    _max_error(result.flow_rate, logged[:, 4], 1e-13, "orifice flow")
    _max_error(result.absorbed_power / 1000, logged[:, 5], 1e-10,
               "orifice power")
    _max_error(result.absorbed_power, -result.force * piston_speed,
               1e-8, "passive orifice energy")


@pytest.mark.skipif(not APPLICATIONS,
                    reason="BEMIO-generated OWC HDF5 absent")
def test_published_owc_wave_and_seventh_excitation_channel():
    hydro = (Path(APPLICATIONS)
             / "OWC/OrificeModel/hydroData/test17a_clean.h5")
    components = _load("components")
    wave = _load("wave")
    rigid = _load("OrificeModel_body1")
    flex = loadmat(Path(REFERENCE) / "OWC_ORIFICE_Flex_out.mat",
                   simplify_cells=True)["Flex_out"]
    flexible = np.asarray(flex["signals"]["values"])
    assert components.shape == (500, 4)
    assert wave.shape == (26001, 2)
    assert rigid.shape == (26001, 25)
    assert flexible.shape == (26001, 10)

    realized = pm_equal_energy_components(
        hydro, significant_height=1, peak_period=4,
        directions=np.array([0.0]), spreading=np.array([1.0]),
        phase=components[:, 3, None],
    )
    _max_error(realized.omega, components[:, 0], 1e-12, "PM frequency")
    _max_error(realized.spectral_amplitude, components[:, 1], 1e-12,
               "PM spectral amplitude")
    _max_error(realized.d_omega, components[:, 2], 1e-12, "PM bin width")
    result = synthesize_irregular_response(
        hydro, realized, dt=0.005, end_time=130, ramp_time=10,
    )
    assert result.excitation_force.shape == (26001, 7)
    _max_error(result.time, wave[:, 0], 1e-10, "wave time")
    _max_error(result.elevation, wave[:, 1], 1e-11, "wave elevation")
    _max_error(result.excitation_force[:, :6], rigid[:, 19:25], 1e-7,
               "rigid excitation")
    _max_error(result.excitation_force[:, 6], flexible[:, 4], 1e-7,
               "flexible excitation")

    body = BodyClass(str(hydro))
    body.bodyNumber = body.bodyTotal = 1
    body.readH5file()
    body.mass = "equilibrium"
    memory = np.arange(3001) * 0.005
    body.hydroForcePre(
        components[:, 0], [0], len(memory), memory, len(components),
        0.005, 1000, 9.81, "irregular", np.zeros((2, 1)), 1, 1, 0, 0, 0,
    )
    force = body.hydroForce
    for key, filename in (
        ("mass_ff", "flexible_effective_mass"),
        ("stiffness", "flexible_stiffness"),
        ("damping", "flexible_damping"),
    ):
        _max_error(np.asarray(force["gbm"][key]), _load(filename), 1e-9,
                   f"flexible {key}")
    _max_error(np.asarray(force["fAddedMass"]), _load("added_mass"),
               1e-9, "added mass")
    _max_error(np.asarray(force["linearHydroRestCoef"]),
               _load("hydrostatic_stiffness"), 1e-9,
               "hydrostatic stiffness")

    # The source flexible state equation applies the orifice reaction in
    # addition to its logged hydrodynamic force balance.
    total = flexible[:, 3]
    hydro_terms = flexible[:, 4:]
    _max_error(total, hydro_terms[:, 0] - hydro_terms[:, 1:].sum(axis=1),
               1e-8, "source flexible hydrodynamic force balance")
    source_orifice = _load("orifice")
    reaction = np.interp(wave[:, 0], source_orifice[:, 0], source_orifice[:, 1])
    mass = float(np.asarray(force["gbm"]["mass_ff"])[0, 0])
    _max_error(mass * flexible[:, 2], total + reaction, 1e-6,
               "source flexible acceleration and orifice balance")

    # Reconstruct the seventh radiation channel from all seven source speeds.
    # The finite-memory convolution is trapezoidal at its first and last
    # samples, matching the source's continuous convolution block.
    velocities = np.column_stack((rigid[:, 7:13], flexible[:, 1]))
    kernel = np.asarray(force["irkb"])[:, 6, :]
    radiation = 0.005 * sum(
        fftconvolve(kernel[:, channel], velocities[:, channel])[:len(wave)]
        for channel in range(7)
    )
    radiation -= 0.0025 * (velocities * kernel[0]).sum(axis=1)
    last_lag = len(kernel) - 1
    radiation[last_lag:] -= 0.0025 * (
        velocities[:-last_lag] * kernel[-1]
    ).sum(axis=1)
    _max_error(radiation, flexible[:, 5], 0.01,
               "source seventh radiation force")
