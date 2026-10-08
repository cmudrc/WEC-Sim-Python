"""Pair a fixed monopile with current against pinned MATLAB WEC-Sim."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import Current, PMWave, WEC
from wecsim.irregularWave import pm_equal_energy_components
from wecsim.morison import finite_depth_wavenumber


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")


def _check(actual, expected, tolerance, label):
    actual, expected = np.asarray(actual), np.asarray(expected)
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = np.max(np.abs(actual - expected))
    assert error < tolerance, f"{label}: {error:.6g} exceeds {tolerance}"


@pytest.mark.skipif(not REFERENCE, reason="paired MATLAB current output not provided")
def test_fixed_monopile_power_law_current_against_matlab(tmp_path):
    source = Path(REFERENCE)
    components = np.loadtxt(source / "MORISON_FIXED_CURRENT_components.csv",
                            delimiter=",")
    wave = np.loadtxt(source / "MORISON_FIXED_CURRENT_wave.csv", delimiter=",")
    force = np.loadtxt(source / "MORISON_FIXED_CURRENT_morison_force.csv",
                       delimiter=",")
    monopile = np.loadtxt(source / "MORISON_FIXED_CURRENT_morisonElement_body1.csv",
                          delimiter=",")
    assert components.shape == (500, 5)
    assert wave.shape == (2001, 2)
    assert force.shape == (2001, 7)
    assert monopile.shape == (2001, 25)
    _check(monopile[:, 13:19], -force[:, 1:7], 1e-8,
           "MATLAB total and Morison force sign")

    phases = tmp_path / "phases.csv"
    np.savetxt(phases, components[:, 4:5], delimiter=",")
    sea = pm_equal_energy_components(
        None, significant_height=2, peak_period=5,
        directions=(30,), spreading=(1,), frequency_range=(.001, 10),
        phase=components[:, 4:5],
    )
    _check(sea.omega, components[:, 0], 1e-12, "frequency")
    _check(sea.spectral_amplitude, components[:, 1], 1e-12, "PM spectrum")
    _check(sea.d_omega, components[:, 2], 1e-12, "frequency width")
    _check(finite_depth_wavenumber(sea.omega, water_depth=30),
           components[:, 3], 1e-12, "wavenumber")

    wec = WEC("fixed monopile with current")
    body = wec.fixed_body(
        "monopile", center_gravity=(0, 0, -15),
        mass="equilibrium", inertia=(1.25e9, 1.25e9, .15e9),
        volume=np.pi * 10**2 * 30,
    )
    wec.fixed_body(
        "tower", center_gravity=(0, 0, 25), mass=1031930,
        inertia=(9.66e8, 9.66e8, .132e8),
    )
    wec.morison_element(
        body, point=body.at(0, 0, 10),
        drag_coefficient=(1, 1, 1),
        added_mass_coefficient=(1, 1, 1),
        area=(300, 300, np.pi * 10**2 / 4),
        volume=np.pi * 10**2 * 30,
    )
    result = wec.run(
        PMWave(2, 5, phase_file=phases, directions=(30,), spreading=(1,),
               frequency_range=(.001, 10), water_depth=30,
               current=Current(.8, 45, "power", 30)),
        dt=.01, end_time=20, ramp_time=10, rho=1025,
    )
    _check(result.time, wave[:, 0], 1e-12, "time")
    _check(result.wave_elevation, wave[:, 1], 1e-10, "wave elevation")
    _check(result.body_forces["monopile"], -force[:, 1:7], 1e-3,
           "six-component current Morison force")
    _check(result.body_forces["tower"], np.zeros((2001, 6)), 1e-12,
           "tower force")
    _check(result.bodies["monopile"].position, monopile[:, 1:7], 1e-12,
           "fixed monopile position")
