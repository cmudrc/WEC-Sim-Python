"""Independent coupled MBARI Cable motion against pinned MATLAB WEC-Sim."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import run_mbari_cable


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
HYDRO = os.environ.get("WEC_SIM_CABLE_H5")
pytestmark = pytest.mark.skipif(
    not REFERENCE or not HYDRO,
    reason="pinned MBARI Cable source output and hydrodynamics absent",
)


def _source(name):
    return np.loadtxt(Path(REFERENCE) / f"CABLE_SOURCE_{name}.csv",
                      delimiter=",", ndmin=2)


def test_published_mbari_three_body_coupled_motion_and_cable():
    source_body = np.stack([_source(f"Cable_body{i}") for i in (1, 2, 3)],
                           axis=1)
    source_cable = _source("cable1")
    source_wave = _source("wave")
    assert source_body.shape == (12_001, 3, 25)
    assert source_cable.shape == (12_001, 37)
    assert source_wave.shape == (12_001, 2)

    result = run_mbari_cable(HYDRO)
    np.testing.assert_allclose(result.time, source_wave[:, 0],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.wave_elevation, source_wave[:, 1],
                               rtol=0, atol=1e-11)
    position_error = np.max(np.abs(
        result.body_position - source_body[:, :, 1:7]), axis=0)
    velocity_error = np.max(np.abs(
        result.body_velocity - source_body[:, :, 7:13]), axis=0)
    position_gate = np.array([
        [.01, 1e-10, .05, 1e-10, .01, 1e-10],
        [.005, 1e-10, .05, 1e-10, .003, 1e-10],
        [.01, 1e-10, .05, 1e-10, .003, 1e-10],
    ])
    velocity_gate = np.array([
        [.03, 1e-10, .45, 1e-10, .05, 1e-10],
        [.015, 1e-10, .35, 1e-10, .01, 1e-10],
        [.04, 1e-10, .45, 1e-10, .02, 1e-10],
    ])
    np.testing.assert_array_less(position_error, position_gate)
    np.testing.assert_array_less(velocity_error, velocity_gate)

    displacement_error = result.cable_displacement - source_cable[:, 3]
    speed_error = result.cable_speed - source_cable[:, 9]
    force_error = result.cable_force_z - source_cable[:, 21]
    assert np.max(np.abs(displacement_error)) < .06
    assert np.max(np.abs(speed_error)) < .8
    assert np.max(np.abs(force_error)) < 45_000
    assert np.sqrt(np.mean(force_error**2)) < 3_500
    assert abs(np.trapezoid(force_error, result.time)) < 500
    assert abs(np.mean(result.cable_force_z != 0)
               - np.mean(source_cable[:, 21] != 0)) < .02
