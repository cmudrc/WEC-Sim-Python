"""Independent WaveBot motion against the pinned 600 s MATLAB application."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import run_wavebot_impedance


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
HYDRO = os.environ.get("WEC_SIM_WAVEBOT_H5")
MULTISINE = os.environ.get("WEC_SIM_WAVEBOT_MULTISINE")
pytestmark = pytest.mark.skipif(
    not (REFERENCE and HYDRO and MULTISINE),
    reason="pinned WaveBot output, hydrodynamics, and multisine are absent",
)


def _source(name):
    return np.loadtxt(
        Path(REFERENCE) / f"WAVEBOT_IMPEDANCE_{name}.csv",
        delimiter=",", ndmin=2,
    )


@pytest.fixture(scope="module")
def comparison():
    return (
        _source("CalcImpedance_body1"),
        _source("body1_forces"),
        _source("mooring"),
        _source("actuation"),
        run_wavebot_impedance(HYDRO, MULTISINE,
                              source_linear_damping=True),
    )


def test_published_wavebot_coupled_motion_and_force_channels(comparison):
    body, force, mooring, actuation, result = comparison
    assert body.shape == (60_001, 25)
    np.testing.assert_allclose(result.time, body[:, 0], rtol=0, atol=1e-8)

    active = [0, 2, 4]
    position_error = np.max(
        np.abs(result.body_position[:, active] - body[:, [1, 3, 5]]), axis=0,
    )
    velocity_error = np.max(
        np.abs(result.body_velocity[:, active] - body[:, [7, 9, 11]]), axis=0,
    )
    np.testing.assert_array_less(position_error, [.0002, .0003, .0004])
    np.testing.assert_array_less(velocity_error, [.0012, .001, .0018])

    sampled_command = np.column_stack([
        np.interp(result.time, actuation[:, 0], actuation[:, i])
        for i in (1, 2, 3)
    ])
    sampled_command[:, 2] *= -1
    np.testing.assert_allclose(result.applied_force, sampled_command,
                               rtol=0, atol=1e-8)
    force_pairs = (
        (result.mooring_force, mooring[:, 13:19], [4, 1.3, .4]),
        (result.linear_damping_force, -force[:, 25:31], [1.3, 1.3, .15]),
        (result.drag_force, -force[:, 19:25], [.2, .8, .5]),
    )
    for python_force, matlab_force, limits in force_pairs:
        error = np.max(np.abs(python_force[:, active]
                              - matlab_force[:, active]), axis=0)
        np.testing.assert_array_less(error, limits)


def test_ordinary_diagonal_damping_stays_the_default(comparison):
    body, _, _, _, source_compatible = comparison
    physical = run_wavebot_impedance(HYDRO, MULTISINE)
    # Each positive diagonal coefficient dissipates power independently.
    damping_power = np.sum(
        physical.linear_damping_force * physical.body_velocity, axis=1,
    )
    assert np.max(damping_power) < 1e-10
    # MATLAB's linearly indexed matrix is different from that physical
    # default. Keep the source convention isolated and report both traces.
    assert np.max(np.abs(physical.body_position[:, 2] - body[:, 3])) > .01
    assert np.max(np.abs(physical.body_position[:, 4] - body[:, 5])) > .005
    assert np.max(np.abs(source_compatible.body_position[:, 2]
                         - body[:, 3])) < .0003
