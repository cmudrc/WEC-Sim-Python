"""Pinned RM3 hydraulic PTO source signals and replayed incident sea.

This is an input/force diagnostic. It does not claim coupled hydraulic
trajectory parity with a Python PTO model.
"""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim.irregularWave import (
    pm_equal_energy_components, synthesize_irregular_response,
)


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
RM3_H5 = os.environ.get("WEC_SIM_RM3_H5")
pytestmark = pytest.mark.skipif(
    not (REFERENCE and RM3_H5),
    reason="pinned RM3 hydraulic MATLAB output and HDF5 absent",
)


def _load(name):
    values = np.loadtxt(Path(REFERENCE) / f"RM3_HYDRAULIC_SOURCE_{name}.csv",
                        delimiter=",")
    assert values.ndim == 2 and np.isfinite(values).all(), name
    return values


def test_hydraulic_source_wave_and_component_traces():
    components = _load("components")
    matlab_wave = _load("wave")
    assert components.shape == (500, 4)
    assert matlab_wave.shape == (40_001, 2)
    generated = pm_equal_energy_components(
        RM3_H5, significant_height=2.5, peak_period=8,
        directions=(0,), spreading=(1,), phase=components[:, 3:4],
    )
    for actual, saved in ((generated.omega, components[:, 0]),
                          (generated.spectral_amplitude, components[:, 1]),
                          (generated.d_omega, components[:, 2])):
        assert np.max(np.abs(actual - saved)) < 1e-12
    replay = synthesize_irregular_response(
        RM3_H5, generated, dt=.01, end_time=400, ramp_time=100,
    )
    np.testing.assert_allclose(replay.time, matlab_wave[:, 0], atol=1e-10,
                               rtol=0)
    assert np.max(np.abs(replay.elevation - matlab_wave[:, 1])) < 1e-10

    for name, columns in (("cylinder", 4), ("valve", 5),
                          ("high_accumulator", 3),
                          ("low_accumulator", 3), ("motor", 5),
                          ("generator", 5)):
        values = _load(name)
        assert values.shape == (40_001, columns), name
        np.testing.assert_allclose(values[:, 0], matlab_wave[:, 0],
                                   atol=1e-10, rtol=0)
    cylinder = _load("cylinder")
    assert np.max(np.abs(cylinder[:, 2])) > 1e3
