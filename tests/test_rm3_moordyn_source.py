"""Validate saved output from the pinned MATLAB RM3 MoorDyn application."""

import os
from pathlib import Path

import numpy as np
import pytest


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="MATLAB MoorDyn baseline output not provided",
)


def _read(name):
    values = np.loadtxt(Path(REFERENCE) / name, delimiter=",")
    assert values.ndim == 2 and np.isfinite(values).all(), name
    return values


def test_published_moordyn_case_exports_motion_and_line_tension():
    prefix = "RM3_MOORDYN"
    wave = _read(f"{prefix}_wave.csv")
    coupling = _read(f"{prefix}_coupling.csv")
    fairlead = _read(f"{prefix}_fairlead_tension.csv")
    assert wave.shape[1] == 2
    assert coupling.shape == (40_001, 19)
    np.testing.assert_allclose(coupling[:, 0], np.arange(40_001) * 0.01,
                               rtol=0, atol=1e-8)
    assert np.max(np.abs(coupling[:, 13:19])) > 1_000
    assert fairlead.shape[1] == 4
    assert fairlead.shape[0] > 100
    assert np.all(np.diff(fairlead[:, 0]) > 0)
    assert fairlead[0, 0] <= 0.01
    assert fairlead[-1, 0] >= 390
    assert np.max(np.abs(fairlead[:, 1:])) > 1_000
    assert np.max(np.abs(wave[:, 1])) > 0.01

    for body in (1, 2):
        motion = _read(f"{prefix}_MoorDyn_body{body}.csv")
        assert motion.shape == (40_001, 25)
        np.testing.assert_allclose(motion[:, 0], np.arange(40_001) * 0.01,
                                   rtol=0, atol=1e-8)
    pto = _read(f"{prefix}_MoorDyn_pto1.csv")
    assert pto.shape == (40_001, 25)
    np.testing.assert_allclose(pto[:, 0], np.arange(40_001) * 0.01,
                               rtol=0, atol=1e-8)
