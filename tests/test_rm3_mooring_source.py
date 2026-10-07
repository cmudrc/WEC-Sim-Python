"""Check the pinned RM3 MooringMatrix source before pairing a Python solver."""

import os
from pathlib import Path

import numpy as np
import pytest


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="paired MATLAB mooring output not provided",
)


def _read(name):
    values = np.loadtxt(Path(REFERENCE) / name, delimiter=",")
    assert np.isfinite(values).all(), name
    return values


def test_mooring_matrix_source_exports_complete_motion_and_force():
    prefix = "RM3_MOORING_MATRIX"
    wave = _read(f"{prefix}_wave.csv")
    mooring = _read(f"{prefix}_mooring1.csv")
    stiffness = _read(f"{prefix}_stiffness.csv")
    assert stiffness.shape == (6, 6)
    np.testing.assert_allclose(stiffness, np.diag([100_000, 0, 0, 0, 0, 0]))
    assert mooring.shape == (40_001, 19)
    np.testing.assert_allclose(mooring[:, 0], np.arange(40_001) * 0.01,
                               rtol=0, atol=1e-8)
    assert wave.shape[1] == 2
    assert wave[0, 0] == 0
    assert wave[-1, 0] >= 400
    assert np.max(np.abs(wave[:, 1])) > 0.01
    assert np.max(np.abs(mooring[:, 13])) > 1_000
    for body in (1, 2):
        values = _read(f"{prefix}_MooringMatrix_body{body}.csv")
        assert values.shape == (40_001, 25)
        np.testing.assert_allclose(values[:, 0], mooring[:, 0],
                                   rtol=0, atol=1e-8)
    pto = _read(f"{prefix}_MooringMatrix_pto1.csv")
    assert pto.shape == (40_001, 25)
    np.testing.assert_allclose(pto[:, 0], mooring[:, 0], rtol=0, atol=1e-8)
