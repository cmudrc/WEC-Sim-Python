"""Compare a production Python free-decay solver with current MATLAB WEC-Sim."""

import os
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from source.objects.linearHeave import solve_heave_free_decay  # noqa: E402

SPHERE_H5 = os.environ.get("WEC_SIM_SPHERE_H5")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (SPHERE_H5 and REFERENCE), reason="live MATLAB Sphere inputs not provided",
)


@pytest.mark.parametrize("case,initial_displacement", [
    ("0m", 0.0), ("1m", 1.0), ("1m-ME", 1.0), ("3m", 3.0), ("5m", 5.0),
])
def test_sphere_heave_free_decay_against_matlab(case, initial_displacement):
    expected = np.loadtxt(
        Path(REFERENCE) / f"Sphere_{case}_body1.csv", delimiter=",",
    )
    response = solve_heave_free_decay(SPHERE_H5, initial_displacement)
    np.testing.assert_allclose(response.time, expected[:, 0], rtol=0, atol=1e-10)
    scale = max(1.0, abs(initial_displacement))
    assert np.max(np.abs(response.position - expected[:, 3])) < 2e-4 * scale
    assert np.max(np.abs(response.velocity - expected[:, 9])) < 2e-4 * scale
    assert np.max(np.abs(response.force_total - expected[:, 15])) < 160 * scale
