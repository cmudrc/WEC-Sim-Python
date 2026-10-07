"""Check the finer-step MATLAB End_Stops diagnostic export."""

import os
from pathlib import Path

import numpy as np
import pytest


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="finer-step MATLAB End_Stops output not provided",
)


def test_finer_step_export_has_force_and_acceleration_history():
    reference = Path(REFERENCE)
    np.testing.assert_array_equal(
        np.loadtxt(reference / "RM3_END_STOPS_STEP_run_settings.csv",
                   delimiter=","),
        [0.05, 120, 100],
    )
    time = np.arange(2401) * 0.05
    pto = np.loadtxt(
        reference / "RM3_END_STOPS_STEP_End_Stops_dt005_pto1.csv",
        delimiter=",",
    )
    assert pto.shape == (2401, 49)
    np.testing.assert_allclose(pto[:, 0], time, rtol=0, atol=1e-10)
    assert np.isfinite(pto).all()
    stroke, speed, force = pto[:, 3], pto[:, 9], pto[:, 15]
    expected_force = (-1_200_000 * speed
                      - 1e8 * np.maximum(stroke - 0.6, 0)
                      + 1e8 * np.maximum(-0.6 - stroke, 0))
    assert np.max(np.abs(force - expected_force)) < 1e-5
    assert np.count_nonzero(np.abs(stroke) > 0.6) > 0
    for body in (1, 2):
        response = np.loadtxt(
            reference / f"RM3_END_STOPS_STEP_End_Stops_dt005_body{body}.csv",
            delimiter=",",
        )
        forces = np.loadtxt(
            reference / f"RM3_END_STOPS_STEP_body{body}_forces.csv",
            delimiter=",",
        )
        assert response.shape == (2401, 25)
        assert forces.shape == (2401, 25)
        np.testing.assert_allclose(response[:, 0], time, rtol=0, atol=1e-10)
        np.testing.assert_allclose(forces[:, 0], time, rtol=0, atol=1e-10)
        assert np.isfinite(response).all() and np.isfinite(forces).all()
