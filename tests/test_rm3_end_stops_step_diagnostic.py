"""Check the finer-step MATLAB End_Stops diagnostic export."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import LinearHardStops


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
MODEL = os.environ.get("WEC_SIM_REFERENCE_MODEL", "RM3_END_STOPS_STEP")
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="finer-step MATLAB End_Stops output not provided",
)


def test_finer_step_export_has_force_and_acceleration_history():
    reference = Path(REFERENCE)
    models = {
        "RM3_END_STOPS_STEP": (0.05, 120, "End_Stops_dt005"),
        "RM3_END_STOPS_STEP_FINE": (0.025, 120, "End_Stops_dt0025"),
        "RM3_END_STOPS_STEP_FINER": (0.0125, 120, "End_Stops_dt00125"),
        "RM3_END_STOPS_FULL_FINE": (0.025, 400, "End_Stops_dt0025_full"),
        "RM3_END_STOPS_FULL_FINER": (0.0125, 400, "End_Stops_dt00125_full"),
    }
    assert MODEL in models
    step, end_time, case = models[MODEL]
    sample_count = round(end_time / step) + 1
    np.testing.assert_array_equal(
        np.loadtxt(reference / f"{MODEL}_run_settings.csv",
                   delimiter=","),
        [step, end_time, 100],
    )
    time = np.arange(sample_count) * step
    pto = np.loadtxt(
        reference / f"{MODEL}_{case}_pto1.csv",
        delimiter=",",
    )
    assert pto.shape == (sample_count, 49)
    np.testing.assert_allclose(pto[:, 0], time, rtol=0, atol=1e-10)
    assert np.isfinite(pto).all()
    stroke, speed, force = pto[:, 3], pto[:, 9], pto[:, 15]
    expected_force = (-1_200_000 * speed
                      + LinearHardStops(-0.6, 0.6, 1e8, 1e8).force(
                          stroke, speed))
    assert np.max(np.abs(force - expected_force)) < 1e-5
    assert np.count_nonzero(np.abs(stroke) > 0.6) > 0
    if MODEL == "RM3_END_STOPS_FULL_FINER":
        # The 400 s trace includes samples inside the source's 0.1 mm
        # transition, so this checks the smoothing law as well as the spring.
        assert np.count_nonzero((np.abs(stroke) > 0.6)
                                & (np.abs(stroke) < 0.6001)) > 0
    for body in (1, 2):
        response = np.loadtxt(
            reference / f"{MODEL}_{case}_body{body}.csv",
            delimiter=",",
        )
        forces = np.loadtxt(
            reference / f"{MODEL}_body{body}_forces.csv",
            delimiter=",",
        )
        assert response.shape == (sample_count, 25)
        assert forces.shape == (sample_count, 25)
        np.testing.assert_allclose(response[:, 0], time, rtol=0, atol=1e-10)
        np.testing.assert_allclose(forces[:, 0], time, rtol=0, atol=1e-10)
        assert np.isfinite(response).all() and np.isfinite(forces).all()
