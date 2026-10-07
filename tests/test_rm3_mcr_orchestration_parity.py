"""Pair RM3 MCR inputs, eight trajectories, and power matrices with MATLAB."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import mcr_grid, mcr_mat_file, mcr_wave_statistics, run_rm3_mcr


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
LABEL = os.environ.get("WEC_SIM_MCR_LABEL", "RM3_MCR_MAT")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="paired MATLAB RM3 MCR output not provided",
)


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: max error {error:.6g} exceeds {limit}"


def test_published_mcr_inputs_and_actual_matlab_mcr_outputs():
    root = Path(APPLICATIONS)
    reference = Path(REFERENCE)
    mcr_root = root / "Multiple_Condition_Runs"
    option_1 = mcr_grid([1.5, 2.5], [6, 8],
                        [1_200_000, 2_400_000])
    option_2 = mcr_wave_statistics(
        mcr_root / "RM3_MCROPT2/waveStatistic.xlsx",
        [1_200_000, 2_400_000],
    )
    option_3 = mcr_mat_file(mcr_root / "RM3_MCROPT3/mcrExample.mat")
    assert option_1 == option_2 == option_3

    assert LABEL in ("RM3_MCR_ARRAY", "RM3_MCR_EXCEL", "RM3_MCR_MAT")
    summary = np.loadtxt(reference / f"{LABEL}_summary.csv", delimiter=",")
    assert summary.shape == (8, 6)
    np.testing.assert_array_equal(
        np.array([tuple(vars(case).values()) for case in option_3]),
        summary[:, :4],
    )
    np.testing.assert_array_equal(summary[:, 2], summary[:, 5])

    conditions = {
        "RM3_MCR_ARRAY": option_1,
        "RM3_MCR_EXCEL": option_2,
        "RM3_MCR_MAT": option_3,
    }[LABEL]
    hydro = root / "_Common_Input_Files/RM3/hydroData/rm3.h5"
    result = run_rm3_mcr(hydro, conditions)
    assert result.conditions == conditions
    assert len(result.traces) == 8

    for index, trace in enumerate(result.traces, start=1):
        response = trace.response
        assert response is not None
        assert response.body_position.shape == (4001, 2, 6)
        np.testing.assert_allclose(trace.time, response.time, rtol=0, atol=0)
        for body in (1, 2):
            expected = np.loadtxt(
                reference / f"{LABEL}_case{index}_body{body}.csv",
                delimiter=",",
            )
            assert expected.shape == (4001, 25)
            np.testing.assert_allclose(response.time, expected[:, 0],
                                       rtol=0, atol=1e-10)
            for dof, axis, position_limit, velocity_limit in (
                (0, "surge", 0.115, 0.050),
                (2, "heave", 0.007, 0.007),
                (4, "pitch", 0.0045, 0.0018),
            ):
                _max_error(response.body_position[:, body - 1, dof],
                           expected[:, 1 + dof], position_limit,
                           f"{LABEL} case {index} body {body} {axis} position")
                _max_error(response.body_velocity[:, body - 1, dof],
                           expected[:, 7 + dof], velocity_limit,
                           f"{LABEL} case {index} body {body} {axis} velocity")

        pto = np.loadtxt(
            reference / f"{LABEL}_case{index}_pto1.csv", delimiter=",",
        )
        assert pto.shape == (4001, 25)
        np.testing.assert_allclose(response.time, pto[:, 0], rtol=0, atol=1e-10)
        center_gap = (response.body_position[0, 0, 2]
                      - response.body_position[0, 1, 2])
        pitch = response.body_position[:, 0, 4]
        stroke = (response.body_position[:, 0, 2]
                  - response.body_position[:, 1, 2]
                  - center_gap * np.cos(pitch))
        _max_error(stroke, pto[:, 3], 0.0115,
                   f"MCR case {index} PTO stroke")
        _max_error(response.pto_velocity, pto[:, 9], 0.0075,
                   f"MCR case {index} PTO speed")
        _max_error(response.pto_force, pto[:, 15], 9_000,
                   f"MCR case {index} PTO force")
        _max_error(trace.absorbed_power, -pto[:, 21], 7_500,
                   f"MCR case {index} PTO absorbed power")
        # MATLAB's userDefinedFunctionsMCR averages rows 2000:end, including
        # the sample at t=199.9 s. Its absorbed-power sign is negative.
        np.testing.assert_allclose(-summary[index - 1, 4],
                                   np.mean(-pto[1999:, 21]), rtol=0, atol=1e-6)

    expected_power = -summary[:, 4]
    _max_error(result.mean_absorbed_power, expected_power, 2_100,
               "eight MATLAB MCR mean powers")
    assert np.max(np.abs(result.mean_absorbed_power - expected_power)
                  / expected_power) < 0.008
    for damping, offset in ((1_200_000, 0), (2_400_000, 4)):
        matrix = result.power_matrix(damping=damping)
        np.testing.assert_array_equal(matrix.periods, [8, 6])
        np.testing.assert_array_equal(matrix.heights, [1.5, 2.5])
        label = "D12" if damping == 1_200_000 else "D24"
        matlab_matrix = np.loadtxt(
            reference / f"{LABEL}_power_matrix_{label}.csv",
            delimiter=",",
        )
        assert matlab_matrix.shape == (2, 2)
        np.testing.assert_allclose(
            matlab_matrix,
            summary[np.array([[1, 3], [0, 2]]) + offset, 4],
            rtol=0, atol=1e-6,
        )
        _max_error(matrix.absorbed_power, -matlab_matrix, 2_100,
                   f"MCR {damping:g} PTO power matrix")
