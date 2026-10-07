"""Offline checks for multiple-condition ordering and power accounting."""

import numpy as np
import pytest

from wecsim import MCRCondition, MCRTrace, mcr_grid, mcr_wave_statistics, run_mcr


def test_wave_statistics_selects_positive_cells_in_matlab_order():
    # Periods are columns, heights are rows. Zero cells do not run.
    statistics = np.array([
        [0, 6, 8],
        [1.5, 0.25, 0],
        [2.5, 0, 0.75],
    ])
    selected = mcr_wave_statistics(statistics, [1_200_000, 2_400_000])
    assert selected == (
        MCRCondition(1.5, 6, 1_200_000),
        MCRCondition(2.5, 8, 1_200_000),
        MCRCondition(1.5, 6, 2_400_000),
        MCRCondition(2.5, 8, 2_400_000),
    )
    with pytest.raises(ValueError, match="incomplete"):
        run_mcr(selected, lambda case: MCRTrace(np.array([0, 1]),
                                                  np.array([0, 1])),
                averaging_start_time=1).power_matrix(damping=1_200_000)


def test_power_matrix_axes_and_signed_spring_power():
    conditions = mcr_grid([1.5, 2.5], [6, 8], [1_200_000])
    power = {
        (1.5, 6): np.array([0, -2, 6]),
        (1.5, 8): np.array([0, 4, 8]),
        (2.5, 6): np.array([0, 10, 14]),
        (2.5, 8): np.array([0, 16, 20]),
    }
    result = run_mcr(
        conditions,
        lambda case: MCRTrace(
            np.array([0, 1, 2], dtype=float),
            power[case.wave_height, case.wave_period],
            response={"case": case},
        ),
        averaging_start_time=1,
    )
    np.testing.assert_array_equal(result.mean_absorbed_power, [2, 6, 12, 18])
    assert result.traces[0].response == {"case": conditions[0]}
    matrix = result.power_matrix(damping=1_200_000)
    np.testing.assert_array_equal(matrix.periods, [8, 6])
    np.testing.assert_array_equal(matrix.heights, [1.5, 2.5])
    np.testing.assert_array_equal(matrix.absorbed_power, [[6, 18], [2, 12]])


def test_mcr_rejects_duplicate_grid_and_nonfinite_trace():
    with pytest.raises(ValueError, match="duplicates"):
        mcr_grid([1.5, 1.5], [6], [1_200_000])
    with pytest.raises(ValueError, match="finite absorbed power"):
        run_mcr(
            [MCRCondition(1.5, 6, 1_200_000)],
            lambda case: MCRTrace(np.array([0, 1]), np.array([0, np.nan])),
            averaging_start_time=0,
        )
