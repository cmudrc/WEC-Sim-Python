"""Gate full-duration source time-step convergence before claiming parity."""

import os
from pathlib import Path

import numpy as np
import pytest


FINE = os.environ.get("WEC_SIM_END_STOPS_FINE_OUTPUT_DIR")
FINER = os.environ.get("WEC_SIM_END_STOPS_FINER_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (FINE and FINER), reason="paired full-duration MATLAB outputs not provided",
)


def _load(directory, model, case, channel):
    return np.loadtxt(Path(directory) / f"{model}_{case}_{channel}.csv",
                      delimiter=",")


def _max_error(left, right, limit, label):
    assert left.shape == right.shape, label
    assert np.isfinite(left).all() and np.isfinite(right).all(), label
    error = np.max(np.abs(left - right))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


def _added_mass_recurrence_residual(directory, model, body, dt):
    forces = np.loadtxt(
        Path(directory) / f"{model}_body{body}_forces.csv", delimiter=",",
    )
    applied = np.loadtxt(
        Path(directory) / f"{model}_body{body}_added_mass_applied.csv",
        delimiter=",",
    )
    acceleration = forces[:, 19:25]
    delayed = (acceleration[1:-1]
               + (1 - 1e-7 / dt) *
               (acceleration[1:-1] - acceleration[:-2]))
    predicted = delayed @ applied.T
    observed = forces[2:, 7:13]
    return np.max(np.abs(predicted[:, [0, 2]]
                         - observed[:, [0, 2]]))


def test_refined_matlab_end_stops_converge_through_400_seconds():
    fine = "RM3_END_STOPS_FULL_FINE"
    finer = "RM3_END_STOPS_FULL_FINER"
    p1 = _load(FINE, fine, "End_Stops_dt0025_full", "pto1")
    p2 = _load(FINER, finer, "End_Stops_dt00125_full", "pto1")
    assert p1.shape == (16001, 49)
    assert p2.shape == (32001, 49)
    np.testing.assert_allclose(p1[:, 0], np.arange(16001) * 0.025,
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(p2[:, 0], np.arange(32001) * 0.0125,
                               rtol=0, atol=1e-10)
    _max_error(p1[:, 3], p2[::2, 3], 0.004, "PTO stroke")
    _max_error(p1[:, 9], p2[::2, 9], 0.025, "PTO speed")
    _max_error(p1[:, 15], p2[::2, 15], 350_000, "PTO force")
    e1 = np.trapezoid(1_200_000 * p1[:, 9]**2, p1[:, 0])
    e2 = np.trapezoid(1_200_000 * p2[:, 9]**2, p2[:, 0])
    assert abs(e1 - e2) / e2 < 0.006

    for body in (1, 2):
        b1 = _load(FINE, fine, "End_Stops_dt0025_full", f"body{body}")
        b2 = _load(FINER, finer, "End_Stops_dt00125_full", f"body{body}")
        assert b1.shape == (16001, 25)
        assert b2.shape == (32001, 25)
        for dof, name, position_limit, speed_limit in (
            (0, "surge", 0.0025, 0.004),
            (2, "heave", 0.004, 0.022),
            (4, "pitch", 8e-5, 1e-4),
        ):
            _max_error(b1[:, 1 + dof], b2[::2, 1 + dof],
                       position_limit, f"body {body} {name} position")
            _max_error(b1[:, 7 + dof], b2[::2, 7 + dof],
                       speed_limit, f"body {body} {name} speed")
        for directory, model, dt in ((FINE, fine, 0.025),
                                     (FINER, finer, 0.0125)):
            assert _added_mass_recurrence_residual(
                directory, model, body, dt,
            ) < 1e-3
