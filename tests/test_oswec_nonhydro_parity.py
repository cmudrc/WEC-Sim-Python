"""Pair the published OSWEC nonhydrodynamic-base case with MATLAB output."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim.caseDynamics import run_case


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="paired MATLAB OSWEC nonhydrodynamic output not provided",
)


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(expected).all(), label
    error = np.max(np.abs(actual - expected))
    assert error < limit, f"{label}: max error {error:.6g} exceeds {limit}"


def test_published_nonhydrodynamic_base_against_matlab():
    hydro = (Path(APPLICATIONS)
             / "_Common_Input_Files/OSWEC/hydroData/oswec.h5").resolve()
    case = {
        "simulation": {"dt": 0.1, "end_time": 400, "ramp_time": 100},
        "wave": {"type": "regular", "height": 2.5, "period": 8},
        "bodies": [
            {"name": "flap", "hydro_file": str(hydro), "hydro_body": 1,
             "mass": 127_000, "pitch_inertia": 1.85e6},
            {"name": "base", "nonhydro": True, "fixed": True,
             "center_gravity": [0, 0, -10.9], "mass": 999,
             "inertia": [1, 1, 1]},
        ],
        "constraint": {"kind": "fixed_hinge", "location": [0, 0, -8.9]},
        "pto": {"kind": "pitch", "damping": 0},
    }
    response = run_case(case)
    reference = Path(REFERENCE)
    flap = np.loadtxt(
        reference / "OSWEC_Nonhydro_Nonhydro_body1.csv", delimiter=",",
    )
    base = np.loadtxt(
        reference / "OSWEC_Nonhydro_Nonhydro_body2.csv", delimiter=",",
    )
    np.testing.assert_allclose(response.time, flap[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(response.time, base[:, 0], rtol=0, atol=1e-10)
    for dof, label, position_limit, velocity_limit in (
        (0, "surge", 0.05, 0.05),
        (2, "heave", 0.02, 0.025),
        (4, "pitch", 0.01, 0.01),
    ):
        _max_error(response.body_position[:, 0, dof], flap[:, 1 + dof],
                   position_limit, f"flap {label} position")
        _max_error(response.body_velocity[:, 0, dof], flap[:, 7 + dof],
                   velocity_limit, f"flap {label} velocity")
    np.testing.assert_allclose(
        response.body_position[:, 1, :], base[:, 1:7], rtol=0, atol=1e-10,
    )
    np.testing.assert_allclose(
        response.body_velocity[:, 1, :], base[:, 7:13], rtol=0, atol=1e-10,
    )
    pto = np.loadtxt(
        reference / "OSWEC_Nonhydro_Nonhydro_pto1.csv", delimiter=",",
    )
    _max_error(response.body_position[:, 0, 4], pto[:, 5], 0.01, "PTO angle")
    _max_error(response.body_velocity[:, 0, 4], pto[:, 11], 0.01,
               "PTO angular speed")
    _max_error(response.pto_force, pto[:, 17], 1e-8, "PTO torque")
