"""Compare published RM3 regularCIC body-to-body cases with MATLAB."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim.caseDynamics import run_case


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="paired MATLAB RM3 regularCIC output not provided",
)


def _max_error(actual, expected, limit, label):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, label
    assert np.isfinite(expected).all(), label
    error = np.max(np.abs(actual - expected))
    assert error < limit, f"{label}: max error {error:.6g} exceeds {limit}"


@pytest.mark.parametrize("case,b2b", [
    ("B2B_Case3", False), ("B2B_Case4", True),
])
def test_rm3_regular_cic_against_matlab(case, b2b):
    hydro = (Path(APPLICATIONS) / "_Common_Input_Files/RM3/hydroData/rm3.h5").resolve()
    config = {
        "simulation": {"dt": 0.1, "end_time": 400, "ramp_time": 100,
                       "radiation_memory": 60,
                       "added_mass_scheme": "simulink_delay"},
        "wave": {"type": "regularCIC", "height": 2.5, "period": 8},
        "bodies": [
            {"hydro_file": str(hydro), "hydro_body": 1,
             "mass": "equilibrium", "pitch_inertia": 21_306_090.66},
            {"hydro_file": str(hydro), "hydro_body": 2,
             "mass": "equilibrium", "pitch_inertia": 94_407_091.24},
        ],
        "constraint": {"kind": "floating_joint", "location": [0, 0, 0]},
        "pto": {"kind": "relative_heave", "damping": 1_200_000},
        "body_to_body": b2b,
    }
    response = run_case(config)
    reference = Path(REFERENCE)
    for body in (1, 2):
        expected = np.loadtxt(
            reference / f"RM3_B2B_{case}_body{body}.csv", delimiter=",",
        )
        np.testing.assert_allclose(response.time, expected[:, 0], rtol=0, atol=1e-10)
        for dof, label, position_limit, velocity_limit in (
            (0, "surge", 0.075, 0.04),
            (2, "heave", 0.005, 0.0035),
            (4, "pitch", 0.0015, 0.0011),
        ):
            _max_error(response.body_position[:, body - 1, dof],
                       expected[:, 1 + dof], position_limit,
                       f"{case} body{body} {label} position")
            _max_error(response.body_velocity[:, body - 1, dof],
                       expected[:, 7 + dof], velocity_limit,
                       f"{case} body{body} {label} velocity")
    pto = np.loadtxt(reference / f"RM3_B2B_{case}_pto1.csv", delimiter=",")
    np.testing.assert_allclose(response.time, pto[:, 0], rtol=0, atol=1e-10)
    # Recover the PTO translation at the floating joint, accounting for
    # the two body-center offsets when the assembly pitches.
    center_gap = (response.body_position[0, 0, 2]
                  - response.body_position[0, 1, 2])
    pitch = response.body_position[:, 0, 4]
    pitch_speed = response.body_velocity[:, 0, 4]
    cosine = np.cos(pitch)
    stroke = ((response.body_position[:, 0, 2]
               - response.body_position[:, 1, 2]) / cosine - center_gap)
    speed = ((response.body_velocity[:, 0, 2]
              - response.body_velocity[:, 1, 2]
              + (center_gap + stroke) * np.sin(pitch) * pitch_speed) / cosine)
    _max_error(stroke, pto[:, 3], 0.008, f"{case} PTO stroke")
    _max_error(speed, pto[:, 9], 0.005, f"{case} PTO speed")
    _max_error(response.pto_force, pto[:, 15], 6_000, f"{case} PTO force")
    _max_error(response.pto_force * speed, pto[:, 21], 7_000,
               f"{case} PTO power")
