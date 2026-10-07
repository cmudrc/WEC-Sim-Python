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


@pytest.mark.parametrize("case,b2b", [
    ("B2B_Case3", False), ("B2B_Case4", True),
    # MATLAB uses a fitted radiation state space for Cases 5 and 6. Python
    # integrates the corresponding impulse-response kernel directly.
    ("B2B_Case5", False), ("B2B_Case6", True),
])
def test_rm3_regular_cic_against_matlab(case, b2b):
    hydro = (Path(APPLICATIONS) / "_Common_Input_Files/RM3/hydroData/rm3.h5").resolve()
    config = {
        "simulation": {"dt": 0.1, "end_time": 400, "ramp_time": 100,
                       "radiation_memory": 60},
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
        for dof, label, limit in ((0, "surge", 0.1), (2, "heave", 0.03),
                                  (4, "pitch", 0.005)):
            position_error = np.max(np.abs(
                response.body_position[:, body - 1, dof] - expected[:, 1 + dof]
            ))
            velocity_error = np.max(np.abs(
                response.body_velocity[:, body - 1, dof] - expected[:, 7 + dof]
            ))
            print(f"{case} body{body} {label}: position {position_error:.6g}, velocity {velocity_error:.6g}")
            assert position_error < limit
            assert velocity_error < limit
    pto = np.loadtxt(reference / f"RM3_B2B_{case}_pto1.csv", delimiter=",")
    np.testing.assert_allclose(response.time, pto[:, 0], rtol=0, atol=1e-10)
    force_error = np.max(np.abs(response.pto_force - pto[:, 15]))
    print(f"{case} PTO force: {force_error:.6g}")
    assert force_error < 30_000
