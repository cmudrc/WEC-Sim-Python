"""Compare the eight published RM3 MCR Option 1 parameter combinations."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim.caseDynamics import run_case


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="paired MATLAB RM3 MCR output not provided",
)


@pytest.mark.parametrize("damping", [1_200_000, 2_400_000])
@pytest.mark.parametrize("period", [6, 8])
@pytest.mark.parametrize("height", [1.5, 2.5])
def test_rm3_mcr_condition_against_matlab(height, period, damping):
    label = f"H{round(10 * height)}_T{period}_D{round(damping / 100_000)}"
    hydro = (Path(APPLICATIONS) / "_Common_Input_Files/RM3/hydroData/rm3.h5").resolve()
    case = {
        "simulation": {"dt": 0.1, "end_time": 400, "ramp_time": 100,
                       "radiation_memory": 60},
        "wave": {"type": "regularCIC", "height": height, "period": period},
        "bodies": [
            {"hydro_file": str(hydro), "hydro_body": 1,
             "mass": "equilibrium", "pitch_inertia": 21_306_090.66},
            {"hydro_file": str(hydro), "hydro_body": 2,
             "mass": "equilibrium", "pitch_inertia": 94_407_091.24},
        ],
        "constraint": {"kind": "floating_joint", "location": [0, 0, 0]},
        "pto": {"kind": "relative_heave", "damping": damping},
    }
    response = run_case(case)
    reference = Path(REFERENCE)
    for body in (1, 2):
        expected = np.loadtxt(reference / f"RM3_MCR_{label}_body{body}.csv",
                              delimiter=",")
        np.testing.assert_allclose(response.time, expected[:, 0], rtol=0, atol=1e-10)
        for dof, name, limit in ((0, "surge", 0.2), (2, "heave", 0.03),
                                 (4, "pitch", 0.01)):
            position_error = np.max(np.abs(
                response.body_position[:, body - 1, dof] - expected[:, 1 + dof]
            ))
            velocity_error = np.max(np.abs(
                response.body_velocity[:, body - 1, dof] - expected[:, 7 + dof]
            ))
            print(f"{label} body{body} {name}: position {position_error:.6g}, velocity {velocity_error:.6g}")
            assert position_error < limit
            assert velocity_error < limit
    pto = np.loadtxt(reference / f"RM3_MCR_{label}_pto1.csv", delimiter=",")
    np.testing.assert_allclose(response.time, pto[:, 0], rtol=0, atol=1e-10)
    force_error = np.max(np.abs(response.pto_force - pto[:, 15]))
    print(f"{label} PTO force: {force_error:.6g}")
    assert force_error < 50_000
