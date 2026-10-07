"""Keep MATLAB's nonpassive RM3 fit behind an explicit compatibility choice."""

from pathlib import Path

import numpy as np
import pytest

from wecsim.caseDynamics import run_case
from wecsim.rm3Regular import solve_rm3_regular


RM3 = (Path(__file__).parent / "test_objects" / "test_bodyclass"
       / "testData" / "hydroData" / "rm3.h5")


@pytest.mark.parametrize("b2b", [False, True])
def test_negative_surge_fit_is_rejected_by_default(b2b):
    with pytest.raises(ValueError, match="negative zero-frequency common-surge damping"):
        solve_rm3_regular(
            RM3, b2b=b2b, state_space=True, wave_height=0, end_time=0,
        )


@pytest.mark.parametrize("b2b", [False, True])
def test_matlab_compatibility_requires_explicit_override(b2b):
    response = solve_rm3_regular(
        RM3, b2b=b2b, state_space=True,
        allow_negative_surge_damping=True,
        wave_height=0, end_time=0,
    )
    assert np.isfinite(response.body_position).all()
    assert np.isfinite(response.body_velocity).all()


def test_override_cannot_be_silently_ignored():
    with pytest.raises(ValueError, match="requires state_space"):
        solve_rm3_regular(
            RM3, allow_negative_surge_damping=True,
            wave_height=0, end_time=0,
        )


def test_case_api_requires_explicit_compatibility_choice():
    case = {
        "simulation": {"dt": 0.1, "end_time": 0, "state_space": True},
        "wave": {"type": "regularCIC", "height": 0, "period": 8},
        "bodies": [
            {"hydro_file": str(RM3.resolve()), "hydro_body": 1,
             "mass": "equilibrium", "pitch_inertia": 21_306_090.66},
            {"hydro_file": str(RM3.resolve()), "hydro_body": 2,
             "mass": "equilibrium", "pitch_inertia": 94_407_091.24},
        ],
        "constraint": {"kind": "floating_joint", "location": [0, 0, 0]},
        "pto": {"kind": "relative_heave", "damping": 1_200_000},
    }
    with pytest.raises(ValueError, match="negative zero-frequency common-surge damping"):
        run_case(case)
    case["constraint"] = {"kind": "linear_subspace", "coordinates": []}
    with pytest.raises(ValueError, match="fitted radiation settings require floating_joint"):
        run_case(case)
