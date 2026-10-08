"""Pair the public Python floating-joint builder with published RM3 MCR."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularCICWave, WEC


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="paired MATLAB RM3 MCR output not provided",
)


def _max_error(actual, expected, limit, label):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


def test_public_floating_joint_tracks_published_mcr_condition():
    root = Path(APPLICATIONS)
    reference = Path(REFERENCE)
    hydro = "_Common_Input_Files/RM3/hydroData/rm3.h5"
    wec = WEC("RM3 MCR condition 1")
    float_body = wec.body("float", hydro,
                          inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", hydro,
                    inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(
        float_body, spar, damping=1_200_000,
        radiation_method="convolution",
        added_mass_scheme="simulink_delay",
    )
    result = wec.run(
        RegularCICWave(1.5, 6), dt=0.1, end_time=400,
        ramp_time=100, radiation_memory=60, base_dir=root,
    )
    assert result.case["constraint"]["kind"] == "floating_joint"
    assert result.time.shape == (4001,)
    for body_index, name in enumerate(("float", "spar"), start=1):
        expected = np.loadtxt(
            reference / f"RM3_MCR_ARRAY_case1_body{body_index}.csv",
            delimiter=",",
        )
        _max_error(result.time, expected[:, 0], 1e-10, "time")
        for dof, position_limit, speed_limit in (
            (0, 0.06, 0.004),
            (2, 0.002, 0.001),
            (4, 0.0001, 0.00005),
        ):
            _max_error(result.bodies[name].position[:, dof],
                       expected[:, 1 + dof], position_limit,
                       f"{name} position DOF {dof}")
            _max_error(result.bodies[name].velocity[:, dof],
                       expected[:, 7 + dof], speed_limit,
                       f"{name} velocity DOF {dof}")
    expected_pto = np.loadtxt(
        reference / "RM3_MCR_ARRAY_case1_pto1.csv", delimiter=",",
    )
    pto = result.ptos["relative_heave"]
    _max_error(pto.stroke, expected_pto[:, 3], 0.002, "PTO stroke")
    _max_error(pto.velocity, expected_pto[:, 9], 0.001, "PTO speed")
    _max_error(pto.force, expected_pto[:, 15], 500, "PTO force")
    _max_error(pto.absorbed_power, -expected_pto[:, 21], 500,
               "PTO absorbed power")
    np.testing.assert_allclose(
        pto.stroke,
        result.coordinates["float_heave"].position
        - result.coordinates["spar_heave"].position,
        rtol=0, atol=1e-12,
    )
