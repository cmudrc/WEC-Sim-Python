"""Pair the published floating OWC's seven-coordinate column joint."""

import os
from pathlib import Path

import h5py
import numpy as np
import pytest

from wecsim import FloatingOwcColumnJoint


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
HYDRO = os.environ.get("WEC_SIM_FLOATING_OWC_H5")


def _read(name):
    return np.loadtxt(Path(REFERENCE) / f"FLOATING_OWC_SOURCE_{name}.csv",
                      delimiter=",", ndmin=2)


def test_column_joint_rotates_its_offset_and_sliding_speed():
    joint = FloatingOwcColumnJoint(10)
    floater = np.array([0, 0, -30, 0, np.pi / 2, 0])
    velocity = np.array([0, 0, 0, 0, 1, 0])
    np.testing.assert_allclose(
        joint.column_pose(floater, 2)[:3], [12, 0, -30], atol=1e-14,
    )
    np.testing.assert_allclose(
        joint.column_velocity(floater, velocity, 2, 3)[:3],
        [3, 0, -12], atol=1e-14,
    )
    with pytest.raises(ValueError, match="separation"):
        FloatingOwcColumnJoint(0)


@pytest.mark.skipif(not (REFERENCE and HYDRO),
                    reason="pinned floating OWC bodies, PTO, and HDF5 absent")
def test_published_column_joint_reconstructs_body_and_chamber_motion():
    floater = _read("FloatingOWC_body1")
    column = _read("FloatingOWC_body2")
    pto = _read("FloatingOWC_pto1")
    with h5py.File(HYDRO) as hydro:
        floater_z = float(np.asarray(hydro["body1/properties/cg"])[2, 0])
        column_z = float(np.asarray(hydro["body2/properties/cg"])[2, 0])
    joint = FloatingOwcColumnJoint(column_z - floater_z)

    np.testing.assert_allclose(pto[:, 0], floater[:, 0], rtol=0, atol=1e-10)
    predicted_pose = joint.column_pose(floater[:, 1:7], pto[:, 3])
    predicted_velocity = joint.column_velocity(
        floater[:, 1:7], floater[:, 7:13], pto[:, 3], pto[:, 9],
    )
    assert floater.shape[0] == column.shape[0] == pto.shape[0] == 50001
    assert np.max(np.abs(predicted_pose - column[:, 1:7])) < 1e-10
    assert np.max(np.abs(predicted_velocity - column[:, 7:13])) < 1e-10

    chamber_heave = _read("x_xOWC")
    chamber_speed = _read("x_signal4")
    np.testing.assert_allclose(chamber_heave[:, 0], floater[:, 0],
                               rtol=0, atol=1e-10)
    assert np.max(np.abs(predicted_pose[:, 2] - column_z
                         - chamber_heave[:, 1])) < 1e-10
    assert np.max(np.abs(predicted_velocity[:, 2]
                         - chamber_speed[:, 1])) < 1e-10
