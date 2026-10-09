"""MOST three-line restoring force against the pinned MATLAB source."""

import os

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import MostStaticMooring


def test_most_static_mooring_symmetry_and_pose_validation():
    mooring = MostStaticMooring()
    force, tensions = mooring.force(np.zeros(6))
    np.testing.assert_allclose(force[[0, 1, 3, 4, 5]], 0, atol=1e-6)
    np.testing.assert_allclose(tensions, np.broadcast_to(tensions[0], (3, 2)),
                               rtol=0, atol=1e-6)
    assert force[2] < 0
    with pytest.raises(ValueError, match="six finite"):
        mooring.force([0, 0, np.nan, 0, 0, 0])


@pytest.mark.skipif(not os.environ.get("WEC_SIM_MOST_MOORING_BASELINE"),
                    reason="pinned MATLAB MOST mooring output not provided")
def test_most_static_mooring_against_matlab_source():
    reference = loadmat(os.environ["WEC_SIM_MOST_MOORING_BASELINE"])
    model = MostStaticMooring()
    poses = reference["poses"]
    assert poses.shape == (10, 6)
    for pose, expected_force, expected_hv in zip(
        poses, reference["forces"], reference["tensions"], strict=True,
    ):
        force, hv = model.force(pose)
        np.testing.assert_allclose(force[:3], expected_force[:3],
                                   rtol=0, atol=50)
        np.testing.assert_allclose(force[3:], expected_force[3:],
                                   rtol=0, atol=3000)
        np.testing.assert_allclose(hv, expected_hv,
                                   rtol=0, atol=50)
