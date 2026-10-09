"""MOST's wind-plane advection against the pinned MATLAB preprocessor."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import MostWindField, read_turbsim_bts


def test_most_wind_advection_sample_order_and_bounds():
    # One component carries its time index so X-plane shifts are observable.
    from wecsim.turbSim import TurbSimWind

    velocity = np.broadcast_to(np.arange(30)[:, None, None, None],
                               (30, 3, 2, 2)).copy()
    wind = TurbSimWind(velocity, np.empty((30, 3, 0)),
                       np.array([-100., 100.]), np.array([0., 1.]),
                       1., 10., 8., "")
    field = MostWindField(wind, speed=10)
    assert field.discarded == 10
    assert field.n_time == 11
    first = field.at_index(0)
    assert first.shape == (3, 10, 2, 2)
    np.testing.assert_array_equal(first[0, :, 0, 0], np.arange(13, 3, -1))
    with pytest.raises(IndexError):
        field.at_index(field.n_time)


@pytest.mark.skipif(not os.environ.get("WEC_SIM_MOST_ADVECTION_DIR"),
                    reason="pinned MATLAB MOST advection output not provided")
def test_most_wind_planes_against_run_turbsim():
    directory = Path(os.environ["WEC_SIM_MOST_ADVECTION_DIR"])
    wind = read_turbsim_bts(directory / "WIND_8mps.bts")
    field = MostWindField(wind)
    reference = loadmat(directory.parent / "matlab-most-advection.mat")
    values = reference["values"]
    assert values.shape == (249, 3, 10, 12, 12)
    assert field.n_time == 249 and field.discarded == 376
    np.testing.assert_allclose(field.time, reference["time"].ravel(), rtol=0, atol=1e-12)
    np.testing.assert_array_equal(field.x, reference["x"].ravel())
    np.testing.assert_array_equal(wind.y, reference["y"].ravel())
    np.testing.assert_array_equal(wind.z, reference["z"].ravel())
    for index in range(field.n_time):
        np.testing.assert_allclose(field.at_index(index), values[index],
                                   rtol=0, atol=1e-12)
