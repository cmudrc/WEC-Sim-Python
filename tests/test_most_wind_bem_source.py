"""Pair TurbSim advection, spatial sampling, and rotor loads with MATLAB MOST."""

import os

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import MostBEM, MostWindField, read_turbsim_bts


@pytest.mark.skipif(
    not all(os.environ.get(name) for name in (
        "WEC_SIM_MOST_BEM_ADVECTION_BASELINE", "WEC_SIM_MOST_ADVECTION_DIR",
        "WEC_SIM_MOST_BLADE_DIR",
    )),
    reason="pinned MATLAB MOST wind/BEM inputs not provided",
)
def test_most_turbsim_wind_to_bem_against_pinned_matlab():
    source = loadmat(os.environ["WEC_SIM_MOST_BEM_ADVECTION_BASELINE"])
    field = MostWindField(read_turbsim_bts(
        os.path.join(os.environ["WEC_SIM_MOST_ADVECTION_DIR"], "WIND_8mps.bts"),
    ))
    rotor = MostBEM.from_iea15mw(os.environ["WEC_SIM_MOST_BLADE_DIR"])
    assert source["loads"].shape == (3, 6, 3)
    assert source["sampled_wind"].shape == (3, 3, 3)
    for case, index in enumerate(source["indices"].ravel()):
        sampler = field.sampler(int(index)-1)
        sampled = np.array([sampler(point) for point in source["points"]])
        np.testing.assert_allclose(
            sampled, source["sampled_wind"][case], rtol=0, atol=1e-12,
        )
        actual = rotor.loads(source["q"][case],
                             source["bladepitch"].ravel()[case], sampler)
        np.testing.assert_allclose(
            actual, source["loads"][case], rtol=0, atol=1e-4,
        )
