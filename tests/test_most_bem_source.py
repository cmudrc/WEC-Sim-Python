"""Check the pinned MOST library BEM source baseline before Python pairing."""

import os

import numpy as np
import pytest
from scipy.io import loadmat


@pytest.mark.skipif(not os.environ.get("WEC_SIM_MOST_BEM_BASELINE"),
                    reason="pinned MATLAB MOST BEM output not provided")
def test_most_bem_source_baseline_is_finite():
    source = loadmat(os.environ["WEC_SIM_MOST_BEM_BASELINE"])
    assert source["q"].shape == (4, 14)
    assert source["loads"].shape == (4, 6, 3)
    assert np.isfinite(source["loads"]).all()
    assert np.max(np.abs(source["loads"])) > 0
