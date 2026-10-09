"""Pair published IEA 15 MW BEM inputs and blade loads with MATLAB MOST."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import MostBEM


@pytest.mark.skipif(
    not (os.environ.get("WEC_SIM_MOST_BEM_BASELINE")
         and os.environ.get("WEC_SIM_MOST_BLADE_DIR")),
    reason="pinned MATLAB MOST BEM and blade inputs not provided",
)
def test_most_bem_against_pinned_library():
    source = loadmat(os.environ["WEC_SIM_MOST_BEM_BASELINE"], simplify_cells=True)
    assert source["q"].shape == (5, 14)
    assert source["loads"].shape == (5, 6, 3)
    assert np.isfinite(source["loads"]).all()
    model = MostBEM.from_iea15mw(Path(os.environ["WEC_SIM_MOST_BLADE_DIR"]))
    bem = source["bem"]
    for attribute, name in (
        ("radius", "bladedata_r"),
        ("radius_interval", "bladedata_r_int"),
        ("twist", "bladedata_twist"),
        ("chord", "bladedata_chord"),
        ("curve_angle", "bladedata_BlCrvAng"),
        ("sweep", "bladedata_BlSwpAC"),
        ("curve", "bladedata_BlCrvAC"),
        ("airfoil", "bladedata_airfoil"),
    ):
        np.testing.assert_allclose(
            getattr(model, attribute), np.asarray(bem[name]).ravel()
            if attribute != "airfoil" else bem[name],
            rtol=0, atol=1e-12,
        )
    np.testing.assert_array_equal(model.airfoil_index + 1,
                                  bem["bladedata_airfoil_index"].ravel())
    np.testing.assert_allclose([model.hub_radius, model.tip_radius],
                               bem["bladedata_RhubRblade"].ravel(),
                               rtol=0, atol=1e-12)
    for i, (state, speed, pitch) in enumerate(zip(
        source["q"], source["wind_speeds"], source["bladepitch"], strict=True,
    )):
        if i == 4:
            wind = lambda point: np.array([
                8 + .002*point[0] + .004*point[1] + .01*(point[2]-150),
                .1 + .001*point[1], .02*point[2]/150,
            ])
        else:
            wind = [speed, 0, 0]
        actual = model.loads(state, pitch, wind)
        np.testing.assert_allclose(actual, source["loads"][i],
                                   rtol=0, atol=1e-4)
