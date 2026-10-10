"""Advance the MOST rotor independently against a pinned coupled source run."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import (
    MostBaselineController, MostRotor, MostWindField, read_turbsim_bts,
)


@pytest.mark.skipif(
    not all(os.environ.get(name) for name in (
        "WEC_SIM_MOST_SHORT_BASELINE", "WEC_SIM_MOST_ADVECTION_DIR",
        "WEC_SIM_MOST_BLADE_DIR", "WEC_SIM_MOST_PROPERTIES",
    )),
    reason="pinned MATLAB MOST short-run inputs not provided",
)
def test_most_rotor_against_pinned_coupled_source():
    source = loadmat(os.environ["WEC_SIM_MOST_SHORT_BASELINE"])
    time = source["turbine_time"].ravel()
    np.testing.assert_array_equal(time, source["body_time"].ravel())
    assert time.shape == (1001,)
    wind = MostWindField(read_turbsim_bts(
        Path(os.environ["WEC_SIM_MOST_ADVECTION_DIR"]) / "WIND_8mps.bts",
    ))
    controller = None
    if os.environ.get("WEC_SIM_MOST_CONTROL"):
        controller = MostBaselineController.from_matlab_files(
            os.environ["WEC_SIM_MOST_CONTROL"],
            os.environ["WEC_SIM_MOST_STEADY_STATES"], wind_speed=8,
        )
    rotor = MostRotor.from_iea15mw(
        os.environ["WEC_SIM_MOST_BLADE_DIR"], controller=controller,
    )
    properties = loadmat(os.environ["WEC_SIM_MOST_PROPERTIES"],
                         simplify_cells=True)["WTcomponents"]
    np.testing.assert_allclose(rotor.inertia,
                               properties["Inertia_Rotor_cogHub"],
                               rtol=0, atol=1e-6)
    sampled_wind = np.array([
        wind.sampler_at(float(t))([0, 0, 150]) for t in time
    ])
    np.testing.assert_allclose(sampled_wind, source["wind_speed"],
                               rtol=0, atol=1.5e-5)

    response = rotor.simulate(
        time, source["body_position"], source["body_velocity"], wind,
    )
    assert response.blade_root_load.shape == (1001, 6, 3)
    np.testing.assert_allclose(
        response.blade_root_load[0], source["blade_aero_load"][0],
        rtol=0, atol=1e-4,
    )
    np.testing.assert_allclose(
        response.rotor_speed * 60/(2*np.pi), source["rotor_speed"].ravel(),
        rtol=0, atol=1e-3,
    )
    np.testing.assert_allclose(response.azimuth, source["azimuth"].ravel(),
                               rtol=0, atol=1e-3)
    np.testing.assert_allclose(
        response.generator_torque, source["generator_torque"].ravel(),
        rtol=0, atol=2e3,
    )
    np.testing.assert_allclose(response.blade_pitch,
                               source["blade_pitch"].ravel(),
                               rtol=0, atol=1e-8)
    assert np.ptp(source["body_position"][:, 0]) > 1
    assert np.ptp(source["rotor_speed"]) > .8
