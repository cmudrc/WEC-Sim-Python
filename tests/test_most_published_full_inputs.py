"""Pair the published 1,000 s MOST sea with a fresh MATLAB source run."""

import os

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import MostBaselineController, MostWindField, read_turbsim_bts
from wecsim.irregularWave import (
    jonswap_equal_energy_components, synthesize_irregular_response,
)


@pytest.mark.skipif(
    not all(os.environ.get(name) for name in (
        "WEC_SIM_MOST_FULL_SOURCE", "WEC_SIM_MOST_H5", "WEC_SIM_MOST_BTS",
        "WEC_SIM_MOST_CONTROL", "WEC_SIM_MOST_STEADY_STATES",
    )),
    reason="pinned 1,000 s MATLAB MOST source inputs not provided",
)
def test_most_published_full_duration_inputs_against_source():
    source = loadmat(os.environ["WEC_SIM_MOST_FULL_SOURCE"])
    hydro = os.environ["WEC_SIM_MOST_H5"]
    sea = jonswap_equal_energy_components(
        hydro, significant_height=4, peak_period=8,
        directions=np.array([0.]), spreading=np.array([1.]),
        seed=1, phase_generator="matlab",
    )
    for actual, expected in (
        (sea.omega, source["wave_omega"].ravel()),
        (sea.spectral_amplitude, source["wave_amplitude"].ravel()),
        (sea.d_omega, source["wave_d_omega"].ravel()),
        (sea.phase, source["wave_phase"]),
    ):
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-11)

    wave = synthesize_irregular_response(
        hydro, sea, dt=.01, end_time=1000, ramp_time=20,
        rho=1025, g=9.80665,
    )
    assert wave.time.size == 100001
    np.testing.assert_allclose(
        wave.time, source["wave_time"].ravel(), rtol=0, atol=1e-10,
    )
    np.testing.assert_allclose(
        wave.time[::10], source["body_time"].ravel(), rtol=0, atol=1e-10,
    )
    np.testing.assert_allclose(
        wave.elevation, source["wave_elevation"].ravel(),
        rtol=0, atol=1e-12,
    )
    np.testing.assert_allclose(
        wave.excitation_force[::10], source["body_force_excitation"],
        rtol=0, atol=1e-4,
    )

    controller = MostBaselineController.from_matlab_files(
        os.environ["WEC_SIM_MOST_CONTROL"],
        os.environ["WEC_SIM_MOST_STEADY_STATES"], wind_speed=8,
    )
    np.testing.assert_allclose(
        controller.initial_omega * 60 / (2 * np.pi),
        source["rotor_speed"][0, 0], rtol=0, atol=1e-12,
    )
    initial_state = np.array([0., controller.initial_omega / 5, 0.])
    np.testing.assert_allclose(
        controller._command(initial_state, controller.initial_pitch)[0],
        source["generator_torque"][0, 0], rtol=0, atol=1e-6,
    )

    wind = MostWindField(read_turbsim_bts(os.environ["WEC_SIM_MOST_BTS"]))
    assert wind.n_time == 19999
    assert 999.8 < wind.time[-1] < 1000
    probe = np.array([0., 0., 150.])
    np.testing.assert_array_equal(
        wind.sampler_at(1000)(probe), wind.sampler_at(float(wind.time[-1]))(probe),
    )
