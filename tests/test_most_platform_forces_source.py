"""Pair active MOST wave and mooring forces on a moving platform."""

import os

import h5py
import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import MostStaticMooring
from wecsim.irregularWave import (
    jonswap_equal_energy_components, synthesize_irregular_response,
)


@pytest.mark.skipif(
    not all(os.environ.get(name) for name in (
        "WEC_SIM_MOST_SHORT_BASELINE", "WEC_SIM_MOST_H5",
        "WEC_SIM_MOST_MASS_PROPERTIES",
    )),
    reason="pinned MATLAB MOST platform force trace not provided",
)
def test_most_platform_wave_and_mooring_against_pinned_source():
    source = loadmat(os.environ["WEC_SIM_MOST_SHORT_BASELINE"])
    h5_file = os.environ["WEC_SIM_MOST_H5"]
    time = source["body_time"].ravel()
    np.testing.assert_allclose(time, source["wave_time"].ravel(),
                               rtol=0, atol=1e-12)
    assert time.shape == (1001,)

    sea = jonswap_equal_energy_components(
        h5_file, significant_height=4, peak_period=8,
        directions=np.array([0.]), spreading=np.array([1.]),
        seed=1, phase_generator="matlab",
    )
    for actual, expected in (
        (sea.omega, source["wave_omega"].ravel()),
        (sea.spectral_amplitude, source["wave_amplitude"].ravel()),
        (sea.d_omega, source["wave_d_omega"].ravel()),
        (sea.phase, source["wave_phase"]),
    ):
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)
    wave = synthesize_irregular_response(
        h5_file, sea, dt=.01, end_time=10, ramp_time=20,
        rho=1025, g=9.80665,
    )
    np.testing.assert_allclose(wave.elevation,
                               source["wave_elevation"].ravel(),
                               rtol=0, atol=1e-11)
    np.testing.assert_allclose(wave.excitation_force,
                               source["body_force_excitation"],
                               rtol=0, atol=1e-5)

    mooring = MostStaticMooring()
    force = np.array([
        mooring.force(pose)[0] for pose in source["mooring_position"]
    ])
    np.testing.assert_allclose(force[:, :3],
                               source["mooring_force"][:, :3],
                               rtol=0, atol=1e-3)
    np.testing.assert_allclose(force[:, 3:],
                               source["mooring_force"][:, 3:],
                               rtol=0, atol=1e-2)

    # Logged WEC-Sim body-force channels use resisting-force signs.
    total = (source["body_force_excitation"]
             - source["body_force_radiation"]
             - source["body_force_added_mass"]
             - source["body_force_restoring"]
             - source["body_force_viscous"]
             - source["body_force_linear_damping"])
    np.testing.assert_allclose(total, source["body_force_total"],
                               rtol=0, atol=1e-5)
    assert np.ptp(source["mooring_position"][:, 0]) > 1

    # The turbine tower-base load is logged in the rotating platform frame.
    # WEC-Sim moves twice the trace of translational infinite-frequency
    # added mass into the Simscape body's scalar mass in this case.
    with h5py.File(os.environ["WEC_SIM_MOST_MASS_PROPERTIES"]) as properties:
        platform_mass = float(properties["Platform/VolturnUS/mass"][0, 0])
    with h5py.File(h5_file) as hydro:
        added_mass = hydro["body1/hydro_coeffs/added_mass/inf_freq"][:]
    adjusted_mass = platform_mass + 2*1025*np.trace(added_mass[:3, :3])
    rotation = np.array([
        mooring._rotation(*pose[3:]) for pose in source["body_position"]
    ])
    tower_force_world = np.einsum(
        "nij,nj->ni", rotation, source["tower_base_load"][:, :3],
    )
    net_force = (source["body_force_total"][:, :3]
                 + source["mooring_force"][:, :3] + tower_force_world)
    np.testing.assert_allclose(
        net_force, adjusted_mass*source["body_acceleration"][:, :3],
        rtol=0, atol=1e-4,
    )
