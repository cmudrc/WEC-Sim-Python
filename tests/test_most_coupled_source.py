"""Pair a source-free reduced MOST platform/turbine trajectory with MATLAB."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import (
    MostCoupled, MostPlatformHydrodynamics, MostRotor,
    MostTowerReaction, MostWindField, read_turbsim_bts,
)
from wecsim.irregularWave import (
    jonswap_equal_energy_components, synthesize_irregular_response,
)


@pytest.mark.skipif(
    not all(os.environ.get(name) for name in (
        "WEC_SIM_MOST_SHORT_BASELINE", "WEC_SIM_MOST_H5",
        "WEC_SIM_MOST_MASS_PROPERTIES", "WEC_SIM_MOST_PROPERTIES",
        "WEC_SIM_MOST_BLADE_DIR", "WEC_SIM_MOST_ADVECTION_DIR",
    )),
    reason="pinned MATLAB MOST coupled-run inputs not provided",
)
def test_most_reduced_coupled_trajectory_against_pinned_source():
    h5_file = os.environ["WEC_SIM_MOST_H5"]
    platform = MostPlatformHydrodynamics.from_volturnus(
        h5_file, os.environ["WEC_SIM_MOST_MASS_PROPERTIES"],
    )
    tower = MostTowerReaction.from_iea15mw(
        os.environ["WEC_SIM_MOST_PROPERTIES"],
        platform_cg=platform.equilibrium_pose[:3],
    )
    rotor = MostRotor.from_iea15mw(os.environ["WEC_SIM_MOST_BLADE_DIR"])
    wind = MostWindField(read_turbsim_bts(
        Path(os.environ["WEC_SIM_MOST_ADVECTION_DIR"]) / "WIND_8mps.bts",
    ))
    sea = jonswap_equal_energy_components(
        h5_file, significant_height=4, peak_period=8,
        directions=np.array([0.]), spreading=np.array([1.]),
        seed=1, phase_generator="matlab",
    )
    wave = synthesize_irregular_response(
        h5_file, sea, dt=.01, end_time=10, ramp_time=20,
        rho=1025, g=9.80665,
    )
    time = wave.time
    assert time.shape == (1001,)

    # The source trajectory is used only after the independent solve.
    result = MostCoupled(platform, rotor, tower).simulate(
        time, wave.excitation_force, wind,
    )
    source = loadmat(os.environ["WEC_SIM_MOST_SHORT_BASELINE"])
    np.testing.assert_allclose(time, source["body_time"].ravel(),
                               rtol=0, atol=1e-12)
    assert result.iterations <= 12
    assert result.position_residual <= 1e-6
    assert result.velocity_residual <= 1e-6
    for axis, position_gate, velocity_gate in (
        (0, 4e-4, 1e-4),    # surge, m and m/s
        (2, 1.5e-4, 7e-5),  # heave, m and m/s
        (4, 4e-5, 1e-5),    # pitch, rad and rad/s
    ):
        np.testing.assert_allclose(
            result.platform.position[:, axis], source["body_position"][:, axis],
            rtol=0, atol=position_gate,
        )
        np.testing.assert_allclose(
            result.platform.velocity[:, axis], source["body_velocity"][:, axis],
            rtol=0, atol=velocity_gate,
        )
    np.testing.assert_allclose(
        result.rotor.rotor_speed*60/(2*np.pi), source["rotor_speed"].ravel(),
        rtol=0, atol=.01,
    )
    np.testing.assert_allclose(
        result.rotor.azimuth, source["azimuth"].ravel(),
        rtol=0, atol=.006,
    )
    np.testing.assert_allclose(
        result.rotor.generator_torque, source["generator_torque"].ravel(),
        rtol=0, atol=25000,
    )
