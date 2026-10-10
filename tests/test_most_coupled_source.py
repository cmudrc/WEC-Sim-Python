"""Pair a source-free six-DOF MOST platform/turbine trajectory with MATLAB."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import (
    MostBaselineController, MostCoupled, MostPlatformHydrodynamics, MostRotor,
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
def test_most_six_dof_coupled_trajectory_against_pinned_source():
    h5_file = os.environ["WEC_SIM_MOST_H5"]
    platform = MostPlatformHydrodynamics.from_volturnus(
        h5_file, os.environ["WEC_SIM_MOST_MASS_PROPERTIES"],
    )
    tower = MostTowerReaction.from_iea15mw(
        os.environ["WEC_SIM_MOST_PROPERTIES"],
        platform_cg=platform.equilibrium_pose[:3],
    )
    controller = None
    if os.environ.get("WEC_SIM_MOST_CONTROL"):
        controller = MostBaselineController.from_matlab_files(
            os.environ["WEC_SIM_MOST_CONTROL"],
            os.environ["WEC_SIM_MOST_STEADY_STATES"], wind_speed=8,
        )
    rotor = MostRotor.from_iea15mw(
        os.environ["WEC_SIM_MOST_BLADE_DIR"], controller=controller,
    )
    wind = MostWindField(read_turbsim_bts(
        Path(os.environ["WEC_SIM_MOST_ADVECTION_DIR"]) / "WIND_8mps.bts",
    ))
    sea = jonswap_equal_energy_components(
        h5_file,
        significant_height=float(os.environ.get("WEC_SIM_MOST_WAVE_HEIGHT", "4")),
        peak_period=8,
        directions=np.array([0.]), spreading=np.array([1.]),
        seed=int(os.environ.get("WEC_SIM_MOST_PHASE_SEED", "1")),
        phase_generator="matlab",
    )
    end_time = float(os.environ.get("WEC_SIM_MOST_END_TIME", "10"))
    assert end_time in (10, 30)
    developed_sea = end_time == 30
    if developed_sea:
        assert controller is not None, "30 s case needs its generated MATLAB controller"
    wave = synthesize_irregular_response(
        h5_file, sea, dt=.01, end_time=end_time, ramp_time=20,
        rho=1025, g=9.80665,
    )
    time = wave.time
    assert time.shape == (round(end_time/.01) + 1,)

    # The source trajectory is used only after the independent solve.
    result = MostCoupled(platform, rotor, tower).simulate(
        time, wave.excitation_force, wind,
    )
    source = loadmat(os.environ["WEC_SIM_MOST_SHORT_BASELINE"])
    if controller is not None:
        np.testing.assert_allclose(
            controller.initial_omega*60/(2*np.pi), source["rotor_speed"][0, 0],
            rtol=0, atol=1e-12,
        )
        initial_state = np.array([0, controller.initial_omega/5, 0])
        np.testing.assert_allclose(
            controller._command(initial_state, controller.initial_pitch)[0],
            source["generator_torque"][0, 0], rtol=0, atol=1e-6,
        )
    np.testing.assert_allclose(time, source["body_time"].ravel(),
                               rtol=0, atol=1e-12)
    for actual, expected in (
        (sea.omega, source["wave_omega"].ravel()),
        (sea.spectral_amplitude, source["wave_amplitude"].ravel()),
        (sea.d_omega, source["wave_d_omega"].ravel()),
        (sea.phase, source["wave_phase"]),
    ):
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-11)
    np.testing.assert_allclose(
        wave.elevation, source["wave_elevation"].ravel(), rtol=0, atol=1e-12,
    )
    np.testing.assert_allclose(
        wave.excitation_force, source["body_force_excitation"],
        rtol=0, atol=1e-4,
    )
    assert result.iterations <= (20 if developed_sea else 12)
    assert result.position_residual <= 1e-6
    assert result.velocity_residual <= 1e-6
    for axis, position_gate, velocity_gate in (
        (0, 6e-4 if developed_sea else 2e-4, 5e-5),  # surge, m and m/s
        (1, 5e-4, 2e-4),    # sway, m and m/s
        (2, 1.5e-4, 7e-5),  # heave, m and m/s
        (3, 4e-5, 2e-5),    # roll, rad and rad/s
        (4, 5e-6, 2e-6),    # pitch, rad and rad/s
        (5, 1.5e-4, 1e-5 if developed_sea else 3e-6),  # yaw, rad and rad/s
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
        rtol=0, atol=.001,
    )
    np.testing.assert_allclose(
        result.rotor.azimuth, source["azimuth"].ravel(),
        rtol=0, atol=.0004 if developed_sea else .0003,
    )
    np.testing.assert_allclose(
        result.rotor.generator_torque, source["generator_torque"].ravel(),
        rtol=0, atol=1000,
    )
