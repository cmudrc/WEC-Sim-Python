"""Check the pinned WaveStar NMPC application's sea and mechanical signals."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim.irregularWave import (
    jonswap_equal_energy_components, synthesize_irregular_response,
)
from wecsim.wavestar import WaveStarLinkage, run_wavestar_published
from wecsim import WaveStarNmpcActuator, WaveStarNmpcObserver


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
HYDRO = os.environ.get("WEC_SIM_WAVESTAR_H5")
pytestmark = pytest.mark.skipif(
    not (REFERENCE and HYDRO),
    reason="pinned WaveStar NMPC output and hydrodynamics are absent",
)


def _source(name):
    values = np.loadtxt(
        Path(REFERENCE) / f"WECCCOMP_NMPC_SOURCE_{name}.csv",
        delimiter=",", ndmin=2,
    )
    assert np.isfinite(values).all()
    return values


def test_nmpc_published_sea_and_excitation():
    saved = _source("components")
    wave = _source("wave")
    source_float = _source("WECCCOMP_Nonlinear_Model_Predictive_body1")
    assert saved.shape == (500, 4)
    assert wave.shape == (4501, 2)
    assert source_float.shape[1] == 25
    np.testing.assert_allclose(wave[:, 0], np.arange(4501) * .05,
                               rtol=0, atol=1e-9)
    components = jonswap_equal_energy_components(
        HYDRO, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=3.3, phase=saved[:, 3, None],
    )
    np.testing.assert_allclose(components.omega, saved[:, 0],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(components.spectral_amplitude, saved[:, 1],
                               rtol=0, atol=1e-11)
    np.testing.assert_allclose(components.d_omega, saved[:, 2],
                               rtol=0, atol=1e-10)
    incident = synthesize_irregular_response(
        HYDRO, components, dt=.05, end_time=225, ramp_time=25,
        rho=1000, g=9.80665,
    )
    np.testing.assert_allclose(incident.elevation, wave[:, 1],
                               rtol=0, atol=1e-10)
    indices = np.rint(source_float[:, 0] / .05).astype(int)
    np.testing.assert_allclose(source_float[:, 0], wave[indices, 0],
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(
        incident.excitation_force[indices], source_float[:, 19:25],
        rtol=0, atol=1e-7,
    )


def test_nmpc_published_linkage_and_control_trace():
    source_float = _source("WECCCOMP_Nonlinear_Model_Predictive_body1")
    arm = _source("WECCCOMP_Nonlinear_Model_Predictive_body2")
    pto = _source("WECCCOMP_Nonlinear_Model_Predictive_pto1")
    assert arm.shape == source_float.shape
    assert pto.shape[1] == 43
    np.testing.assert_allclose(arm[:, 0], source_float[:, 0],
                               rtol=0, atol=1e-9)
    linkage = WaveStarLinkage()
    angle = source_float[:, 5]
    speed = source_float[:, 11]
    for source in (source_float, arm):
        neutral = source[0, [1, 3]]
        np.testing.assert_allclose(
            linkage.point_position(neutral, angle),
            source[:, [1, 3]], rtol=0, atol=1e-10,
        )
        np.testing.assert_allclose(
            linkage.point_velocity(neutral, angle, speed),
            source[:, [7, 9]], rtol=0, atol=1e-9,
        )
    np.testing.assert_allclose(linkage.pto_stroke(angle), pto[:, 3],
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(linkage.pto_speed(angle, speed), pto[:, 9],
                               rtol=0, atol=1e-9)
    assert np.max(np.abs(pto[pto[:, 0] > 15, 33])) > 1

    controller = loadmat(
        Path(REFERENCE) / "WECCCOMP_NMPC_SOURCE_controller.mat",
        simplify_cells=True,
    )
    for name in ("cmd_ptoM", "Output_power", "Output_energy",
                 "estimated_states", "AR_excM_pred", "excM_wecSim",
                 "motor_displacement"):
        assert name in controller, f"source controller trace {name} missing"


def test_nmpc_actuator_on_source_command_and_stroke():
    pto = _source("WECCCOMP_Nonlinear_Model_Predictive_pto1")
    controller = loadmat(
        Path(REFERENCE) / "WECCCOMP_NMPC_SOURCE_controller.mat",
        simplify_cells=True,
    )
    command = np.asarray(controller["cmd_ptoM"]["signals"]["values"])
    stroke = np.asarray(
        controller["motor_displacement"]["signals"]["values"]
    )
    time = np.asarray(controller["cmd_ptoM"]["time"])
    assert command.shape == stroke.shape == time.shape == (4501,)
    np.testing.assert_allclose(time, np.arange(4501) * .05,
                               rtol=0, atol=1e-9)
    assert np.max(np.abs(command[time < 10])) < 1e-12
    assert np.max(np.abs(command)) <= 12 + 1e-10
    assert np.max(np.abs(command[time >= 15])) > 1
    estimated = np.asarray(
        controller["estimated_states"]["signals"]["values"]
    )
    assert estimated.shape == (4501, 5)
    resistive = (time >= 10) & (time < 15)
    np.testing.assert_allclose(
        command[resistive], np.clip(-19.4 * estimated[resistive, 1], -12, 12),
        rtol=0, atol=1e-12,
    )
    actuator = WaveStarNmpcActuator()
    force = np.array([
        actuator.step(request, position)
        for request, position in zip(command, stroke)
    ])
    indices = np.rint(pto[:, 0] / .05).astype(int)
    np.testing.assert_allclose(stroke[indices], pto[:, 3],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(force[indices], pto[:, 33],
                               rtol=0, atol=1e-7)


def test_nmpc_observer_on_source_command_and_stroke():
    controller = loadmat(
        Path(REFERENCE) / "WECCCOMP_NMPC_SOURCE_controller.mat",
        simplify_cells=True,
    )
    command = np.asarray(controller["cmd_ptoM"]["signals"]["values"])
    stroke = np.asarray(
        controller["motor_displacement"]["signals"]["values"]
    )
    saved = np.asarray(
        controller["estimated_states"]["signals"]["values"]
    )
    assert command.shape == stroke.shape == (4501,)
    assert saved.shape == (4501, 5)
    observer = WaveStarNmpcObserver()
    replay = np.array([
        observer.step(position, command[index - 1] if index else 0.)
        for index, position in enumerate(stroke)
    ])
    np.testing.assert_allclose(replay, saved, rtol=0, atol=1e-10)


def test_nmpc_unforced_plant_motion_diagnostic():
    """Bound the pre-control solver difference without claiming NMPC parity."""
    source_float = _source("WECCCOMP_Nonlinear_Model_Predictive_body1")
    source_pto = _source("WECCCOMP_Nonlinear_Model_Predictive_pto1")
    source_forces = _source("body1_forces")
    assert source_float.shape == (4501, 25)
    assert source_pto.shape == (4501, 43)
    assert source_forces.shape == (4501, 37)
    components = _source("components")
    sea = jonswap_equal_energy_components(
        HYDRO, significant_height=.1042, peak_period=1.836,
        directions=np.array([0.]), spreading=np.array([1.]),
        gamma=3.3, phase=components[:, 3, None],
    )
    # The published source uses ode8 at 0.05 s. A finer Python step keeps
    # the fitted radiation system stable; the numerical paths still differ.
    response = run_wavestar_published(
        HYDRO, sea, dt=.001, end_time=9.95, ramp_time=25,
        g=9.80665, output_stride=50,
    )
    count = len(response.time)
    np.testing.assert_allclose(response.time, source_float[:count, 0],
                               rtol=0, atol=1e-9)
    assert np.max(np.abs(source_pto[:count, 31:37])) < 1e-12
    assert np.max(np.abs(response.angle - source_float[:count, 5])) < .004
    assert np.max(np.abs(
        response.angular_speed - source_float[:count, 11]
    )) < .04
    assert np.max(np.abs(
        response.float_position[:, [0, 2]]
        - source_float[:count][:, [1, 3]]
    )) < .002
    assert np.max(np.abs(
        response.pto_stroke - source_pto[:count, 3]
    )) < .001
    assert np.max(np.abs(
        response.radiation_force - source_forces[:count, 1:7]
    )) < .2
