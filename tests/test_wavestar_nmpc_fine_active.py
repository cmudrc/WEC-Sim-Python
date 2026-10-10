"""Keep the fine-step WaveStar NMPC activation diagnostic reproducible."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import WaveStarNmpcController, WaveStarNmpcObserver


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="fine-step WaveStar MATLAB trace is absent",
)
PREFIX = "WECCCOMP_NMPC_FINE_ACTIVE_DIAG_"


def test_fine_step_nmpc_activation_diagnostic():
    source = Path(REFERENCE)
    log = loadmat(source / f"{PREFIX}controller.mat", simplify_cells=True)
    time = np.asarray(log["cmd_ptoM"]["time"])
    command = np.asarray(log["cmd_ptoM"]["signals"]["values"])
    stroke = np.asarray(log["motor_displacement"]["signals"]["values"])
    estimated = np.asarray(log["estimated_states"]["signals"]["values"])
    forecast = np.asarray(log["AR_excM_pred"]["signals"]["values"]).T
    assert time.shape == command.shape == stroke.shape == (15201,)
    assert estimated.shape == (15201, 5)
    assert forecast.shape == (15201, 40)
    np.testing.assert_allclose(time, np.arange(len(time)) * .001,
                               rtol=0, atol=1e-9)

    observer = WaveStarNmpcObserver(dt=.001)
    replay_state = np.array([
        observer.step(position, command[index - 1] if index else 0.)
        for index, position in enumerate(stroke)
    ])
    np.testing.assert_allclose(replay_state, estimated, rtol=0, atol=1e-10)
    control = WaveStarNmpcController(dt=.001)
    replay_command = np.array([
        control.step(estimated[index], forecast[index], instant)
        for index, instant in enumerate(time)
    ])
    np.testing.assert_allclose(replay_command, command, rtol=0, atol=1e-4)

    # This probes MATLAB's ar/forecast on the saved estimate window, rather
    # than treating the ill-conditioned derived run as a motion parity gate.
    probe = np.loadtxt(source / f"{PREFIX}forecast_probe.csv", delimiter=",")
    assert probe.shape == (6, 81)
    np.testing.assert_allclose(
        probe[:, 0], [15, 15.001, 15.01, 15.05, 15.1, 15.2],
        rtol=0, atol=1e-9,
    )
    np.testing.assert_allclose(probe[:, 1:41], probe[:, 41:81],
                               rtol=0, atol=1e-10)
    assert np.max(np.abs(forecast[time >= 15])) > 5e4
    moment = estimated[15001-180:15001, -1]
    order = 18
    forward = np.array([
        moment[t-order:t][::-1] for t in range(order, len(moment))
    ])
    backward = np.array([
        moment[t-order+1:t+1] for t in range(order, len(moment))
    ])
    condition = np.linalg.cond(np.vstack((forward, backward)))
    assert condition > 1e13
