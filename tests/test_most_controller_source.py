"""Pair the published MOST controller with its pinned Simulink block."""

import os

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import MostBaselineController


@pytest.mark.skipif(not os.environ.get("WEC_SIM_MOST_CONTROLLER_BASELINE"),
                    reason="pinned Simulink MOST controller trace not provided")
def test_most_baseline_controller_against_pinned_block():
    source = loadmat(os.environ["WEC_SIM_MOST_CONTROLLER_BASELINE"],
                     simplify_cells=True)
    control = source["control"]
    model = MostBaselineController.iea15mw()
    filt = control["omegaFilter"]
    np.testing.assert_array_equal(filt["A"], [[-5, -5], [1, 0]])
    np.testing.assert_array_equal(filt["B"], [1, 0])
    np.testing.assert_array_equal(filt["C"], [0, 5])
    assert filt["D"] == 0
    for actual, expected in (
        (model.omega_gen, control["omega_gen"]),
        (model.generator_torque, control["Cgen"]),
        (model.sensitivity, [control["c1"], control["c2"], control["c3"]]),
        (model.omega_max, control["omegaMax"]),
        (model.torque_max_rate, control["torqueMaxRate"]),
        (model.pitch_max_rate, control["thetaMaxRate"]),
        (model.kp, control["KP"]),
        (model.ki, control["KI"]),
        (model.initial_omega, source["omega0"]),
        (model.initial_pitch, source["pitch0"]),
    ):
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-9)
    time = source["time"]
    assert time.shape == source["speed"].shape == (4001,)
    torque, pitch = model.simulate(time, source["speed"])
    np.testing.assert_allclose(torque, source["torque"], rtol=0, atol=1e-4)
    np.testing.assert_allclose(pitch, source["pitch"], rtol=0, atol=1e-9)
    assert np.max(np.abs(np.diff(source["torque"])/np.diff(time))) > 4.49e6
    assert np.max(np.abs(np.diff(source["pitch"])/np.diff(time))) > 0.122
