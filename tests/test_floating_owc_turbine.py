"""Pair floating OWC turbine control, power, and speed with MATLAB."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from wecsim import FloatingOwcChamber, FloatingOwcTurbine


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")


def test_turbine_load_limits_and_overspeed():
    turbine = FloatingOwcTurbine()
    stopped_air = turbine.evaluate(0, 150)
    assert stopped_air.turbine_torque == 0
    assert stopped_air.speed_derivative < 0
    assert turbine.evaluate(1000, 1200).control_torque == 216.5
    with pytest.raises(ValueError, match="rotor speed"):
        turbine.evaluate(0, 0)


@pytest.mark.skipif(not REFERENCE, reason="pinned floating OWC turbine traces absent")
def test_published_turbine_outputs_and_coupled_airtrain_on_source_motion():
    def read(name):
        return np.loadtxt(
            Path(REFERENCE) / f"FLOATING_OWC_SOURCE_{name}.csv",
            delimiter=",", ndmin=2,
        )

    pressure = read("x_deltaP")
    speed = read("x_vTurb")
    displacement = read("x_xOWC")
    # Bus input four is the water-column heave speed fed to velOWC.
    column_speed = read("x_signal4")
    time = pressure[:, 0]
    assert len(time) == 50001 and time[-1] == 500
    for signal in (speed, displacement, column_speed):
        np.testing.assert_allclose(signal[:, 0], time, rtol=0, atol=1e-10)
    turbine = FloatingOwcTurbine()
    on_source_states = turbine.evaluate(pressure[:, 1], speed[:, 1])
    source_outputs = {}
    for actual, name, limit in (
        (on_source_states.control_torque, "u", 1e-10),
        (on_source_states.load_power, "P_turb", 1e-8),
        (on_source_states.pneumatic_power, "P_pneumatic", 1e-7),
        (on_source_states.efficiency, "eta_turb", 1e-12),
    ):
        logged = read(name)
        np.testing.assert_allclose(logged[:, 0], time, rtol=0, atol=1e-10)
        source_outputs[name] = logged[:, 1]
        assert np.max(np.abs(actual - logged[:, 1])) < limit, name

    area = np.pi * 5.89**2 / 4
    chamber = FloatingOwcChamber(
        area, area * 4.5, 1.4, 101325, 1.25, 0.75, 0.775)

    def derivative(t, state):
        return [
            chamber.pressure_derivative(
                state[0], np.interp(t, time, displacement[:, 1]),
                np.interp(t, time, column_speed[:, 1]), state[1],
            ),
            turbine.evaluate(state[0], state[1]).speed_derivative,
        ]

    predicted = solve_ivp(
        derivative, (time[0], time[-1]), [pressure[0, 1], speed[0, 1]],
        t_eval=time, rtol=1e-8, atol=1e-8, max_step=0.02,
    )
    assert predicted.success
    assert np.isfinite(predicted.y).all()
    pressure_error = np.abs(predicted.y[0] - pressure[:, 1])
    speed_error = np.abs(predicted.y[1] - speed[:, 1])
    assert np.max(pressure_error) < 0.05
    assert np.sqrt(np.mean(pressure_error**2)) < 0.01
    assert np.max(speed_error) < 0.002
    assert np.sqrt(np.mean(speed_error**2)) < 0.001
    coupled = turbine.evaluate(predicted.y[0], predicted.y[1])
    for actual, name, limit in (
        (coupled.control_torque, "u", 2e-4),
        (coupled.load_power, "P_turb", 0.1),
        (coupled.pneumatic_power, "P_pneumatic", 0.5),
        (coupled.efficiency, "eta_turb", 5e-5),
    ):
        assert np.max(np.abs(actual - source_outputs[name])) < limit, name
