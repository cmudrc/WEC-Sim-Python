"""Pair the published floating OWC chamber state with a Python pressure ODE."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from wecsim import FloatingOwcChamber


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")


def _published_chamber():
    diameter = 5.89
    area = np.pi * diameter**2 / 4
    return FloatingOwcChamber(
        area=area,
        initial_volume=area * 4.5,
        gamma=1.4,
        ambient_pressure=101325,
        ambient_density=1.25,
        turbine_diameter=0.75,
        turbine_kappa=0.775,
    )


def test_chamber_restoring_pressure_and_turbine_venting():
    chamber = _published_chamber()
    assert chamber.pressure_derivative(0, 0, 0, 150) == 0
    assert chamber.pressure_derivative(0, 0, 0.1, 150) > 0
    assert chamber.pressure_derivative(0, 0, -0.1, 150) < 0
    assert chamber.pressure_derivative(1000, 0, 0, 150) < 0
    assert chamber.force_on_column(1000) < 0
    with pytest.raises(ValueError, match="volume"):
        chamber.volume(5)
    with pytest.raises(ValueError, match="nonzero turbine speed"):
        chamber.pressure_derivative(0, 0, 0, 0)


@pytest.mark.skipif(not REFERENCE, reason="pinned floating OWC state bus absent")
def test_published_chamber_pressure_on_source_motion_and_turbine_speed():
    def read(name):
        return np.loadtxt(
            Path(REFERENCE) / f"FLOATING_OWC_SOURCE_x_{name}.csv",
            delimiter=",", ndmin=2,
        )

    displacement = read("xOWC")
    speed = read("vOWC")
    turbine_speed = read("vTurb")
    pressure = read("deltaP")
    time = pressure[:, 0]
    for signal in (displacement, speed, turbine_speed):
        np.testing.assert_allclose(signal[:, 0], time, rtol=0, atol=1e-10)
    assert len(time) == 50001 and time[-1] == 500
    chamber = _published_chamber()

    def derivative(t, state):
        return [chamber.pressure_derivative(
            state[0],
            np.interp(t, time, displacement[:, 1]),
            np.interp(t, time, speed[:, 1]),
            np.interp(t, time, turbine_speed[:, 1]),
        )]

    predicted = solve_ivp(
        derivative, (time[0], time[-1]), [pressure[0, 1]],
        t_eval=time, rtol=1e-8, atol=1e-8, max_step=0.02,
    )
    assert predicted.success
    assert np.isfinite(predicted.y).all()
    error = np.abs(predicted.y[0] - pressure[:, 1])
    assert np.max(error) < 100, f"chamber pressure error {np.max(error):.6g} Pa"
