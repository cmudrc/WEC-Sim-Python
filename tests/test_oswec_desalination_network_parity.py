"""Independent high-pressure network against pinned 300 s Desalination."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import (
    DynamicPressureReliefValve, GasChargedAccumulator,
    ReverseOsmosisHydraulicNetwork, ReverseOsmosisMembrane,
)


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="pinned OSWEC Desalination output absent",
)


def _source(name):
    return np.loadtxt(
        Path(REFERENCE) / f"OSWEC_DESALINATION_SOURCE_{name}.csv",
        delimiter=",",
    )


def _network():
    return ReverseOsmosisHydraulicNetwork(
        accumulator=GasChargedAccumulator(
            initial_gas_volume=4, precharge_pressure=3e6,
            atmospheric_pressure=101_325, dead_gas_volume=4e-5,
            hard_stop_stiffness=1e7, hard_stop_damping=1e7,
        ),
        membrane=ReverseOsmosisMembrane(
            resistance=.6023e8, osmotic_pressure=3e6,
            regulation_range=5e4, max_area=.3, leakage_area=1e-12,
            discharge_coefficient=.7, fluid_density=850,
        ),
        relief=DynamicPressureReliefValve(
            set_pressure=5.6e6, regulation_range=1e5,
            max_area=.0267, leakage_area=1e-12,
            discharge_coefficient=.7, fluid_density=850,
            opening_time_constant=.1,
        ),
        cylinder_area=.26,
        pump_to_motor_displacement_ratio=.95,
        pump_outlet_resistance=.6023e8 / 22,
    )


def test_published_relief_and_pressure_recovery_laws():
    network = _network()
    hydraulic, relief = _source("simout1"), _source("relief")
    pressure, time = hydraulic[:, 6], hydraulic[:, 0]
    area = np.empty_like(pressure)
    area[0] = relief[0, 2]
    for index in range(1, len(pressure)):
        area[index] = network.relief.next_area(
            pressure[index], area[index - 1], time[index] - time[index - 1],
        )
    np.testing.assert_allclose(area, relief[:, 2], rtol=0, atol=1e-12)
    actual_relief = np.array([
        network.relief.flow(p, a) for p, a in zip(pressure, area)
    ])
    np.testing.assert_allclose(actual_relief, relief[:, 3], rtol=0, atol=1e-10)
    ratio = network.pump_to_motor_displacement_ratio
    brine = pressure * (1 - ratio) / (
        ratio * ratio * network.pump_outlet_resistance
    )
    np.testing.assert_allclose(brine, hydraulic[:, 3], rtol=0, atol=1e-12)
    np.testing.assert_allclose(ratio * brine, hydraulic[:, 4],
                               rtol=0, atol=1e-12)


def test_published_network_from_source_feed_flow():
    hydraulic = _source("simout1")
    network = _network()
    state = network.initial_state(pressure=hydraulic[0, 6])
    predicted = np.empty(6001)
    predicted[0] = state.pressure
    for index in range(1, len(predicted)):
        state = network.step_from_feed_flow(
            hydraulic[index, 1], state,
            hydraulic[index, 0] - hydraulic[index - 1, 0],
        )
        predicted[index] = state.pressure
    np.testing.assert_allclose(predicted, hydraulic[:len(predicted), 6],
                               rtol=0, atol=1e-2)


def test_published_network_from_rod_speed_alone():
    hydraulic, cylinder, relief = (
        _source("simout1"), _source("cylinder"), _source("relief")
    )
    assert hydraulic.shape == (30_001, 8)
    assert cylinder.shape == (30_001, 8)
    network = _network()
    state = network.initial_state()
    predicted = np.empty((len(hydraulic), 6))
    predicted[0] = 0
    for index in range(1, len(hydraulic)):
        state = network.step(
            cylinder[index, 2], state,
            hydraulic[index, 0] - hydraulic[index - 1, 0],
        )
        predicted[index] = (
            state.pressure, state.permeate_flow, state.brine_flow,
            state.recovered_feed_flow, state.accumulator_flow,
            state.relief_flow,
        )
    references = np.column_stack((
        hydraulic[:, 6], hydraulic[:, 2], hydraulic[:, 3],
        hydraulic[:, 4], hydraulic[:, 7], relief[:, 3],
    ))
    limits = [500, 1e-5, 1e-5, 1e-5, 3e-3, 3e-3]
    assert np.all(np.max(np.abs(predicted - references), axis=0) < limits)
    assert np.max(predicted[:, 0]) > 5e6
