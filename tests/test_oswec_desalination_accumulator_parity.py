"""Published Desalination accumulator pressure from its measured inlet flow."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import brentq

from wecsim import GasChargedAccumulator


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="pinned OSWEC Desalination output absent",
)


def test_published_accumulator_pressure_from_flow():
    measured = np.loadtxt(
        Path(REFERENCE) / "OSWEC_DESALINATION_SOURCE_simout1.csv",
        delimiter=",",
    )
    assert measured.shape == (30_001, 8)
    time, inlet_pressure, accumulator_flow = (
        measured[:, 0], measured[:, 6], measured[:, 7]
    )
    np.testing.assert_allclose(np.diff(time), .01, rtol=0, atol=1e-10)
    accumulator = GasChargedAccumulator(
        initial_gas_volume=4,
        precharge_pressure=30e5,
        exponent=1.4,
        atmospheric_pressure=101_325,
        dead_gas_volume=4e-5,
        hard_stop_stiffness=1e7,
        hard_stop_damping=1e7,
    )
    # The source starts its connected network near zero pressure. Its finite
    # lower hard stop therefore permits a small negative liquid volume even
    # though the block requests a zero-volume initialization target.
    initial_volume = brentq(
        lambda volume: accumulator.pressure(
            volume, accumulator_flow[0]) - inlet_pressure[0],
        -1, 0,
    )
    assert -.3 < initial_volume < -.2
    liquid_volume = initial_volume + np.r_[
        0, np.cumsum(np.diff(time) * accumulator_flow[1:]),
    ]
    actual = accumulator.pressure(liquid_volume, accumulator_flow)
    np.testing.assert_allclose(actual, inlet_pressure, rtol=0, atol=1e-3)
    assert np.max(inlet_pressure) > 5e6
