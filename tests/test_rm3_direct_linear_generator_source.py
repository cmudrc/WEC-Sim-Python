"""Validate the pinned RM3 direct linear-generator reference export."""

import os
from pathlib import Path

import numpy as np
import pytest


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="paired MATLAB RM3 direct-generator output absent",
)


def test_published_direct_generator_source_signals():
    folder = Path(REFERENCE)
    bodies = [
        np.loadtxt(folder / f"RM3_DD_PTO_RM3_DD_PTO_body{body}.csv",
                   delimiter=",")
        for body in (1, 2)
    ]
    pto = np.loadtxt(folder / "RM3_DD_PTO_RM3_DD_PTO_pto1.csv",
                     delimiter=",")
    drive = np.loadtxt(folder / "RM3_DD_PTO_drive.csv", delimiter=",")
    assert all(values.shape[0] == 40_001 for values in (*bodies, pto, drive))
    assert drive.shape == (40_001, 12)
    time = np.arange(40_001) * 0.01
    for values in (*bodies, pto, drive):
        assert np.isfinite(values).all()
        np.testing.assert_allclose(values[:, 0], time, rtol=0, atol=1e-8)

    power, force, friction = drive[:, 1:4].T
    current = drive[:, 4:7]
    voltage = drive[:, 7:10]
    velocity, electrical_power = drive[:, 10:12].T
    np.testing.assert_allclose(friction, -100 * velocity,
                               rtol=1e-9, atol=1e-6)
    np.testing.assert_allclose(power, -force * velocity,
                               rtol=1e-9, atol=1e-5)
    np.testing.assert_allclose(voltage, -117.6471 * current,
                               rtol=1e-9, atol=1e-6)
    np.testing.assert_allclose(electrical_power,
                               -np.sum(voltage * current, axis=1),
                               rtol=1e-9, atol=1e-4)
    np.testing.assert_allclose(pto[:, 3],
                               bodies[0][:, 3] - bodies[1][:, 3]
                               - (bodies[0][0, 3] - bodies[1][0, 3]),
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(pto[:, 9],
                               bodies[0][:, 9] - bodies[1][:, 9],
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(velocity, pto[:, 9], rtol=0, atol=1e-9)
    np.testing.assert_allclose(pto[:, 15], 0, rtol=0, atol=1e-7)
    np.testing.assert_allclose(pto[:, 21], 0, rtol=0, atol=1e-7)
