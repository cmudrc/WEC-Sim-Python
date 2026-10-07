"""Compare a modified official Sphere case with body-local PTO configuration."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularWave, WEC, WorldPoint


SPHERE_H5 = os.environ.get("WEC_SIM_SPHERE_H5")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (SPHERE_H5 and REFERENCE),
    reason="paired MATLAB configured Sphere PTO output not provided",
)


def test_configured_sphere_pto_against_matlab():
    reference = Path(REFERENCE)
    body = np.loadtxt(
        reference / "Sphere_PTO_Config_Configured_body1.csv", delimiter=",",
    )
    pto = np.loadtxt(
        reference / "Sphere_PTO_Config_Configured_pto1.csv", delimiter=",",
    )
    controller = np.loadtxt(
        reference / "Sphere_PTO_Config_controller.csv", delimiter=",",
    )
    wec = WEC("Configured Sphere PTO")
    sphere = wec.body(
        "sphere", Path(SPHERE_H5).resolve(),
        inertia=(20_907_301, 21_306_090.66, 37_085_481.11),
    )
    wec.coordinate("heave", sphere.move("heave"))
    wec.pto("main", sphere.at(1, 0, 2), WorldPoint(1, 0, 0),
            axis=(0, 0, 1), stiffness=50_000, damping=960_870)
    response = wec.run(
        RegularWave(2.5, 9.6664), dt=0.02, end_time=400, ramp_time=10,
    )
    np.testing.assert_allclose(response.time, body[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(response.time, pto[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(response.time, controller[:, 0], rtol=0, atol=1e-10)
    position = response.bodies["sphere"].position[:, 2]
    velocity = response.bodies["sphere"].velocity[:, 2]
    for name, actual, expected in (
        ("position", position, body[:, 3]),
        ("velocity", velocity, body[:, 9]),
        ("PTO force", response.ptos["main"].force, pto[:, 15]),
        ("controller force", response.ptos["main"].force, controller[:, 3]),
        ("combined force", response.ptos["main"].force, pto[:, 15] + controller[:, 3]),
        ("PTO power", response.ptos["main"].force * response.ptos["main"].velocity,
         pto[:, 21] + controller[:, 9]),
    ):
        difference = np.max(np.abs(actual - expected))
        opposite = np.max(np.abs(-actual - expected))
        print(f"{name}: max error {difference:.6g}, opposite sign {opposite:.6g}")
    assert np.max(np.abs(position - body[:, 3])) < 0.1
    assert np.max(np.abs(velocity - body[:, 9])) < 0.1
