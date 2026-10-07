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


def _max_error(actual, expected, limit, label):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, label
    assert np.isfinite(expected).all(), label
    error = np.max(np.abs(actual - expected))
    assert error < limit, f"{label}: max error {error:.6g} exceeds {limit}"


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
    stroke = response.ptos["main"].stroke
    speed = response.ptos["main"].velocity
    # MATLAB's PTO and controller force/stroke channels use the opposite
    # sign to Python's from-to axis. Their power channels have the same sign.
    native_force = -(50_000 * stroke + 100_000 * speed)
    controller_force = -860_870 * speed
    _max_error(position, body[:, 3], 1e-4, "sphere heave position")
    _max_error(velocity, body[:, 9], 1e-4, "sphere heave velocity")
    _max_error(-stroke, pto[:, 3], 1e-4, "PTO stroke")
    _max_error(-speed, pto[:, 9], 1e-4, "PTO speed")
    _max_error(-native_force, pto[:, 15], 15, "native PTO force")
    _max_error(-controller_force, controller[:, 3], 100, "controller force")
    _max_error(-response.ptos["main"].force,
               pto[:, 15] + controller[:, 3], 120, "combined force")
    _max_error(native_force * speed, pto[:, 21], 8, "native PTO power")
    _max_error(controller_force * speed, controller[:, 9], 40,
               "controller power")
    _max_error(response.ptos["main"].force * speed,
               pto[:, 21] + controller[:, 9], 50, "combined power")
