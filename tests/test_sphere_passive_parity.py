"""Compare the published Sphere passive-controller case with Python dynamics."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularWave, WEC, WorldPoint


SPHERE_H5 = os.environ.get("WEC_SIM_SPHERE_H5")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (SPHERE_H5 and REFERENCE),
    reason="paired MATLAB Sphere passive-controller baseline not provided",
)


def test_sphere_passive_controller_against_matlab():
    reference = Path(REFERENCE)
    body = np.loadtxt(reference / "Sphere_Passive_Passive_P_body1.csv", delimiter=",")
    controller = np.loadtxt(reference / "Sphere_Passive_controller.csv", delimiter=",")
    wec = WEC("Sphere passive controller")
    sphere = wec.body(
        "sphere", Path(SPHERE_H5).resolve(),
        inertia=(20_907_301, 21_306_090.66, 37_085_481.11),
    )
    wec.coordinate("heave", sphere.move("heave"))
    wec.pto("passive", sphere.at(0, 0, 2), WorldPoint(0, 0, 0),
            axis=(0, 0, 1), damping=860_870)
    response = wec.run(
        RegularWave(2.5, 9.6664), dt=0.02, end_time=400, ramp_time=10,
    )
    np.testing.assert_allclose(response.time, body[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(response.time, controller[:, 0], rtol=0, atol=1e-10)
    position = response.bodies["sphere"].position[:, 2]
    velocity = response.bodies["sphere"].velocity[:, 2]
    # The Python PTO reports force along its from-to axis and positive
    # absorbed power. MATLAB's controller channels use the opposite signs.
    force = -response.ptos["passive"].force
    power = -response.ptos["passive"].absorbed_power
    for name, actual, expected, limit in (
        ("position", position, body[:, 3], 1.5e-4),
        ("velocity", velocity, body[:, 9], 1.5e-4),
        ("controller force", force, controller[:, 3], 150),
        ("controller power", power, controller[:, 9], 100),
    ):
        assert actual.shape == expected.shape, name
        assert np.isfinite(expected).all(), name
        error = np.max(np.abs(actual - expected))
        assert error < limit, f"{name}: max error {error:.6g} exceeds {limit}"
