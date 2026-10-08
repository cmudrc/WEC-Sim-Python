"""Paired force law for the published MBARI cable application."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import PlanarCableAttachment, WecSimCableTension


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="pinned MBARI cable MATLAB output absent",
)


def _load(name):
    return np.loadtxt(Path(REFERENCE) / f"CABLE_SOURCE_{name}.csv",
                      delimiter=",", ndmin=2)


def test_published_mbari_cable_force_from_relative_motion():
    parameters = _load("parameters").ravel()
    np.testing.assert_allclose(parameters, (1e6, 100, 17.8, 18),
                               rtol=0, atol=1e-10)
    cable = _load("cable1")
    wave = _load("wave")
    bodies = [_load(f"Cable_body{number}") for number in (1, 2, 3)]
    assert cable.shape == (12_001, 37)
    assert wave.shape == (12_001, 2)
    for body in bodies:
        assert body.shape == (12_001, 25)
        np.testing.assert_allclose(body[:, 0], cable[:, 0],
                                   rtol=0, atol=1e-10)

    time = cable[:, 0]
    ramp = np.where(time >= 40, 1,
                    (1 + np.cos(np.pi + np.pi * time / 40)) / 2)
    expected_wave = .5 * ramp * np.cos(2 * np.pi * time / 8)
    np.testing.assert_allclose(wave[:, 0], time, rtol=0, atol=1e-10)
    np.testing.assert_allclose(wave[:, 1], expected_wave,
                               rtol=0, atol=1e-11)

    law = WecSimCableTension(*parameters)
    force = law.force_z(cable[:, 3], cable[:, 9])
    source_force = cable[:, 21]
    assert np.max(np.abs(source_force)) > 1e5
    np.testing.assert_allclose(force, source_force, rtol=0, atol=1e-5)
    assert np.all(force[np.abs(cable[:, 3] + parameters[3])
                        <= parameters[2]] == 0)


def test_published_cable_attachment_from_body_local_points():
    cable = _load("cable1")
    base = _load("Cable_body2")
    follower = _load("Cable_body3")
    attachment = PlanarCableAttachment(
        base_offset=(0, 1.95), follower_offset=(0, -5.2),
        initial_length=18,
    )

    def pose_and_rate(body):
        return (body[:, (1, 3, 5)], body[:, (7, 9, 11)])

    displacement, speed = attachment.motion(
        *pose_and_rate(base), *pose_and_rate(follower),
    )
    np.testing.assert_allclose(displacement, cable[:, 3],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(speed, cable[:, 9],
                               rtol=0, atol=1e-11)
    force = WecSimCableTension(1e6, 100, 17.8, 18).force_z(
        displacement, speed,
    )
    np.testing.assert_allclose(force, cable[:, 21], rtol=0, atol=1e-5)
