"""Pair the published Sphere latching controller with MATLAB dynamics."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import LatchingControl, RegularWave, WEC, WorldPoint


SPHERE_H5 = os.environ.get("WEC_SIM_SPHERE_H5")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (SPHERE_H5 and REFERENCE),
    reason="paired MATLAB Sphere latching output not provided",
)


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = np.max(np.abs(actual - expected))
    assert error < limit, f"{label}: max error {error:.6g} exceeds {limit}"


def _latched(force, velocity):
    normal = np.abs(force + 49_181 * velocity)
    latched = np.abs(force + 37_308_296 * velocity)
    return latched < normal


def _switches(mode):
    starts = np.flatnonzero(~mode[:-1] & mode[1:]) + 1
    ends = np.flatnonzero(mode[:-1] & ~mode[1:]) + 1
    return starts, ends


def test_published_sphere_latching_against_matlab():
    reference = Path(REFERENCE)
    body = np.loadtxt(
        reference / "Sphere_Latching_Latching_body1.csv", delimiter=",",
    )
    controller = np.loadtxt(
        reference / "Sphere_Latching_controller.csv", delimiter=",",
    )
    assert body.shape == (40_001, 25)
    assert controller.shape == (40_001, 13)

    wec = WEC("Sphere latching")
    sphere = wec.body(
        "sphere", Path(SPHERE_H5).resolve(),
        inertia=(20_907_301, 21_306_090.66, 37_085_481.11),
    )
    wec.coordinate("heave", sphere.move("heave"))
    wec.pto(
        "latch", sphere.at(0, 0, 2), WorldPoint(0, 0, 0),
        axis=(0, 0, 1),
        control=LatchingControl(
            gain=49_181, latch_damping=37_308_296, latch_time=2.4,
        ),
    )
    response = wec.run(
        RegularWave(2.5, 9.6664), dt=0.01, end_time=400, ramp_time=100,
    )
    np.testing.assert_allclose(response.time, body[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(response.time, controller[:, 0],
                               rtol=0, atol=1e-10)

    position = response.bodies["sphere"].position[:, 2]
    velocity = response.bodies["sphere"].velocity[:, 2]
    force = -response.ptos["latch"].force
    power = -response.ptos["latch"].absorbed_power
    matlab_velocity = body[:, 9]
    matlab_force = controller[:, 3]
    matlab_power = controller[:, 9]
    matlab_mode = _latched(matlab_force, matlab_velocity)
    python_mode = _latched(force, velocity)

    # Both force laws dissipate energy. A large finite damping represents the
    # latch in the published Simulink model; it is not a rigid lock.
    _max_error(matlab_force,
               -np.where(matlab_mode, 37_308_296, 49_181) * matlab_velocity,
               1e-6, "MATLAB latching force law")
    _max_error(force, -np.where(python_mode, 37_308_296, 49_181) * velocity,
               1e-8, "Python latching force law")
    _max_error(matlab_power, matlab_force * matlab_velocity,
               1e-6, "MATLAB controller power law")
    _max_error(power, force * velocity, 1e-8,
               "Python controller power law")
    assert np.max(matlab_power) < 1e-6
    assert np.min(response.ptos["latch"].absorbed_power) >= -1e-8
    _max_error(response.raw.pto_generalized_force[:, 0], force, 1e-8,
               "applied generalized controller force")

    matlab_starts, matlab_ends = _switches(matlab_mode)
    python_starts, python_ends = _switches(python_mode)
    assert len(matlab_starts) == len(python_starts) == 83
    assert len(matlab_ends) == len(python_ends) == 82  # final latch continues
    _max_error(response.time[python_starts], response.time[matlab_starts],
               0.06000001, "latch event time")
    for starts, ends in ((matlab_starts, matlab_ends),
                         (python_starts, python_ends)):
        # The source's repeated floating-point timer addition yields 241
        # latched samples per complete 2.4 s interval at dt = 0.01 s.
        _max_error(response.time[ends] - response.time[starts[:len(ends)]],
                   np.full(len(ends), 2.41), 1e-10, "latch duration")
    mismatch = np.count_nonzero(matlab_mode != python_mode)
    assert mismatch <= 600, f"{mismatch} mismatched controller samples"

    # Switching is sensitive to a one-step velocity reversal. Gate the
    # trajectory, envelope, event times, and integrated energy separately.
    _max_error(position, body[:, 3], 0.42, "heave position")
    _max_error(velocity, matlab_velocity, 0.65, "heave velocity")
    _max_error(np.array([position.min(), position.max()]),
               np.array([body[:, 3].min(), body[:, 3].max()]),
               0.05, "heave envelope")
    python_energy = np.trapezoid(
        response.ptos["latch"].absorbed_power, response.time,
    )
    matlab_energy = -np.trapezoid(matlab_power, response.time)
    assert abs(python_energy - matlab_energy) < 3_000_000, (
        f"absorbed energy differs by {python_energy - matlab_energy:.1f} J"
    )
