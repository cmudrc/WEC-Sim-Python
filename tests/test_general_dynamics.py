"""Check the shared dynamics stage against independent oscillator solutions."""

import numpy as np
from scipy.integrate import solve_ivp

from wecsim.generalDynamics import (
    BodyMotion, DynamicBody, GeneralizedDynamics,
)


def _oscillator(*, mass=2.0, stiffness=8.0, radiation=None,
                pto_damping=0.0, pto_stiffness=0.0, pto_equilibrium=0.0):
    jacobian = np.zeros((6, 1))
    jacobian[2, 0] = 1
    rigid_mass = np.zeros((6, 6))
    rigid_mass[2, 2] = mass
    restoring = np.zeros((6, 6))
    restoring[2, 2] = stiffness

    def motion(q, v):
        displacement = np.zeros(6)
        displacement[2] = q[0]
        return BodyMotion(displacement, jacobian, np.zeros(6))

    body = DynamicBody(
        rigid_mass=rigid_mass,
        added_mass=(np.zeros((6, 6)),),
        damping=(np.zeros((6, 6)),),
        restoring=restoring,
        static_force=np.zeros(6),
        reference_position=np.zeros(6),
        motion=motion,
        excitation=lambda t: np.zeros(6),
        radiation_kernel=radiation,
    )
    return GeneralizedDynamics(
        (body,), 1, pto_damping=np.array([[pto_damping]]),
        pto_stiffness=np.array([[pto_stiffness]]),
        pto_equilibrium=np.array([pto_equilibrium]),
    )


def test_regular_rk4_matches_damped_oscillator_solution():
    mass, stiffness, damping = 2.0, 8.0, 0.5
    system = _oscillator(mass=mass, stiffness=stiffness, pto_damping=damping)
    response = system.integrate(
        dt=0.01, end_time=10, initial_coordinate=np.array([1.0]),
    )
    decay = damping / (2 * mass)
    frequency = np.sqrt(stiffness / mass - decay**2)
    expected = np.exp(-decay * response.time) * (
        np.cos(frequency * response.time)
        + decay / frequency * np.sin(frequency * response.time)
    )
    assert np.max(np.abs(response.coordinate[:, 0] - expected)) < 1e-7
    assert np.max(np.abs(response.body_position[:, 0, 2] - expected)) < 1e-7


def test_pto_equilibrium_shift_changes_force_and_motion():
    system = _oscillator(stiffness=0, pto_stiffness=8,
                         pto_equilibrium=1)
    response = system.integrate(dt=0.01, end_time=2)
    expected = 1 - np.cos(2 * response.time)
    assert np.max(np.abs(response.coordinate[:, 0] - expected)) < 1e-8
    np.testing.assert_allclose(
        system.acceleration(0, np.array([0.0]), np.array([0.0])),
        [4.0], rtol=0, atol=1e-12,
    )


def test_radiation_memory_matches_augmented_state_equation():
    dt, duration, mass, stiffness = 0.01, 5.0, 2.0, 8.0
    strength, decay = 1.5, 0.8
    lag = np.arange(round(duration / dt) + 1) * dt
    kernel = np.zeros((len(lag), 6, 6))
    kernel[:, 2, 2] = strength * np.exp(-decay * lag)
    system = _oscillator(mass=mass, stiffness=stiffness, radiation=kernel)
    response = system.integrate(
        dt=dt, end_time=duration, initial_coordinate=np.array([1.0]),
    )

    def derivative(t, state):
        q, v, radiation = state
        return [v, (-stiffness * q - radiation) / mass,
                strength * v - decay * radiation]

    expected = solve_ivp(
        derivative, (0, duration), [1, 0, 0], t_eval=response.time,
        rtol=1e-11, atol=1e-13,
    )
    assert expected.success
    assert np.max(np.abs(response.coordinate[:, 0] - expected.y[0])) < 3e-4
    assert np.max(np.abs(response.speed[:, 0] - expected.y[1])) < 6e-4
