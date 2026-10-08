"""Regular-wave floating body with rigid surge/heave/pitch and flexible modes.

The generalized mass is solved implicitly. MATLAB's sampled acceleration
feedback is a source-numerics diagnostic, not the default physical model.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .bodyClass import BodyClass


@dataclass(frozen=True)
class FloatingGBMResponse:
    time: np.ndarray
    body_position: np.ndarray
    body_velocity: np.ndarray
    mode_position: np.ndarray
    mode_velocity: np.ndarray
    mode_acceleration: np.ndarray
    wave_elevation: np.ndarray


def solve_floating_gbm_regular(
    hydro_file: str | Path, *, dt: float, end_time: float,
    height: float, period: float, pitch_inertia: float,
    mass: str | float = "equilibrium", ramp_time: float = 100,
    rho: float = 1000, g: float = 9.81,
) -> FloatingGBMResponse:
    """Solve one body with a 3-DOF floating joint and HDF5 flexible modes.

    The joint is at the body's origin; the equilibrium center of gravity must
    coincide with it. No PTO, mooring, drag, or body-to-body coupling is used.
    """
    body = BodyClass(str(hydro_file))
    body.bodyNumber = body.bodyTotal = 1
    body.readH5file()
    dof = int(np.asarray(body.dof).item())
    if dof <= 6 or int(np.asarray(body.dof_gbm).item()) != dof - 6:
        raise ValueError("floating GBM needs an HDF5 body with flexible modes")
    if not np.allclose(np.asarray(body.cg).ravel(), 0, atol=1e-10):
        raise ValueError("floating GBM currently needs the body CG at the joint origin")
    if mass == "equilibrium":
        physical_mass = rho * float(np.asarray(body.dispVol).item())
    else:
        physical_mass = float(mass)
    if physical_mass <= 0 or pitch_inertia <= 0:
        raise ValueError("body mass and pitch inertia must be positive")
    body.mass = mass
    omega = 2 * np.pi / period
    body.hydroForcePre(
        omega, [0], 1, np.array([0.0]), [], dt, rho, g,
        "regular", np.zeros((2, 1)), 1, 1, 0, 0, 0,
    )
    force = body.hydroForce
    active = np.r_[0, 2, 4, np.arange(6, dof)]
    n = len(active)
    mass_matrix = np.asarray(force["fAddedMass"])[np.ix_(active, active)].copy()
    mass_matrix[:3, :3] += np.diag([physical_mass, physical_mass, pitch_inertia])
    mass_matrix[3:, 3:] += force["gbm"]["mass_ff"]
    damping = np.asarray(force["fDamping"])[np.ix_(active, active)].copy()
    damping[3:, 3:] += force["gbm"]["damping"]
    stiffness = np.asarray(force["linearHydroRestCoef"])[np.ix_(active, active)].copy()
    stiffness[3:, 3:] += force["gbm"]["stiffness"]
    real = np.asarray(force["fExt"]["re"])[active]
    imaginary = np.asarray(force["fExt"]["im"])[active]
    inverse_mass = np.linalg.solve(mass_matrix, np.eye(n))

    def excitation(at_time):
        ramp = (1 if ramp_time == 0 or at_time >= ramp_time else
                (1 - np.cos(np.pi * at_time / ramp_time)) / 2)
        return height / 2 * ramp * (
            real * np.cos(omega * at_time)
            - imaginary * np.sin(omega * at_time)
        )

    def derivative(at_time, state):
        q, speed = state[:n], state[n:]
        acceleration = inverse_mass @ (
            excitation(at_time) - damping @ speed - stiffness @ q
        )
        return np.r_[speed, acceleration]

    time = np.arange(round(end_time / dt) + 1) * dt
    if not np.isclose(time[-1], end_time, atol=1e-10):
        raise ValueError("end_time must be an integer multiple of dt")
    state = np.zeros((len(time), 2 * n))
    for index in range(len(time) - 1):
        start = time[index]
        value = state[index]
        k1 = derivative(start, value)
        k2 = derivative(start + dt / 2, value + dt * k1 / 2)
        k3 = derivative(start + dt / 2, value + dt * k2 / 2)
        k4 = derivative(start + dt, value + dt * k3)
        state[index + 1] = value + dt * (k1 + 2*k2 + 2*k3 + k4) / 6
    acceleration = np.array([
        derivative(at_time, value)[n:]
        for at_time, value in zip(time, state)
    ])
    body_position = np.zeros((len(time), 6))
    body_velocity = np.zeros_like(body_position)
    body_position[:, [0, 2, 4]] = state[:, :3]
    body_velocity[:, [0, 2, 4]] = state[:, n:n+3]
    ramp = np.ones_like(time)
    if ramp_time:
        early = time < ramp_time
        ramp[early] = (1 - np.cos(np.pi * time[early] / ramp_time)) / 2
    return FloatingGBMResponse(
        time, body_position, body_velocity,
        state[:, 3:n], state[:, n+3:], acceleration[:, 3:],
        height / 2 * ramp * np.cos(omega * time),
    )
