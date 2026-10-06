"""Linear, single-body heave free decay using WEC-Sim hydrodynamic data.

This solver covers the published Sphere ``noWaveCIC`` free-decay cases. It
uses the production BodyClass preprocessing and a fixed-step trapezoidal
integration of the Cummins radiation convolution. Other WEC-Sim degrees of
freedom, excitation, constraints, PTOs, and nonlinear forces are outside its
scope.

The solved equation is ``(m + A_inf) q'' + C q + integral(K(tau)
q'(t - tau) d tau) = 0`` in heave, with zero initial velocity.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .bodyClass import BodyClass


@dataclass(frozen=True)
class HeaveResponse:
    time: np.ndarray
    position: np.ndarray
    displacement: np.ndarray
    velocity: np.ndarray
    force_total: np.ndarray


@dataclass(frozen=True)
class TwoBodyHeaveResponse:
    """Heave-only response of two bodies joined by a linear PTO.

    Body-indexed arrays have shape ``(time_steps, 2)``. ``pto_force`` is the
    force exerted by the PTO on body 1, with the opposite force on body 2.
    """

    time: np.ndarray
    position: np.ndarray
    displacement: np.ndarray
    velocity: np.ndarray
    excitation_force: np.ndarray
    pto_force: np.ndarray


def _sample_count(duration, dt):
    count = round(duration / dt)
    if count < 0 or not np.isclose(count * dt, duration, rtol=0, atol=1e-10):
        raise ValueError("durations must be nonnegative multiples of dt")
    return count


def solve_two_body_regular_heave(
    h5_file: str | Path,
    *,
    wave_height: float,
    wave_period: float,
    pto_damping: float,
    pto_stiffness: float = 0.0,
    dt: float = 0.1,
    end_time: float = 400.0,
    ramp_time: float = 100.0,
    rho: float = 1000.0,
    g: float = 9.81,
    wave_direction: float = 0.0,
    body_masses: tuple[float | str, float | str] = ("equilibrium", "equilibrium"),
) -> TwoBodyHeaveResponse:
    """Solve the linear heave subsystem for a regular-wave, two-body WEC.

    This uses each body's frequency-dependent added mass and radiation
    damping at the incident frequency, hydrostatic heave stiffness, and
    complex excitation coefficients from ``BodyClass.hydroForcePre``. The
    relative-motion PTO supplies equal and opposite damping/stiffness forces.
    It assumes no body-to-body hydrodynamic interaction or other active DOFs.
    Integration is fixed-step classical RK4, matching the published RM3
    example's ``ode4`` setting. It does not represent surge, pitch, or the
    full Simscape constraint and PTO mechanics.
    """
    inputs = [wave_height, wave_period, pto_damping, pto_stiffness,
              dt, end_time, ramp_time, rho, g, wave_direction]
    if not np.isfinite(inputs).all():
        raise ValueError("solver inputs must be finite")
    if (wave_height < 0 or wave_period <= 0 or pto_damping < 0
            or pto_stiffness < 0 or dt <= 0 or end_time < 0
            or ramp_time < 0 or rho <= 0 or g <= 0):
        raise ValueError("wave period, dt, rho, and g must be positive; other inputs must be nonnegative")
    if len(body_masses) != 2:
        raise ValueError("body_masses must contain two values")

    steps = _sample_count(end_time, dt)
    time = np.arange(steps + 1) * dt
    omega = 2 * np.pi / wave_period
    bodies = []
    for index, mass_setting in enumerate(body_masses, start=1):
        body = BodyClass(str(h5_file))
        body.bodyNumber = index
        body.bodyTotal = 2
        body.readH5file()
        if int(np.asarray(body.dof).item()) != 6:
            raise ValueError("each hydrodynamic body must have six DOFs")
        if mass_setting != "equilibrium":
            mass_setting = float(mass_setting)
            if not np.isfinite(mass_setting) or mass_setting <= 0:
                raise ValueError("body masses must be positive and finite")
        body.mass = mass_setting
        body.hydroStiffness = np.zeros((6, 6))
        body.viscDrag = {
            "Drag": np.zeros((6, 6)),
            "cd": np.zeros(6),
            "characteristicArea": np.zeros(6),
        }
        body.linearDamping = np.zeros((6, 6))
        body.hydroForcePre(
            omega, [wave_direction], 1, np.array([0.0]), [], dt,
            rho, g, "regular", np.vstack((time, np.zeros_like(time))),
            index, 2, 0, 0, 0,
        )
        bodies.append(body)

    mass = np.array([
        float(np.asarray(body.mass).item())
        + body.hydroForce["fAddedMass"][2, 2]
        for body in bodies
    ])
    if not np.isfinite(mass).all() or np.any(mass <= 0):
        raise ValueError("body mass plus added mass must be positive and finite")
    stiffness = np.diag([
        body.hydroForce["linearHydroRestCoef"][2, 2] for body in bodies
    ])
    damping = np.diag([
        body.hydroForce["fDamping"][2, 2] for body in bodies
    ])
    relative = np.array([[1.0, -1.0], [-1.0, 1.0]])
    stiffness += pto_stiffness * relative
    damping += pto_damping * relative
    excitation_re = np.array([
        body.hydroForce["fExt"]["re"][2] for body in bodies
    ])
    excitation_im = np.array([
        body.hydroForce["fExt"]["im"][2] for body in bodies
    ])

    def excitation(at_time):
        ramp = (1.0 if ramp_time == 0 or at_time >= ramp_time
                else (1.0 - np.cos(np.pi * at_time / ramp_time)) / 2.0)
        return (wave_height / 2.0) * ramp * (
            excitation_re * np.cos(omega * at_time)
            - excitation_im * np.sin(omega * at_time)
        )

    def derivative(at_time, state):
        displacement = state[:2]
        velocity = state[2:]
        acceleration = (excitation(at_time) - stiffness @ displacement
                        - damping @ velocity) / mass
        return np.concatenate((velocity, acceleration))

    state = np.zeros((steps + 1, 4))
    for step in range(steps):
        t = time[step]
        y = state[step]
        k1 = derivative(t, y)
        k2 = derivative(t + dt / 2, y + dt * k1 / 2)
        k3 = derivative(t + dt / 2, y + dt * k2 / 2)
        k4 = derivative(t + dt, y + dt * k3)
        state[step + 1] = y + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6

    displacement = state[:, :2]
    velocity = state[:, 2:]
    equilibrium_z = np.array([
        float(body.hydroData["properties"]["cg"][0, 2])
        for body in bodies
    ])
    return TwoBodyHeaveResponse(
        time=time,
        position=equilibrium_z + displacement,
        displacement=displacement,
        velocity=velocity,
        excitation_force=np.stack([excitation(t) for t in time]),
        pto_force=(-pto_stiffness * (displacement[:, 0] - displacement[:, 1])
                   - pto_damping * (velocity[:, 0] - velocity[:, 1])),
    )


def solve_heave_free_decay(
    h5_file: str | Path,
    initial_displacement: float,
    *,
    dt: float = 0.01,
    end_time: float = 40.0,
    cic_end_time: float = 15.0,
    rho: float = 1000.0,
    g: float = 9.81,
) -> HeaveResponse:
    """Return a linear heave response relative to the HDF5 equilibrium center.

    The model assumes one unconstrained heave DOF, zero incident waves, zero
    initial velocity, and no PTO, mooring, viscous drag, or body interaction.
    ``position`` is the absolute vertical center of gravity in meters;
    ``displacement`` is relative to its equilibrium value.
    """
    if not np.isfinite([initial_displacement, dt, end_time, cic_end_time, rho, g]).all():
        raise ValueError("solver inputs must be finite")
    if dt <= 0 or cic_end_time <= 0 or rho <= 0 or g <= 0:
        raise ValueError("dt, cic_end_time, rho, and g must be positive")
    steps = _sample_count(end_time, dt)
    memory_steps = _sample_count(cic_end_time, dt)
    time = np.arange(steps + 1) * dt
    memory_time = np.arange(memory_steps + 1) * dt

    body = BodyClass(str(h5_file))
    body.bodyNumber = 1
    body.readH5file()
    if int(np.asarray(body.dof).item()) != 6:
        raise ValueError("the heave solver requires one six-DOF hydrodynamic body")
    irf_time = body.hydroData["hydro_coeffs"]["radiation_damping"]["impulse_response_fun"]["t"]
    if memory_time[-1] > np.max(irf_time) + 1e-10:
        raise ValueError("cic_end_time exceeds the radiation kernel in the HDF5 file")
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {
        "Drag": np.zeros((6, 6)),
        "cd": np.zeros(6),
        "characteristicArea": np.zeros(6),
    }
    body.linearDamping = np.zeros((6, 6))
    body.mass = "equilibrium"
    body.hydroForcePre(
        [], [0], len(memory_time), memory_time, [], dt, rho, g,
        "noWaveCIC", np.vstack((time, np.zeros_like(time))),
        1, [], 0, 0, 0,
    )

    stiffness = float(body.hydroForce["linearHydroRestCoef"][2, 2])
    mass = float(np.asarray(body.mass).item() + body.hydroForce["fAddedMass"][2, 2])
    kernel = np.asarray(body.hydroForce["irkb"][:, 2, 2])
    equilibrium_z = float(body.hydroData["properties"]["cg"][0, 2])
    if mass <= 0:
        raise ValueError("body mass plus infinite-frequency added mass must be positive")

    displacement = np.zeros(steps + 1)
    velocity = np.zeros(steps + 1)
    acceleration = np.zeros(steps + 1)
    displacement[0] = initial_displacement
    acceleration[0] = -stiffness * initial_displacement / mass
    denominator = 1 + dt / (2 * mass) * (stiffness * dt / 2 + dt * kernel[0] / 2)

    for step in range(1, steps + 1):
        memory = min(step, memory_steps)
        known_radiation = dt * np.dot(
            kernel[1:memory + 1], velocity[step - memory:step][::-1],
        )
        velocity[step] = (
            velocity[step - 1] + dt * acceleration[step - 1] / 2
            - dt / (2 * mass) * (
                stiffness * displacement[step - 1]
                + stiffness * dt * velocity[step - 1] / 2
                + known_radiation
            )
        ) / denominator
        displacement[step] = (
            displacement[step - 1]
            + dt * (velocity[step - 1] + velocity[step]) / 2
        )
        acceleration[step] = (
            -stiffness * displacement[step]
            - known_radiation
            - dt * kernel[0] * velocity[step] / 2
        ) / mass

    return HeaveResponse(
        time=time,
        position=equilibrium_z + displacement,
        displacement=displacement,
        velocity=velocity,
        force_total=mass * acceleration,
    )
