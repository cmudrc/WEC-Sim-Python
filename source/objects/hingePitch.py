"""Forced pitch response of one hydrodynamic body on a horizontal hinge.

The solver accepts a six-component excitation force history in the body's
equilibrium global axes. Its hydrodynamic coefficients come from the production
``BodyClass`` preprocessing. A separate wave model must supply excitation;
this module does not synthesize irregular waves.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .bodyClass import BodyClass


@dataclass(frozen=True)
class HingePitchResponse:
    time: np.ndarray
    angle: np.ndarray
    angular_velocity: np.ndarray
    center_position: np.ndarray
    center_velocity: np.ndarray
    excitation_torque: np.ndarray
    pto_torque: np.ndarray


def solve_hinged_pitch_from_excitation(
    h5_file: str | Path,
    excitation_force: np.ndarray,
    *,
    hinge_z: float,
    body_mass: float,
    pitch_inertia: float,
    pto_damping: float,
    pto_stiffness: float = 0.0,
    dt: float = 0.1,
    memory_time: float = 30.0,
    rho: float = 1000.0,
    g: float = 9.81,
) -> HingePitchResponse:
    """Integrate hinge pitch with hydrodynamic radiation and a linear PTO.

    ``excitation_force`` is an ``(N, 6)`` array with columns surge, sway,
    heave, roll, pitch, yaw at uniform ``dt``. The body starts at its HDF5
    equilibrium position, zero angle, and zero velocity. The hinge is parallel
    to the global y axis. The model includes nonlinear rigid-body kinematics,
    hydrostatic restoring, infinite-frequency added mass, and the radiation
    impulse-response convolution. It assumes a fixed hinge and no body-to-body
    hydrodynamic interaction, Morison forces, or other DOFs.

    The fixed-step trapezoidal update uses a finite memory window. It matches
    the published OSWEC example's 0.1 s step and 30 s radiation memory when
    supplied with that run's excitation forces.
    """
    force = np.asarray(excitation_force, dtype=float)
    if force.ndim != 2 or force.shape[1] != 6 or len(force) < 2:
        raise ValueError("excitation_force must have shape (N, 6), N >= 2")
    if not np.isfinite(force).all():
        raise ValueError("excitation_force must be finite")
    values = [hinge_z, body_mass, pitch_inertia, pto_damping,
              pto_stiffness, dt, memory_time, rho, g]
    if not np.isfinite(values).all():
        raise ValueError("solver parameters must be finite")
    if (body_mass <= 0 or pitch_inertia <= 0 or pto_damping < 0
            or pto_stiffness < 0 or dt <= 0 or memory_time <= 0
            or rho <= 0 or g <= 0):
        raise ValueError("mass, inertia, dt, memory time, rho, and g must be positive")
    memory_steps = round(memory_time / dt)
    if not np.isclose(memory_steps * dt, memory_time, rtol=0, atol=1e-10):
        raise ValueError("memory_time must be an integer multiple of dt")

    body = BodyClass(str(h5_file))
    body.bodyNumber = 1
    body.readH5file()
    if int(np.asarray(body.dof).item()) != 6:
        raise ValueError("the hinged-pitch solver requires one six-DOF body")
    cg_z = float(body.hydroData["properties"]["cg"][0, 2])
    lever = cg_z - hinge_z
    if lever <= 0:
        raise ValueError("the hinge must be below the center of gravity")
    irf_time = body.hydroData["hydro_coeffs"]["radiation_damping"]["impulse_response_fun"]["t"]
    if memory_time > np.max(irf_time) + 1e-10:
        raise ValueError("memory_time exceeds the radiation kernel in the HDF5 file")

    count = len(force)
    time = np.arange(count) * dt
    convolution_time = np.arange(memory_steps + 1) * dt
    body.mass = body_mass
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {
        "Drag": np.zeros((6, 6)),
        "cd": np.zeros(6),
        "characteristicArea": np.zeros(6),
    }
    body.linearDamping = np.zeros((6, 6))
    body.hydroForcePre(
        [], [0], len(convolution_time), convolution_time, [], dt, rho, g,
        "noWaveCIC", np.vstack((time, np.zeros_like(time))),
        1, 1, 0, 0, 0,
    )

    added_mass = np.asarray(body.hydroForce["fAddedMass"])
    restoring = np.asarray(body.hydroForce["linearHydroRestCoef"])
    kernel = np.asarray(body.hydroForce["irkb"])
    displaced_volume = float(np.asarray(body.dispVol).item())
    static_force = np.array([
        0.0, 0.0, (rho * displaced_volume - body_mass) * g, 0.0, 0.0, 0.0,
    ])
    rigid_inertia = body_mass * lever**2 + pitch_inertia

    def jacobian(angle):
        return np.array([
            lever * np.cos(angle), 0.0, -lever * np.sin(angle),
            0.0, 1.0, 0.0,
        ])

    def jacobian_prime(angle):
        return np.array([
            -lever * np.sin(angle), 0.0, -lever * np.cos(angle),
            0.0, 0.0, 0.0,
        ])

    def acceleration(angle, speed, index, known_radiation):
        j = jacobian(angle)
        jp = jacobian_prime(angle)
        displacement = np.array([
            lever * np.sin(angle), 0.0,
            lever * (np.cos(angle) - 1.0), 0.0, angle, 0.0,
        ])
        mass = rigid_inertia + j @ added_mass @ j
        if mass <= 0 or not np.isfinite(mass):
            raise ValueError("effective hinge inertia must be positive and finite")
        hydrostatic_torque = j @ (static_force - restoring @ displacement)
        excitation_torque = j @ force[index]
        radiation_torque = j @ (
            known_radiation + dt / 2 * (kernel[0] @ j) * speed
        )
        added_mass_curvature = (j @ added_mass @ jp) * speed**2
        return (hydrostatic_torque + excitation_torque
                - pto_damping * speed - pto_stiffness * angle
                - radiation_torque - added_mass_curvature) / mass

    angle = np.zeros(count)
    speed = np.zeros(count)
    angular_acceleration = np.zeros(count)
    six_velocity = np.zeros((count, 6))
    angular_acceleration[0] = acceleration(0.0, 0.0, 0, np.zeros(6))
    for step in range(1, count):
        memory = min(step, memory_steps)
        known_radiation = dt * np.einsum(
            "tij,tj->i", kernel[1:memory + 1],
            six_velocity[step - memory:step][::-1],
        )
        trial_speed = speed[step - 1]
        for _ in range(12):
            trial_angle = angle[step - 1] + dt * (speed[step - 1] + trial_speed) / 2
            trial_acceleration = acceleration(
                trial_angle, trial_speed, step, known_radiation,
            )
            next_speed = speed[step - 1] + dt * (
                angular_acceleration[step - 1] + trial_acceleration
            ) / 2
            if abs(next_speed - trial_speed) < 1e-12:
                trial_speed = next_speed
                break
            trial_speed = next_speed
        else:
            raise RuntimeError("hinge-pitch step did not converge; reduce dt")
        speed[step] = trial_speed
        angle[step] = angle[step - 1] + dt * (
            speed[step - 1] + speed[step]
        ) / 2
        angular_acceleration[step] = acceleration(
            angle[step], speed[step], step, known_radiation,
        )
        six_velocity[step] = jacobian(angle[step]) * speed[step]

    center_position = np.column_stack((
        lever * np.sin(angle), np.zeros(count),
        hinge_z + lever * np.cos(angle),
    ))
    center_velocity = six_velocity[:, :3]
    excitation_torque = np.einsum(
        "ij,ij->i", np.stack([jacobian(a) for a in angle]), force,
    )
    return HingePitchResponse(
        time=time,
        angle=angle,
        angular_velocity=speed,
        center_position=center_position,
        center_velocity=center_velocity,
        excitation_torque=excitation_torque,
        pto_torque=-pto_damping * speed - pto_stiffness * angle,
    )
