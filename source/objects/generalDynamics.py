"""Generalized-coordinate dynamics for supported linear-hydrodynamic devices.

The mechanical layout is supplied through body kinematics. Each body maps
generalized coordinates into its six WEC-Sim rigid-body coordinates, while
the engine assembles rigid inertia, cross-body added mass, restoring,
excitation, radiation, and linear PTO forces. Constant-frequency radiation
uses RK4; an impulse-response kernel uses an implicit trapezoidal step.
"""

from dataclasses import dataclass
from typing import Callable

import numpy as np


@dataclass(frozen=True)
class BodyMotion:
    """Six-DOF displacement and its generalized-coordinate derivatives."""

    displacement: np.ndarray
    jacobian: np.ndarray
    bias_acceleration: np.ndarray  # J-dot(q, v) times v


@dataclass(frozen=True)
class DynamicBody:
    rigid_mass: np.ndarray
    added_mass: tuple[np.ndarray, ...]
    damping: tuple[np.ndarray, ...]
    restoring: np.ndarray
    static_force: np.ndarray
    reference_position: np.ndarray
    motion: Callable[[np.ndarray, np.ndarray], BodyMotion]
    excitation: Callable[[float], np.ndarray]
    radiation_kernel: np.ndarray | None = None  # (lag, 6, 6 * body_count)


@dataclass(frozen=True)
class DynamicsResponse:
    time: np.ndarray
    coordinate: np.ndarray
    speed: np.ndarray
    acceleration: np.ndarray
    body_position: np.ndarray
    body_velocity: np.ndarray


class GeneralizedDynamics:
    """Integrate a constrained device described by independent coordinates.

    Body kinematics encode the constraints: the engine never integrates an
    unconstrained coordinate and then projects it back onto a joint. This
    avoids artificial constraint drift and accepts new layouts without
    changing the hydrodynamic force assembly.
    """

    def __init__(
        self,
        bodies: tuple[DynamicBody, ...],
        coordinate_count: int,
        *,
        pto_stiffness: np.ndarray | None = None,
        pto_damping: np.ndarray | None = None,
        pto_equilibrium: np.ndarray | None = None,
        pto_bias: np.ndarray | None = None,
    ):
        if not bodies or coordinate_count < 1:
            raise ValueError("a device needs bodies and independent coordinates")
        self.bodies = tuple(bodies)
        self.coordinate_count = coordinate_count
        n = coordinate_count
        self.pto_stiffness = (
            np.zeros((n, n)) if pto_stiffness is None
            else np.asarray(pto_stiffness, dtype=float)
        )
        self.pto_damping = (
            np.zeros((n, n)) if pto_damping is None
            else np.asarray(pto_damping, dtype=float)
        )
        self.pto_equilibrium = (
            np.zeros(n) if pto_equilibrium is None
            else np.asarray(pto_equilibrium, dtype=float)
        )
        self.pto_bias = (
            np.zeros(n) if pto_bias is None else np.asarray(pto_bias, dtype=float)
        )
        if (self.pto_stiffness.shape != (n, n)
                or self.pto_damping.shape != (n, n)
                or self.pto_equilibrium.shape != (n,)
                or self.pto_bias.shape != (n,)
                or not np.isfinite(self.pto_stiffness).all()
                or not np.isfinite(self.pto_damping).all()
                or not np.isfinite(self.pto_equilibrium).all()
                or not np.isfinite(self.pto_bias).all()):
            raise ValueError("PTO matrices and force vectors must match the coordinate count")
        for body in self.bodies:
            if (len(body.added_mass) != len(bodies)
                    or len(body.damping) != len(bodies)):
                raise ValueError("each body needs one hydrodynamic block per body")
            for matrix in (body.rigid_mass, body.restoring,
                           *body.added_mass, *body.damping):
                if np.shape(matrix) != (6, 6) or not np.isfinite(matrix).all():
                    raise ValueError("body mass, restoring, and radiation blocks must be finite 6x6 matrices")
            for vector in (body.static_force, body.reference_position):
                if np.shape(vector) != (6,) or not np.isfinite(vector).all():
                    raise ValueError("body static force and reference position must be finite six-vectors")
            kernel = body.radiation_kernel
            if kernel is not None and (
                np.ndim(kernel) != 3 or np.shape(kernel)[0] < 1
                or np.shape(kernel)[1:] != (6, 6 * len(bodies))
                or not np.isfinite(kernel).all()
            ):
                raise ValueError("radiation kernel must have shape (lags, 6, 6 * bodies)")

    def acceleration(
        self,
        at_time: float,
        coordinate: np.ndarray,
        speed: np.ndarray,
        *,
        known_radiation: tuple[np.ndarray, ...] | None = None,
        dt: float | None = None,
    ) -> np.ndarray:
        """Assemble M(q) and generalized force, then solve M(q) q'' = F."""
        if any(body.radiation_kernel is not None for body in self.bodies):
            if known_radiation is None or dt is None or dt <= 0:
                raise ValueError("radiation-memory acceleration needs history and dt")
        if known_radiation is not None and len(known_radiation) != len(self.bodies):
            raise ValueError("radiation history needs one force vector per body")
        n = self.coordinate_count
        motions = tuple(body.motion(coordinate, speed) for body in self.bodies)
        mass = np.zeros((n, n))
        force = (-self.pto_stiffness @ (coordinate - self.pto_equilibrium)
                 - self.pto_damping @ speed + self.pto_bias)
        velocities = tuple(motion.jacobian @ speed for motion in motions)
        for i, body in enumerate(self.bodies):
            motion = motions[i]
            j = motion.jacobian
            body_force = (body.static_force + body.excitation(at_time)
                          - body.restoring @ motion.displacement
                          - body.rigid_mass @ motion.bias_acceleration)
            mass += j.T @ body.rigid_mass @ j
            for k, other in enumerate(motions):
                a = body.added_mass[k]
                mass += j.T @ a @ other.jacobian
                body_force -= a @ other.bias_acceleration
                body_force -= body.damping[k] @ velocities[k]
            if known_radiation is not None:
                body_force -= known_radiation[i]
                if body.radiation_kernel is not None:
                    current_velocity = np.concatenate(velocities)
                    body_force -= dt / 2 * (body.radiation_kernel[0] @ current_velocity)
            force += j.T @ body_force
        try:
            answer = np.linalg.solve(mass, force)
        except np.linalg.LinAlgError as exc:
            raise ValueError("the constrained device has a singular mass matrix") from exc
        if not np.isfinite(answer).all():
            raise ValueError("the dynamics produced nonfinite acceleration")
        return answer

    def integrate(
        self,
        *,
        dt: float,
        end_time: float,
        initial_coordinate: np.ndarray | None = None,
        initial_speed: np.ndarray | None = None,
    ) -> DynamicsResponse:
        """Integrate using the radiation representation supplied by the bodies."""
        if not np.isfinite([dt, end_time]).all() or dt <= 0 or end_time < 0:
            raise ValueError("dt must be positive and end_time nonnegative")
        steps = round(end_time / dt)
        if not np.isclose(steps * dt, end_time, rtol=0, atol=1e-10):
            raise ValueError("end_time must be an integer multiple of dt")
        n = self.coordinate_count
        q = np.zeros((steps + 1, n))
        v = np.zeros_like(q)
        a = np.zeros_like(q)
        if initial_coordinate is not None:
            q[0] = np.asarray(initial_coordinate, dtype=float)
        if initial_speed is not None:
            v[0] = np.asarray(initial_speed, dtype=float)
        if not np.isfinite(q[0]).all() or not np.isfinite(v[0]).all():
            raise ValueError("initial state must be finite")
        time = np.arange(steps + 1) * dt
        memory = any(body.radiation_kernel is not None for body in self.bodies)
        if memory:
            self._integrate_memory(time, q, v, a, dt)
        else:
            self._integrate_regular(time, q, v, a, dt)
        positions = np.zeros((steps + 1, len(self.bodies), 6))
        velocities = np.zeros_like(positions)
        for step in range(steps + 1):
            for i, body in enumerate(self.bodies):
                motion = body.motion(q[step], v[step])
                positions[step, i] = body.reference_position + motion.displacement
                velocities[step, i] = motion.jacobian @ v[step]
        return DynamicsResponse(time, q, v, a, positions, velocities)

    def _integrate_regular(self, time, q, v, a, dt):
        n = self.coordinate_count

        def derivative(at_time, state):
            coordinate, speed = state[:n], state[n:]
            return np.concatenate((
                speed, self.acceleration(at_time, coordinate, speed),
            ))

        for step in range(len(time) - 1):
            t = time[step]
            state = np.concatenate((q[step], v[step]))
            k1 = derivative(t, state)
            k2 = derivative(t + dt / 2, state + dt * k1 / 2)
            k3 = derivative(t + dt / 2, state + dt * k2 / 2)
            k4 = derivative(t + dt, state + dt * k3)
            next_state = state + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            q[step + 1], v[step + 1] = next_state[:n], next_state[n:]
        for step, t in enumerate(time):
            a[step] = self.acceleration(t, q[step], v[step])

    def _integrate_memory(self, time, q, v, a, dt):
        count = len(self.bodies)
        kernels = [body.radiation_kernel for body in self.bodies]
        if any(kernel is None for kernel in kernels):
            raise ValueError("all bodies must use the same radiation representation")
        velocity_history = np.zeros((len(time), count * 6))

        def known_radiation(step):
            result = []
            for kernel in kernels:
                memory = min(step, len(kernel) - 1)
                known = dt * np.einsum(
                    "tij,tj->i", kernel[1:memory + 1],
                    velocity_history[step - memory:step][::-1],
                )
                result.append(known)
            return tuple(result)

        zeros = tuple(np.zeros(6) for _ in self.bodies)
        a[0] = self.acceleration(time[0], q[0], v[0],
                                 known_radiation=zeros, dt=dt)
        velocity_history[0] = np.concatenate([
            body.motion(q[0], v[0]).jacobian @ v[0] for body in self.bodies
        ])
        for step in range(1, len(time)):
            known = known_radiation(step)
            trial_speed = v[step - 1].copy()
            for _ in range(12):
                trial_coordinate = q[step - 1] + dt * (v[step - 1] + trial_speed) / 2
                trial_acceleration = self.acceleration(
                    time[step], trial_coordinate, trial_speed,
                    known_radiation=known, dt=dt,
                )
                next_speed = v[step - 1] + dt * (
                    a[step - 1] + trial_acceleration
                ) / 2
                if np.max(np.abs(next_speed - trial_speed)) < 1e-12:
                    trial_speed = next_speed
                    break
                trial_speed = next_speed
            else:
                raise RuntimeError("radiation-memory step did not converge; reduce dt")
            v[step] = trial_speed
            q[step] = q[step - 1] + dt * (v[step - 1] + v[step]) / 2
            a[step] = self.acceleration(
                time[step], q[step], v[step], known_radiation=known, dt=dt,
            )
            velocity_history[step] = np.concatenate([
                body.motion(q[step], v[step]).jacobian @ v[step]
                for body in self.bodies
            ])
