"""Hydrodynamic force laws for the published MOST VolturnUS platform."""

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
from scipy.interpolate import CubicSpline

from .generalDynamics import BodyMotion, DynamicBody, GeneralizedDynamics
from .mostMooring import MostStaticMooring


@dataclass(frozen=True)
class MostPlatformResponse:
    time: np.ndarray
    position: np.ndarray  # world-frame center-of-gravity pose, m and rad
    velocity: np.ndarray  # world-frame center-of-gravity velocity, m/s and rad/s
    acceleration: np.ndarray  # independent surge, heave, pitch acceleration


@dataclass(frozen=True)
class MostPlatformHydrodynamics:
    """Evaluate source-case restoring, drag, and radiation on platform motion.

    Forces follow WEC-Sim's logged resisting-force convention. ``pose`` is
    the six-component position of the platform center of gravity in world
    coordinates. ``simulate`` advances the dominant surge, heave, and pitch
    coordinates when a tower-base reaction history is supplied.
    """

    mass: float
    inertia: np.ndarray
    equilibrium_pose: np.ndarray
    restoring_stiffness: np.ndarray
    static_restoring: np.ndarray
    drag: np.ndarray
    radiation_time: np.ndarray
    radiation_kernel: np.ndarray
    added_mass: np.ndarray
    tower_base_height: float = 15.0  # world height at zero platform motion

    @classmethod
    def from_volturnus(cls, h5_file: str | Path,
                       mass_properties_file: str | Path) -> "MostPlatformHydrodynamics":
        """Read the pinned MOST hydrodynamics and platform mass properties."""
        rho, gravity = 1025.0, 9.80665
        with h5py.File(h5_file) as hydro:
            body = hydro["body1"]
            volume = float(body["properties/disp_vol"][0, 0])
            center_gravity = body["properties/cg"][:].ravel()
            center_buoyancy = body["properties/cb"][:].ravel()
            coefficients = body["hydro_coeffs"]
            stiffness = coefficients["linear_restoring_stiffness"][:] * rho * gravity
            added_mass = coefficients["added_mass/inf_freq"][:] * rho
            irf = coefficients["radiation_damping/impulse_response_fun"]
            radiation_time = irf["t"][:].ravel()
            radiation_kernel = np.moveaxis(irf["K"][:], -1, 0) * rho
        with h5py.File(mass_properties_file) as properties:
            platform = properties["Platform/VolturnUS"]
            mass = float(platform["mass"][0, 0])
            inertia = platform["I_COG"][:]

        buoyancy = rho * gravity * volume
        restoring = np.zeros(6)
        restoring[2] = mass * gravity - buoyancy
        restoring[3:] = -np.cross(
            center_buoyancy - center_gravity, [0, 0, buoyancy],
        )
        drag = np.array([
            [9.23e5, 0, 0, 0, -8.92e6, 0],
            [0, 9.23e5, 0, 8.92e6, 0, 0],
            [0, 0, 2.30e6, 0, 0, 0],
            [0, 8.92e6, 0, 1.68e10, 0, 0],
            [-8.92e6, 0, 0, 0, 1.68e10, 0],
            [0, 0, 0, 0, 0, 4.80e10],
        ])
        return cls(
            mass=mass, inertia=inertia,
            equilibrium_pose=np.r_[center_gravity, np.zeros(3)],
            restoring_stiffness=stiffness, static_restoring=restoring,
            drag=drag, radiation_time=radiation_time,
            radiation_kernel=radiation_kernel, added_mass=added_mass,
        )

    def restoring_force(self, pose: np.ndarray) -> np.ndarray:
        """Return the logged resisting hydrostatic force on saved or live poses."""
        pose = np.asarray(pose, dtype=float)
        if pose.shape[-1:] != (6,) or not np.isfinite(pose).all():
            raise ValueError("platform pose must have six finite coordinates")
        return ((pose - self.equilibrium_pose) @ self.restoring_stiffness.T
                + self.static_restoring)

    def drag_force(self, velocity: np.ndarray) -> np.ndarray:
        """Return the logged quadratic viscous force on saved or live speeds."""
        velocity = np.asarray(velocity, dtype=float)
        if velocity.shape[-1:] != (6,) or not np.isfinite(velocity).all():
            raise ValueError("platform velocity must have six finite coordinates")
        return (velocity * np.abs(velocity)) @ self.drag.T

    def radiation_force(self, velocity: np.ndarray, dt: float,
                        memory_time: float = 60.0) -> np.ndarray:
        """Convolve platform velocity with the WEC-Sim radiation kernel."""
        velocity = np.asarray(velocity, dtype=float)
        if (velocity.ndim != 2 or velocity.shape[1] != 6
                or not np.isfinite(velocity).all()
                or not np.isfinite([dt, memory_time]).all()
                or dt <= 0 or memory_time <= 0
                or memory_time > self.radiation_time[-1]):
            raise ValueError("radiation needs finite six-DOF history and valid time steps")
        lags = min(len(velocity) - 1, int(round(memory_time / dt)))
        kernel = CubicSpline(self.radiation_time, self.radiation_kernel,
                             axis=0)(np.arange(lags + 1) * dt)
        force = np.zeros_like(velocity)
        for index in range(1, len(velocity)):
            last = min(index, lags)
            force[index] = dt * np.einsum(
                "tij,tj->i", kernel[:last + 1],
                velocity[index - last:index + 1][::-1],
            )
            force[index] -= dt / 2 * (
                kernel[0] @ velocity[index]
                + kernel[last] @ velocity[index - last]
            )
        return force

    def simulate(self, time: np.ndarray, wave_excitation: np.ndarray,
                 tower_base_load: np.ndarray, *,
                 mooring: MostStaticMooring | None = None,
                 radiation_memory: float = 60.0) -> MostPlatformResponse:
        """Advance published surge, heave, and pitch with a prescribed tower load.

        ``wave_excitation`` is a Python-generated six-component HDF5 force.
        ``tower_base_load`` is a six-component force/moment history in the
        platform frame. The tower history is an external input; this method
        does not model turbine-to-platform feedback or the three inactive
        platform coordinates.
        """
        time = np.asarray(time, dtype=float)
        wave_excitation = np.asarray(wave_excitation, dtype=float)
        tower_base_load = np.asarray(tower_base_load, dtype=float)
        if (time.ndim != 1 or len(time) < 2 or not np.isfinite(time).all()
                or not np.isclose(time[0], 0, rtol=0, atol=1e-12)
                or wave_excitation.shape != (len(time), 6)
                or tower_base_load.shape != (len(time), 6)
                or not np.isfinite(wave_excitation).all()
                or not np.isfinite(tower_base_load).all()):
            raise ValueError("MOST platform needs aligned finite time and six-component loads")
        dt = time[1] - time[0]
        if (dt <= 0 or not np.allclose(np.diff(time), dt, rtol=0, atol=1e-12)
                or not np.isfinite(radiation_memory)
                or radiation_memory <= 0):
            raise ValueError("MOST platform needs uniform positive time and radiation memory")
        memory = min(radiation_memory, time[-1])
        if memory > self.radiation_time[-1]:
            raise ValueError("radiation memory exceeds the hydrodynamic IRF")
        kernel = CubicSpline(self.radiation_time, self.radiation_kernel,
                             axis=0)(np.arange(round(memory / dt) + 1) * dt)
        mooring = MostStaticMooring() if mooring is None else mooring
        if not isinstance(mooring, MostStaticMooring):
            raise TypeError("mooring must be MostStaticMooring")

        jacobian = np.zeros((6, 3))
        jacobian[0, 0] = jacobian[2, 1] = jacobian[4, 2] = 1
        mooring_offset = -self.equilibrium_pose[:3]
        tower_offset = mooring_offset + np.array([0, 0, self.tower_base_height])

        def motion(coordinate, speed):
            return BodyMotion(jacobian @ coordinate, jacobian, np.zeros(6))

        def live_force(at_time, coordinate, speed):
            index = min(len(time) - 1, max(0, int(round(at_time / dt))))
            pose = self.equilibrium_pose + jacobian @ coordinate
            velocity = jacobian @ speed
            rotation = mooring._rotation(*pose[3:])
            arm_mooring = rotation @ mooring_offset
            arm_tower = rotation @ tower_offset
            mooring_pose = np.r_[pose[:3] + arm_mooring, pose[3:]]
            mooring_load = mooring.force(mooring_pose)[0]
            mooring_wrench = np.r_[
                mooring_load[:3],
                mooring_load[3:] + np.cross(arm_mooring, mooring_load[:3]),
            ]
            tower = tower_base_load[index]
            tower_force = rotation @ tower[:3]
            tower_wrench = np.r_[
                tower_force,
                rotation @ tower[3:] + np.cross(arm_tower, tower_force),
            ]
            return (wave_excitation[index] + mooring_wrench + tower_wrench
                    - self.restoring_force(pose) - self.drag_force(velocity))

        rigid_mass = np.diag([self.mass] * 3 + [0.] * 3)
        rigid_mass[3:, 3:] = self.inertia
        body = DynamicBody(
            rigid_mass=rigid_mass, added_mass=(self.added_mass,),
            damping=(np.zeros((6, 6)),), restoring=np.zeros((6, 6)),
            static_force=np.zeros(6), reference_position=self.equilibrium_pose,
            motion=motion, excitation=lambda _: np.zeros(6),
            radiation_kernel=kernel, state_excitation=live_force,
        )
        response = GeneralizedDynamics(
            (body,), 3, added_mass_delay=1e-7,
        ).integrate(dt=dt, end_time=float(time[-1]))
        return MostPlatformResponse(
            response.time, response.body_position[:, 0],
            response.body_velocity[:, 0], response.acceleration,
        )
