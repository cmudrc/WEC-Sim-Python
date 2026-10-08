"""Instantaneous free-surface forces on a triangulated heaving body.

This implements WEC-Sim's nonlinearHydro=2 regular-wave force construction
for a body constrained to heave. Mesh triangles are in body coordinates about
the center of gravity, as required by WEC-Sim's geometry import.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
import trimesh


@dataclass(frozen=True)
class HeaveMeshHydro:
    centers: np.ndarray
    area_vectors: np.ndarray
    center_z: float
    rho: float
    gravity: float
    depth: float
    omega: float
    wave_number: float
    amplitude: float
    ramp_time: float
    mass: float
    drag_coefficient: float
    drag_area: float

    @classmethod
    def from_stl(cls, path: str | Path, *, center_z: float, rho: float,
                 gravity: float, depth: float, period: float, height: float,
                 ramp_time: float, mass: float | None,
                 drag_coefficient: float = 0, drag_area: float = 0):
        mesh = trimesh.load_mesh(path, process=False)
        if (not isinstance(mesh, trimesh.Trimesh) or len(mesh.faces) == 0
                or not np.isfinite(mesh.vertices).all()
                or not np.isfinite(mesh.area_faces).all()):
            raise ValueError("geometry_file must contain a finite triangle mesh")
        centers = mesh.triangles_center
        areas = mesh.face_normals * mesh.area_faces[:, None]
        displaced_mass = rho * np.sum(
            np.minimum(centers[:, 2] + center_z, 0) * areas[:, 2]
        )
        if not np.isfinite(displaced_mass) or displaced_mass <= 0:
            raise ValueError("geometry_file must produce positive equilibrium mass")
        omega = 2 * np.pi / period
        wave_number = brentq(
            lambda k: gravity * k * np.tanh(k * depth) - omega**2,
            np.finfo(float).eps, max(1.0, 2 * omega**2 / gravity),
        )
        return cls(centers, areas, center_z, rho, gravity, depth,
                   omega, wave_number, height / 2, ramp_time,
                   displaced_mass if mass is None else mass,
                   drag_coefficient, drag_area)

    def ramp(self, at_time: float) -> float:
        return (1.0 if self.ramp_time == 0 or at_time >= self.ramp_time
                else (1 - np.cos(np.pi * at_time / self.ramp_time)) / 2)

    def forces(self, at_time: float, heave: float, speed: float):
        """Return buoyancy-minus-weight, FK correction, and applied drag."""
        ramp = self.ramp(at_time)
        eta = (self.amplitude * ramp
               * np.cos(self.wave_number * self.centers[:, 0]
                        - self.omega * at_time))
        mean_z = self.centers[:, 2] + self.center_z
        moved_z = mean_z + heave

        # WEC-Sim clips pressure above the local incident free surface.
        submerged_z = np.where(moved_z > eta, 0, moved_z)
        buoyancy = self.rho * self.gravity * np.sum(
            submerged_z * self.area_vectors[:, 2]
        ) - self.mass * self.gravity

        mean_pressure = self.rho * self.gravity * eta * (
            np.cosh(self.wave_number * (mean_z + self.depth))
            / np.cosh(self.wave_number * self.depth)
        )
        mean_pressure = np.where(mean_z > 0, 0, mean_pressure)
        stretched_z = (moved_z - eta) * self.depth / (self.depth + eta)
        moved_pressure = self.rho * self.gravity * eta * (
            np.cosh(self.wave_number * (stretched_z + self.depth))
            / np.cosh(self.wave_number * self.depth)
        )
        moved_pressure = np.where(stretched_z > 0, 0, moved_pressure)
        fk_correction = -np.sum(
            (moved_pressure - mean_pressure) * self.area_vectors[:, 2]
        )
        # The source body block ramps its nonlinear FK output after waveClass
        # has already ramped the incident elevation.
        fk_correction *= ramp
        drag = (-0.5 * self.rho * self.drag_coefficient
                * self.drag_area * speed * abs(speed))
        return buoyancy, fk_correction, drag
