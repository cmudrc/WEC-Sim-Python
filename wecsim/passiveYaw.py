"""Directional excitation for a body that yaws with respect to the waves."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class PassiveYawExcitation:
    """Regular-wave force from BEM headings at the body's current yaw angle.

    The hydrodynamic force is interpolated at incident heading minus body yaw
    and rotated from body axes into world axes. Heading interpolation is
    continuous; the pinned MATLAB block instead reuses coefficients until yaw
    changes by its configured threshold.
    """

    directions: np.ndarray
    real: np.ndarray  # six DOFs by heading, in N/m or Nm/m
    imaginary: np.ndarray
    incident_direction: float
    omega: float
    amplitude: float
    ramp_time: float

    @classmethod
    def from_hydro_data(cls, hydro_data, *, omega, incident_direction,
                        amplitude, ramp_time, rho, g):
        parameters = hydro_data["simulation_parameters"]
        headings = np.asarray(parameters["wave_dir"], dtype=float).ravel()
        frequency = np.asarray(parameters["w"], dtype=float).ravel()
        if (headings.size < 3 or frequency.size < 2
                or not np.isfinite(headings).all()
                or not np.isfinite(frequency).all()
                or not np.all(np.diff(headings) > 0)
                or not np.all(np.diff(frequency) > 0)
                or np.max(np.diff(np.r_[headings, headings[0] + 360])) > 15
                or omega < frequency[0] or omega > frequency[-1]):
            raise ValueError(
                "passive yaw needs full-circle BEM headings and an in-range frequency"
            )
        excitation = hydro_data["hydro_coeffs"]["excitation"]
        real = np.asarray(excitation["re"], dtype=float)
        imaginary = np.asarray(excitation["im"], dtype=float)
        if (real.shape != (6, headings.size, frequency.size)
                or imaginary.shape != real.shape
                or not np.isfinite(real).all() or not np.isfinite(imaginary).all()):
            raise ValueError("passive-yaw excitation coefficients are incomplete")
        real = np.array([
            [np.interp(omega, frequency, row) for row in dof]
            for dof in real
        ]) * rho * g
        imaginary = np.array([
            [np.interp(omega, frequency, row) for row in dof]
            for dof in imaginary
        ]) * rho * g
        return cls(headings, real, imaginary, incident_direction,
                   omega, amplitude, ramp_time)

    def force(self, time: float, yaw: float) -> np.ndarray:
        """Return the six-component excitation in the world frame."""
        relative_heading = (self.incident_direction - np.degrees(yaw)) % 360
        real = np.array([
            np.interp(relative_heading, self.directions, row, period=360)
            for row in self.real
        ])
        imaginary = np.array([
            np.interp(relative_heading, self.directions, row, period=360)
            for row in self.imaginary
        ])
        ramp = (1.0 if self.ramp_time == 0 or time >= self.ramp_time
                else (1 - np.cos(np.pi * time / self.ramp_time)) / 2)
        local = self.amplitude * ramp * (
            real * np.cos(self.omega * time)
            - imaginary * np.sin(self.omega * time)
        )
        world = local.copy()
        c, s = np.cos(yaw), np.sin(yaw)
        for first in (0, 3):
            world[first] = c * local[first] - s * local[first + 1]
            world[first + 1] = s * local[first] + c * local[first + 1]
        return world
