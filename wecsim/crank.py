"""Rotary-to-linear linkages in the published OSWEC PTO-Sim applications."""

from dataclasses import dataclass

import numpy as np


def _initial_horizontal(crank: float, offset: float, rod_length: float) -> float:
    parameters = np.asarray((crank, offset, rod_length), dtype=float)
    if (not np.isfinite(parameters).all() or crank <= 0 or rod_length <= 0
            or rod_length <= abs(crank - offset)):
        raise ValueError("crank, offset, and initial rod geometry are invalid")
    return float(np.sqrt(rod_length**2 - (crank - offset)**2))


@dataclass(frozen=True)
class FixedRodCrank:
    """Published slider-crank position measured from its zero-angle pose."""

    crank: float
    offset: float
    rod_length: float

    def __post_init__(self):
        _initial_horizontal(self.crank, self.offset, self.rod_length)

    def stroke_and_jacobian(self, angle):
        theta = np.asarray(angle, dtype=float)
        horizontal = self.crank * np.cos(theta) - self.offset
        under_root = self.rod_length**2 - horizontal**2
        if not np.isfinite(theta).all() or np.any(under_root <= 0):
            raise ValueError("crank angle is outside the slider geometry")
        rod_horizontal = np.sqrt(under_root)
        stroke = (self.crank * np.sin(theta) + rod_horizontal
                  - _initial_horizontal(self.crank, self.offset,
                                        self.rod_length))
        jacobian = (self.crank * np.cos(theta)
                    + self.crank * horizontal * np.sin(theta)
                    / rod_horizontal)
        return stroke, jacobian

    def torque(self, piston_force, angle):
        """Published fixed-crank PTO torque from signed cylinder force."""
        return -np.asarray(piston_force, dtype=float) * self.stroke_and_jacobian(angle)[1]


@dataclass(frozen=True)
class AdjustableRodCrank:
    """Published adjustable-rod length change from its zero-angle length."""

    crank: float
    offset: float
    initial_rod_length: float

    def __post_init__(self):
        _initial_horizontal(self.crank, self.offset,
                            self.initial_rod_length)

    def stroke_and_jacobian(self, angle):
        theta = np.asarray(angle, dtype=float)
        if not np.isfinite(theta).all():
            raise ValueError("crank angle must be finite")
        initial_horizontal = _initial_horizontal(
            self.crank, self.offset, self.initial_rod_length,
        )
        dx = self.crank * np.cos(theta) - self.offset
        dy = initial_horizontal - self.crank * np.sin(theta)
        length = np.hypot(dx, dy)
        if np.any(length <= 0):
            raise ValueError("adjustable rod collapsed to zero length")
        stroke = length - self.initial_rod_length
        jacobian = -self.crank * (dx * np.sin(theta)
                                  + dy * np.cos(theta)) / length
        return stroke, jacobian

    def torque(self, piston_force, angle):
        """Published adjustable-rod PTO torque from signed cylinder force."""
        return np.asarray(piston_force, dtype=float) * self.stroke_and_jacobian(angle)[1]
