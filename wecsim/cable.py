"""Axial cable tension law used by the published WEC-Sim Cable application."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class WecSimCableTension:
    """Reproduce ``calcCableTens`` in the cable's local z coordinate.

    ``relative_position`` and ``relative_velocity`` are the cable block's
    displacement and speed relative to its initial endpoint separation.
    The returned force uses the source block's signed actuation convention.
    The source does not clamp the damping term after extension, so a rapidly
    contracting stretched cable can report positive force.
    """

    stiffness: float
    damping: float
    length: float
    initial_length: float

    def __post_init__(self):
        parameters = np.asarray((
            self.stiffness, self.damping, self.length, self.initial_length,
        ), dtype=float)
        if (not np.isfinite(parameters).all()
                or np.any(parameters < 0)):
            raise ValueError("cable coefficients and lengths must be finite and nonnegative")

    def force_z(self, relative_position, relative_velocity):
        position = np.asarray(relative_position, dtype=float)
        velocity = np.asarray(relative_velocity, dtype=float)
        if not np.isfinite(position).all() or not np.isfinite(velocity).all():
            raise ValueError("cable position and velocity must be finite")
        stretched_length = np.abs(position + self.initial_length)
        force = np.where(
            stretched_length > self.length,
            -self.stiffness * (stretched_length - self.length)
            - self.damping * velocity,
            0.,
        )
        return float(force) if np.ndim(force) == 0 else force
