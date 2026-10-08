"""Equivalent-circuit electric generator used by the RM3 PTO-Sim case."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class EquivalentCircuitGenerator:
    armature_resistance: float
    armature_inductance: float
    torque_constant: float
    rotor_inertia: float
    shaft_damping: float

    def __post_init__(self):
        values = np.asarray((self.armature_resistance,
                             self.armature_inductance, self.torque_constant,
                             self.rotor_inertia, self.shaft_damping), dtype=float)
        if not np.isfinite(values).all() or (values[:4] <= 0).any() or values[4] < 0:
            raise ValueError("generator electrical and inertial parameters are invalid")

    def electromagnetic_torque(self, current):
        """Source-signed electromagnetic torque in N m."""
        return self.torque_constant * np.asarray(current, dtype=float)

    def current_rate(self, angular_speed, current, load_voltage):
        """Armature-current rate in A/s for speed in rad/s."""
        return -(self.torque_constant * np.asarray(angular_speed, dtype=float)
                 + self.armature_resistance * np.asarray(current, dtype=float)
                 + np.asarray(load_voltage, dtype=float)) / self.armature_inductance

    def speed_rate(self, angular_speed, current, drive_torque):
        """Shaft acceleration in rad/s² for positive hydraulic drive torque."""
        return (np.asarray(drive_torque, dtype=float)
                + self.electromagnetic_torque(current)
                - self.shaft_damping * np.asarray(angular_speed, dtype=float)) / self.rotor_inertia
