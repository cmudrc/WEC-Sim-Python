"""Compressible hydraulic cylinder equations used by WEC-Sim PTO-Sim.

This component maps prescribed chamber pressures and flows to cylinder force
and pressure rates. It is not a coupled hydraulic PTO network or a WEC runner.
"""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class CompressibleCylinder:
    area_a: float
    area_b: float
    bulk_modulus: float
    stroke: float
    offset: float

    def __post_init__(self):
        values = np.asarray((self.area_a, self.area_b, self.bulk_modulus,
                             self.stroke, self.offset), dtype=float)
        if not np.isfinite(values).all() or np.any(values <= 0):
            raise ValueError("cylinder areas, bulk modulus, stroke, and offset must be positive")

    def force(self, pressure_a, pressure_b):
        """Force in the pinned PTO-Sim output sign convention, in newtons."""
        return (self.area_b * np.asarray(pressure_b, dtype=float)
                - self.area_a * np.asarray(pressure_a, dtype=float))

    def pressure_rates(self, position, velocity, flow_a, flow_b):
        """Return A/B pressure rates for position, speed, and port flows.

        The volume floors match the pinned Simulink cylinder blocks. Inputs
        are SI metres, m/s, and m³/s; outputs are Pa/s.
        """
        position = np.asarray(position, dtype=float)
        velocity = np.asarray(velocity, dtype=float)
        flow_a = np.asarray(flow_a, dtype=float)
        flow_b = np.asarray(flow_b, dtype=float)
        length_a = np.maximum(self.offset - position, 1e-10)
        volume_b = np.maximum((position - self.offset + self.stroke)
                              * self.area_b, 1e-10)
        rate_a = (self.bulk_modulus * (self.area_a * velocity + flow_a)
                  / (self.area_a * length_a))
        rate_b = (self.bulk_modulus * (flow_b - self.area_b * velocity)
                  / volume_b)
        return rate_a, rate_b
