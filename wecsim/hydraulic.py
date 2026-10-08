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


@dataclass(frozen=True)
class RectifyingCheckValve:
    """Four-port PTO-Sim valve with smooth pressure-dependent check openings."""

    discharge_coefficient: float
    area_max: float
    area_min: float
    pressure_max: float
    pressure_min: float
    density: float
    switch_gain: float
    opening_gain: float

    def __post_init__(self):
        values = np.asarray((self.discharge_coefficient, self.area_max,
                             self.area_min, self.pressure_max,
                             self.pressure_min, self.density,
                             self.switch_gain, self.opening_gain), dtype=float)
        if not np.isfinite(values).all() or (values[[0, 1, 5, 6, 7]] <= 0).any():
            raise ValueError("valve coefficient, maximum area, density, and gains must be positive")
        if (self.area_min < 0 or self.area_min > self.area_max
                or self.pressure_max <= self.pressure_min):
            raise ValueError("valve area and pressure bounds are invalid")

    def check_flow(self, pressure_drop):
        """Nonnegative flow through one source check valve, in m³/s."""
        drop = np.asarray(pressure_drop, dtype=float)
        midpoint = (self.pressure_max + self.pressure_min) / 2
        area = self.area_min + (self.area_max - self.area_min) / 2 * (
            1 + np.tanh(self.opening_gain * (drop - midpoint)))
        speed = np.sqrt(2 * drop * np.tanh(self.switch_gain * drop)
                        / self.density)
        return self.discharge_coefficient * area * speed

    def flows(self, pressure_a, pressure_b, pressure_high, pressure_low):
        """Return signed flows at cylinder A/B and high/low pressure ports."""
        a, b, high, low = np.broadcast_arrays(
            np.asarray(pressure_a, dtype=float),
            np.asarray(pressure_b, dtype=float),
            np.asarray(pressure_high, dtype=float),
            np.asarray(pressure_low, dtype=float),
        )
        a_to_high = self.check_flow(a - high)
        b_to_high = self.check_flow(b - high)
        low_to_b = self.check_flow(low - b)
        low_to_a = self.check_flow(low - a)
        return (low_to_a - a_to_high, low_to_b - b_to_high,
                a_to_high + b_to_high, -low_to_b - low_to_a)


@dataclass(frozen=True)
class GasChargedAccumulator:
    """Polytropic gas pressure versus integrated liquid inflow."""

    initial_gas_volume: float
    precharge_pressure: float
    exponent: float = 1.4

    def __post_init__(self):
        values = np.asarray((self.initial_gas_volume, self.precharge_pressure,
                             self.exponent), dtype=float)
        if not np.isfinite(values).all() or (values <= 0).any():
            raise ValueError("accumulator volume, precharge, and exponent must be positive")

    def pressure(self, liquid_inflow_volume):
        """Pressure in Pa; positive accumulated inlet volume compresses gas."""
        gas_volume = self.initial_gas_volume - np.asarray(
            liquid_inflow_volume, dtype=float)
        if np.any(gas_volume <= 0):
            raise ValueError("accumulator gas volume must remain positive")
        return (self.precharge_pressure
                * (self.initial_gas_volume / gas_volume) ** self.exponent)
