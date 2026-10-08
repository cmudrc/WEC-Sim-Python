"""Hydraulic components in the published OSWEC reverse-osmosis application."""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq


@dataclass(frozen=True)
class ReverseOsmosisMembrane:
    """Osmotic relief valve and hydraulic resistance in series.

    Pressures are gauge pressure relative to the permeate outlet. The valve
    follows the published hydraulic relief-block opening and laminar
    transition equations. ``permeate_flow`` solves the pressure split across
    the valve and linear resistance; it does not determine network pressure.
    """

    resistance: float
    osmotic_pressure: float
    regulation_range: float
    max_area: float
    leakage_area: float
    discharge_coefficient: float
    fluid_density: float
    laminar_pressure_ratio: float = .999
    atmospheric_pressure: float = 101_325.

    def __post_init__(self):
        values = np.asarray((
            self.resistance, self.osmotic_pressure, self.regulation_range,
            self.max_area, self.leakage_area, self.discharge_coefficient,
            self.fluid_density, self.laminar_pressure_ratio,
            self.atmospheric_pressure,
        ), dtype=float)
        if (not np.isfinite(values).all() or np.any(values[:7] <= 0)
                or self.leakage_area > self.max_area
                or not 0 < self.laminar_pressure_ratio < 1
                or self.atmospheric_pressure <= 0):
            raise ValueError("reverse-osmosis valve and resistance settings are invalid")

    def valve_flow(self, pressure_drop: float) -> float:
        """Signed flow through the pressure-controlled osmotic valve."""
        dp = float(pressure_drop)
        if not np.isfinite(dp):
            raise ValueError("valve pressure drop must be finite")
        opening = np.clip(
            (dp - self.osmotic_pressure) / self.regulation_range, 0, 1,
        )
        area = self.leakage_area + (self.max_area - self.leakage_area) * opening
        critical = ((self.atmospheric_pressure + dp / 2)
                    * (1 - self.laminar_pressure_ratio))
        return (self.discharge_coefficient * area
                * np.sqrt(2 / self.fluid_density) * dp
                / (dp * dp + critical * critical) ** .25)

    def permeate_flow(self, inlet_pressure):
        """Solve permeate flow from the measured or modeled inlet pressure."""
        pressure = np.asarray(inlet_pressure, dtype=float)
        if not np.isfinite(pressure).all() or np.any(pressure < 0):
            raise ValueError("membrane inlet pressure must be finite and nonnegative")
        flows = np.empty_like(pressure)
        for index, value in np.ndenumerate(pressure):
            if value == 0:
                flows[index] = 0
                continue
            drop = brentq(
                lambda dp: dp + self.resistance * self.valve_flow(dp) - value,
                0, float(value), xtol=1e-8,
            )
            flows[index] = self.valve_flow(drop)
        return float(flows) if flows.ndim == 0 else flows
