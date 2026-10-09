"""Air-chamber dynamics for the published floating oscillating-water-column case."""

from dataclasses import dataclass
from numbers import Real

import numpy as np


@dataclass(frozen=True)
class FloatingOwcChamber:
    """Compressible chamber driven by relative water-column motion.

    Displacement and speed are positive when the water column rises into the
    chamber. ``pressure`` is gauge pressure in Pa, and ``turbine_speed`` is in
    rad/s. The returned force acts on the water column in the opposite direction
    from positive displacement. This is the chamber law in the published
    ``OWC/FloatingOWC`` Simulink application; body and turbine dynamics are
    separate.
    """

    area: float
    initial_volume: float
    gamma: float
    ambient_pressure: float
    ambient_density: float
    turbine_diameter: float
    turbine_kappa: float

    def __post_init__(self) -> None:
        values = tuple(self.__dict__.values())
        if (any(isinstance(value, bool) or not isinstance(value, Real)
                for value in values)
                or not np.isfinite(values).all()
                or any(value <= 0 for value in values)):
            raise ValueError("floating OWC chamber parameters must be finite and positive")

    def volume(self, displacement: float | np.ndarray) -> np.ndarray:
        volume = self.initial_volume - self.area * np.asarray(displacement, dtype=float)
        if not np.isfinite(volume).all() or np.any(volume <= 0):
            raise ValueError("floating OWC chamber volume must remain positive")
        return volume

    def pressure_derivative(
        self,
        pressure: float | np.ndarray,
        displacement: float | np.ndarray,
        speed: float | np.ndarray,
        turbine_speed: float | np.ndarray,
    ) -> np.ndarray:
        """Return the gauge-pressure derivative in Pa/s."""
        gauge = np.asarray(pressure, dtype=float)
        speed = np.asarray(speed, dtype=float)
        omega = np.asarray(turbine_speed, dtype=float)
        absolute = self.ambient_pressure + gauge
        if (not np.isfinite(gauge).all() or not np.isfinite(speed).all()
                or not np.isfinite(omega).all() or np.any(absolute <= 0)
                or np.any(omega == 0)):
            raise ValueError("pressure, water-column speed, and nonzero turbine speed must be physical")
        density = self.ambient_density * (
            absolute / self.ambient_pressure) ** (1 / self.gamma)
        turbine_volume_flow = (
            gauge * self.turbine_diameter / (self.turbine_kappa * omega * density)
        )
        return self.gamma * absolute / self.volume(displacement) * (
            self.area * speed - turbine_volume_flow
        )

    def force_on_column(self, pressure: float | np.ndarray) -> np.ndarray:
        """Return the chamber's axial force on the water column in N."""
        gauge = np.asarray(pressure, dtype=float)
        if not np.isfinite(gauge).all():
            raise ValueError("pressure must be finite")
        return -self.area * gauge
