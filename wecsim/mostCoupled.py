"""Bidirectionally coupled MOST platform and turbine motion."""

from dataclasses import dataclass

import numpy as np

from .mostMooring import MostStaticMooring
from .mostPlatform import MostPlatformHydrodynamics, MostPlatformResponse
from .mostRotor import MostRotor, MostRotorResponse
from .mostTower import MostTowerReaction
from .turbSim import MostWindField


@dataclass(frozen=True)
class MostCoupledResponse:
    platform: MostPlatformResponse
    rotor: MostRotorResponse
    iterations: int
    position_residual: float  # m or rad, maximum change on last pass
    velocity_residual: float  # m/s or rad/s, maximum change on last pass


@dataclass(frozen=True)
class MostCoupled:
    """Advance the published MOST turbine with six-coordinate platform motion.

    Each pass advances the Python rotor from the current platform trajectory,
    then solves platform motion with the resulting blade loads and implicit
    turbine inertia. Passes continue until both trajectories are consistent.
    ``full_six_dof=False`` retains the earlier surge/heave/pitch reduction.
    """

    platform: MostPlatformHydrodynamics
    rotor: MostRotor
    tower: MostTowerReaction

    def simulate(self, time, wave_excitation, wind: MostWindField, *,
                 position_tolerance: float = 1e-6,
                 velocity_tolerance: float = 1e-6,
                 max_iterations: int = 12,
                 mooring: MostStaticMooring | None = None,
                 full_six_dof: bool = True) -> MostCoupledResponse:
        """Use only wave and wind inputs to advance the coupled case."""
        if (not isinstance(self.platform, MostPlatformHydrodynamics)
                or not isinstance(self.rotor, MostRotor)
                or not isinstance(self.tower, MostTowerReaction)):
            raise TypeError("MOST coupling needs platform, rotor, and tower models")
        if (not np.isfinite([position_tolerance, velocity_tolerance]).all()
                or position_tolerance <= 0 or velocity_tolerance <= 0
                or not isinstance(max_iterations, int)
                or max_iterations < 1):
            raise ValueError("MOST coupling needs positive tolerances and iterations")
        time = np.asarray(time, dtype=float)
        if time.ndim != 1 or time.size < 2:
            raise ValueError("MOST coupling needs a time history")
        position = np.tile(self.platform.equilibrium_pose, (time.size, 1))
        velocity = np.zeros_like(position)
        axes = tuple(range(6)) if full_six_dof else (0, 2, 4)
        for iteration in range(1, max_iterations + 1):
            rotor_response = self.rotor.simulate(
                time, position, velocity, wind,
            )
            platform_response = self.platform.simulate_with_turbine(
                time, wave_excitation, self.tower,
                rotor_response.rotor_speed, rotor_response.azimuth,
                rotor_response.generator_torque,
                rotor_response.blade_root_load, mooring=mooring,
                full_six_dof=full_six_dof,
            )
            position_residual = float(np.max(np.abs(
                platform_response.position[:, axes] - position[:, axes]
            )))
            velocity_residual = float(np.max(np.abs(
                platform_response.velocity[:, axes] - velocity[:, axes]
            )))
            if (position_residual <= position_tolerance
                    and velocity_residual <= velocity_tolerance):
                return MostCoupledResponse(
                    platform_response, rotor_response, iteration,
                    position_residual, velocity_residual,
                )
            position = platform_response.position
            velocity = platform_response.velocity
        raise RuntimeError(
            "MOST platform/rotor coupling did not converge in "
            f"{max_iterations} passes (position {position_residual:.3g}, "
            f"velocity {velocity_residual:.3g})"
        )
