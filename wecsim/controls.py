"""Discrete PTO control laws for published WEC-Sim application cases."""

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class DeclutchingState:
    """Controller memory before the next velocity sample."""

    elapsed_off: float = 0.0
    elapsed_on: float = 0.0
    previous_velocity: float = 0.0


@dataclass(frozen=True)
class DeclutchingControl:
    """Zero PTO force for a timed interval after a velocity sign change.

    This is the discrete law in the pinned Sphere Declutching application.
    It samples the PTO velocity once per simulation step. ``gain`` is the
    proportional damping applied while engaged.
    """

    gain: float
    declutch_time: float
    minimum_on_time: float = 0.2

    def __post_init__(self):
        if not all(math.isfinite(value) for value in (
                self.gain, self.declutch_time, self.minimum_on_time)):
            raise ValueError("declutching settings must be finite")
        if self.gain < 0 or self.declutch_time <= 0 or self.minimum_on_time < 0:
            raise ValueError("declutching settings are outside their supported range")

    def sample(self, velocity: float, state: DeclutchingState,
               dt: float) -> tuple[float, DeclutchingState]:
        """Return signed force and next memory for one major simulation step."""
        if not isinstance(state, DeclutchingState):
            raise TypeError("state must be DeclutchingState")
        if not (math.isfinite(velocity) and math.isfinite(dt) and dt > 0):
            raise ValueError("velocity must be finite and dt positive")
        reversal = (self._sign(velocity) != self._sign(state.previous_velocity)
                    and state.elapsed_off == 0
                    and state.elapsed_on > self.minimum_on_time)
        interval_active = 0 < state.elapsed_off < self.declutch_time
        if reversal or interval_active:
            return 0.0, DeclutchingState(
                state.elapsed_off + dt, 0.0, velocity,
            )
        return -self.gain * velocity, DeclutchingState(
            0.0, state.elapsed_on + dt, velocity,
        )

    @staticmethod
    def _sign(value: float) -> int:
        return (value > 0) - (value < 0)
