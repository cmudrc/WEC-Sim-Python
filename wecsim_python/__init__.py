"""Compatibility import for the public :mod:`wecsim` package."""

from wecsim import (
    Body, BodyPoint, Coordinate, LinearPTO, Motion, MotionHistory, NoWave,
    PTOHistory, RegularWave, WEC, WECResult, WorldPoint,
)

__all__ = [
    "Body", "BodyPoint", "Coordinate", "LinearPTO", "Motion",
    "MotionHistory", "NoWave", "PTOHistory", "RegularWave", "WEC",
    "WECResult", "WorldPoint",
]
