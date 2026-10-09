"""Read the TurbSim full-field wind input used by the MOST application."""

from dataclasses import dataclass
from math import ceil, floor
from pathlib import Path
import struct

import numpy as np
from scipy.interpolate import RegularGridInterpolator


_HEADER = struct.Struct("<h4i12fi")


@dataclass(frozen=True)
class TurbSimWind:
    """Decoded wind; velocity axes are time, component (U/V/W), y, z."""

    velocity: np.ndarray
    tower_velocity: np.ndarray
    y: np.ndarray
    z: np.ndarray
    dt: float
    hub_height: float
    mean_wind_speed: float
    description: str


class MostWindField:
    """MOST's frozen-turbulence X planes, evaluated without a large 5-D copy.

    ``at_index(i)`` has axes component, X, Y, Z, matching one time step of
    MOST's ``Wind.SpatialDiscrUVW`` from ``RunTurbsim.m``.
    """

    def __init__(self, wind: TurbSimWind, speed: float = 8.0):
        if not np.isfinite(speed) or speed <= 0:
            raise ValueError("wind advection speed must be positive and finite")
        self.wind = wind
        self.speed = float(speed)
        self.x = np.arange(-30.0, 61.0, 10.0)
        self.discarded = ceil(wind.y[-1] / (self.speed * wind.dt))
        self.n_time = len(wind.velocity) - 2 * self.discarded + 1
        if self.n_time < 1:
            raise ValueError("TurbSim record is too short for MOST advection")
        first = floor(self.discarded + 1 - self.x[-1] / (self.speed * wind.dt) + .5)
        last = floor(self.discarded + self.n_time - self.x[0] / (self.speed * wind.dt) + .5)
        if first < 1 or last > len(wind.velocity):
            raise ValueError("MOST advection samples outside the TurbSim record")

    @property
    def time(self) -> np.ndarray:
        return np.arange(self.n_time) * self.wind.dt

    def at_index(self, index: int) -> np.ndarray:
        """Return all ten X planes at a zero-based output time index."""
        if not 0 <= index < self.n_time:
            raise IndexError("MOST wind time index is out of range")
        source = np.floor(
            self.discarded + 1 + index - self.x / (self.speed * self.wind.dt) + .5,
        ).astype(int) - 1
        return self.wind.velocity[source].transpose(1, 0, 2, 3)

    def sampler(self, index: int):
        """Return MOST's trilinear world-position wind sampler at one time.

        Outside the published X/Y/Z grid the sampler returns NaN, as the
        MATLAB BEM block's ``interp3`` call does.
        """
        values = self.at_index(index).transpose(1, 2, 3, 0)
        interpolator = RegularGridInterpolator(
            (self.x, self.wind.y, self.wind.z), values,
            bounds_error=False, fill_value=np.nan,
        )
        def sample(position):
            return interpolator(np.asarray(position, dtype=float).reshape(1, 3))[0]
        return sample

    def sampler_at(self, time: float):
        """Sample the spatial field at a simulation time in seconds.

        MOST's time-series wind input is linearly interpolated between the
        TurbSim advection frames before the BEM spatial interpolation.
        """
        if not np.isfinite(time) or not 0 <= time <= (self.n_time - 1) * self.wind.dt:
            raise ValueError("MOST wind time is outside the available record")
        fractional = time / self.wind.dt
        before = min(int(np.floor(fractional)), self.n_time - 1)
        after = min(before + 1, self.n_time - 1)
        weight = fractional - before
        if after == before or weight == 0:
            return self.sampler(before)
        values = ((1 - weight) * self.at_index(before)
                  + weight * self.at_index(after)).transpose(1, 2, 3, 0)
        interpolator = RegularGridInterpolator(
            (self.x, self.wind.y, self.wind.z), values,
            bounds_error=False, fill_value=np.nan,
        )
        def sample(position):
            return interpolator(np.asarray(position, dtype=float).reshape(1, 3))[0]
        return sample


def read_turbsim_bts(path: str | Path) -> TurbSimWind:
    """Read a little-endian TurbSim int16 ``.bts`` file.

    The returned grid ordering and scaling follow MOST's ``readfile_BTS.m``.
    Tower points are decoded even though the MOST reader does not return them.
    """
    data = Path(path).read_bytes()
    if len(data) < _HEADER.size:
        raise ValueError("TurbSim file has an incomplete header")
    (format_id, nz, ny, ntwr, nt, dz, dy, dt, mean_speed, hub_height, z1,
     slope_u, offset_u, slope_v, offset_v, slope_w, offset_w,
     description_length) = _HEADER.unpack_from(data)
    if format_id not in (7, 8):
        raise ValueError(f"unsupported TurbSim format identifier {format_id}")
    if min(nz, ny, nt) <= 0 or ntwr < 0 or not 0 <= description_length <= 4096:
        raise ValueError("invalid TurbSim grid dimensions or description length")
    slopes = np.array([slope_u, slope_v, slope_w], dtype=float)
    offsets = np.array([offset_u, offset_v, offset_w], dtype=float)
    if not np.isfinite(slopes).all() or np.any(slopes == 0):
        raise ValueError("invalid TurbSim velocity slopes")
    payload_start = _HEADER.size + description_length
    points_per_step = 3 * (ny * nz + ntwr)
    expected_size = payload_start + 2 * nt * points_per_step
    if len(data) != expected_size:
        raise ValueError(f"TurbSim file size is {len(data)}, expected {expected_size}")
    description = data[_HEADER.size:payload_start].decode("ascii", errors="replace")
    raw = np.frombuffer(data, dtype="<i2", count=nt * points_per_step,
                        offset=payload_start).reshape(nt, points_per_step)
    grid = raw[:, :3 * ny * nz].reshape(nt, nz, ny, 3).transpose(0, 3, 2, 1)
    velocity = (grid - offsets[None, :, None, None]) / slopes[None, :, None, None]
    tower = raw[:, 3 * ny * nz:].reshape(nt, ntwr, 3).transpose(0, 2, 1)
    tower_velocity = (tower - offsets[None, :, None]) / slopes[None, :, None]
    y = np.arange(ny) * dy - dy * (ny - 1) / 2
    z = np.arange(nz) * dz + z1
    return TurbSimWind(velocity, tower_velocity, y, z, dt, hub_height,
                       mean_speed, description)
