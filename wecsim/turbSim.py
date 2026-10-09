"""Read the TurbSim full-field wind input used by the MOST application."""

from dataclasses import dataclass
from pathlib import Path
import struct

import numpy as np


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
