"""Morison force laws paired to pinned WEC-Sim source behavior.

The public WEC runner uses Cartesian option 1 for fixed elements in a
directional irregular sea. Each incident heading contributes its own
nonlinear drag before the forces are summed, as in ``irregWaveMorison.m``.
The public runner also couples axial elements to heave or to a body's
surge/heave/pitch motion in regular waves. The separate full six-DOF
``regular_morison_source_force`` retains source-specific conventions only
for diagnostics.
"""

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .irregularWave import IrregularComponents


@dataclass(frozen=True)
class MorisonElement:
    point: tuple[float, float, float]  # body-local, relative to center of gravity
    drag_coefficient: tuple[float, float, float]
    added_mass_coefficient: tuple[float, float, float]
    area: tuple[float, float, float]
    volume: float
    phase_mode: str = "directional"  # "matlab_shared" replays the pinned source loop


@dataclass(frozen=True)
class FixedMorisonResponse:
    time: np.ndarray
    wave_elevation: np.ndarray
    force: np.ndarray  # physical force and moment at the body's center of gravity


def finite_depth_wavenumber(
    omega: np.ndarray, *, water_depth: float, gravity: float = 9.81,
) -> np.ndarray:
    """Use WEC-Sim's 100 fixed-point dispersion iterations."""
    frequency = np.asarray(omega, dtype=float)
    if (frequency.ndim != 1 or len(frequency) == 0
            or not np.isfinite(frequency).all() or np.any(frequency <= 0)
            or not np.isfinite([water_depth, gravity]).all()
            or water_depth <= 0 or gravity <= 0):
        raise ValueError("wavenumber needs positive frequency, depth, and gravity")
    k = frequency**2 / gravity
    for _ in range(100):
        k = frequency**2 / gravity / np.tanh(k * water_depth)
    return k


def regular_morison_source_force(
    elements: Sequence[MorisonElement], *, time: float,
    position: Sequence[float], velocity: Sequence[float],
    acceleration: Sequence[float], wave_height: float, wave_period: float,
    direction: float, water_depth: float, ramp_time: float,
    rho: float = 1025.0, g: float = 9.81,
) -> np.ndarray:
    """Evaluate pinned ``regWaveMorison.m`` option 1 at a moving-body state.

    The six state components are surge, sway, heave, roll, pitch, and yaw.
    The current speed is zero. This source-law diagnostic does not solve a
    coupled WEC trajectory. It retains the source's rotation and local-point
    angular kinematics so that
    a comparison can expose, rather than conceal, source-specific behavior.
    """
    state = [np.asarray(value, dtype=float) for value in
             (position, velocity, acceleration)]
    if any(value.shape != (6,) or not np.isfinite(value).all()
           for value in state):
        raise ValueError("Morison position, velocity, and acceleration need finite six-vectors")
    if (not np.isfinite([time, wave_height, wave_period, direction,
                         water_depth, ramp_time, rho, g]).all()
            or time < 0 or wave_height < 0 or wave_period <= 0
            or water_depth <= 0 or ramp_time < 0 or rho <= 0 or g <= 0):
        raise ValueError("regular Morison wave and fluid settings are invalid")
    if not elements:
        raise ValueError("regular Morison force needs at least one element")

    pose, speed, accel = state
    roll, pitch, yaw = pose[3:]
    c4, s4 = np.cos(roll), np.sin(roll)
    c5, s5 = np.cos(pitch), np.sin(pitch)
    c6, s6 = np.cos(yaw), np.sin(yaw)
    # Keep the pinned MATLAB source's matrix, including its first-row term.
    rotation = np.array([
        [c5 * c6, c4 * s6 + s4 * s5 * c6, s4 * s6 - c4 * s5 * s6],
        [-c5 * s6, c4 * c6 - s4 * s5 * s6, s4 * c6 + c4 * s5 * s6],
        [s5, -s4 * c5, c4 * c5],
    ])
    omega = 2 * np.pi / wave_period
    k = finite_depth_wavenumber(
        np.array([omega]), water_depth=water_depth, gravity=g,
    )[0]
    heading = np.deg2rad(direction)
    wave_axis = np.array([np.cos(heading), np.sin(heading)])
    amplitude = wave_height / 2
    if ramp_time and time < ramp_time:
        amplitude *= (1 - np.cos(np.pi * time / ramp_time)) / 2
    result = np.zeros(6)
    for element in elements:
        if not isinstance(element, MorisonElement):
            raise TypeError("elements must be MorisonElement values")
        point, cd, ca, area = (
            np.asarray(value, dtype=float) for value in (
                element.point, element.drag_coefficient,
                element.added_mass_coefficient, element.area,
            )
        )
        if (any(value.shape != (3,) or not np.isfinite(value).all()
                for value in (point, cd, ca, area))
                or np.any(cd < 0) or np.any(ca < 0) or np.any(area < 0)
                or not np.isfinite(element.volume) or element.volume <= 0):
            raise ValueError("Morison element coefficients and geometry are invalid")
        rotated_point = rotation @ point
        world = pose[:3] + rotated_point
        if world[2] > 0:
            continue
        angular_cross_point = np.cross(speed[3:], point)
        body_velocity = speed[:3] + angular_cross_point
        body_acceleration = (accel[:3] + np.cross(accel[3:], point)
                             + np.cross(speed[3:], angular_cross_point))
        kh, kz = k * water_depth, k * world[2]
        if kh > np.pi:
            horizontal = vertical = np.exp(kz)
        else:
            horizontal = np.cosh(kz + kh) / np.cosh(kh)
            vertical = np.sinh(kz + kh) / np.cosh(kh)
        phase = omega * time - k * (world[:2] @ wave_axis)
        horizontal_velocity = amplitude * horizontal * np.cos(phase) * g * k / omega
        vertical_velocity = -amplitude * vertical * np.sin(phase) * g * k / omega
        horizontal_acceleration = -amplitude * horizontal * np.sin(phase) * g * k
        vertical_acceleration = -amplitude * vertical * np.cos(phase) * g * k
        fluid_velocity = np.r_[horizontal_velocity * wave_axis, vertical_velocity]
        fluid_acceleration = np.r_[horizontal_acceleration * wave_axis,
                                   vertical_acceleration]
        relative_velocity = fluid_velocity - body_velocity
        # These row-vector coefficient transforms are the source block's
        # option-1 convention; they are not a rotated drag tensor.
        area_rot = np.abs(area @ rotation)
        cd_rot = np.abs(cd) @ rotation
        ca_rot = np.abs(ca @ rotation)
        force = (0.5 * rho * cd_rot * area_rot * relative_velocity
                 * np.abs(relative_velocity)
                 + rho * element.volume * (
                     fluid_acceleration
                     + ca_rot * (fluid_acceleration - body_acceleration)))
        result[:3] += force
        result[3:] += np.cross(rotated_point, force)
    return result


def _axial_heave_element(element: MorisonElement):
    if not isinstance(element, MorisonElement):
        raise TypeError("elements must be MorisonElement values")
    point = np.asarray(element.point, dtype=float)
    cd = np.asarray(element.drag_coefficient, dtype=float)
    ca = np.asarray(element.added_mass_coefficient, dtype=float)
    area = np.asarray(element.area, dtype=float)
    if (any(value.shape != (3,) or not np.isfinite(value).all()
            for value in (point, cd, ca, area))
            or not np.allclose(point[:2], 0, rtol=0, atol=1e-12)
            or any(not np.allclose(value[:2], 0, rtol=0, atol=1e-12)
                   for value in (cd, ca, area))
            or np.any(cd < 0) or np.any(ca < 0) or np.any(area < 0)
            or not np.isfinite(element.volume) or element.volume <= 0):
        raise ValueError("heave Morison elements need axial, nonnegative geometry")
    return point[2], cd[2], ca[2], area[2], element.volume


def no_wave_heave_morison_terms(
    elements: Sequence[MorisonElement], *, center_z: float,
    heave: float, speed: float, rho: float,
) -> tuple[float, float]:
    """Return drag force and added mass for submerged axial heave elements.

    In still water, Cartesian option 1 reduces to ``-Cd*A*rho*v*|v|/2``
    plus ``-Ca*V*rho*a``. The latter belongs on the equation's mass side.
    Submergence follows WEC-Sim's mean-waterline switch at the element point.
    """
    if (not np.isfinite([center_z, heave, speed, rho]).all() or rho <= 0):
        raise ValueError("Morison heave state and density must be finite")
    drag_force = 0.0
    added_mass = 0.0
    for element in elements:
        point_z, cd, ca, area, volume = _axial_heave_element(element)
        if center_z + heave + point_z > 0:
            continue
        drag_force -= 0.5 * rho * cd * area * speed * abs(speed)
        added_mass += rho * volume * ca
    return drag_force, added_mass


def regular_wave_heave_morison_terms(
    elements: Sequence[MorisonElement], *, center_z: float,
    heave: float, speed: float, time: float, wave_height: float,
    wave_period: float, ramp_time: float, water_depth: float,
    rho: float, g: float = 9.81,
) -> tuple[float, float]:
    """Return axial regular-wave Morison force without body acceleration and mass.

    The element follows the moving heave coordinate. The fluid vertical
    velocity and acceleration are evaluated at its current submerged point.
    The returned added mass multiplies body acceleration in the dynamics.
    """
    if (not np.isfinite([center_z, heave, speed, time, wave_height,
                         wave_period, ramp_time, rho, g]).all()
            or time < 0 or wave_height < 0 or wave_period <= 0
            or ramp_time < 0 or rho <= 0 or g <= 0
            or not (np.isfinite(water_depth) or np.isposinf(water_depth))
            or water_depth <= 0):
        raise ValueError("regular-wave Morison heave settings are invalid")
    omega = 2 * np.pi / wave_period
    k = (omega * omega / g if np.isposinf(water_depth) else
         finite_depth_wavenumber(
             np.array([omega]), water_depth=water_depth, gravity=g,
         )[0])
    amplitude = wave_height / 2
    if ramp_time and time < ramp_time:
        amplitude *= (1 - np.cos(np.pi * time / ramp_time)) / 2
    force = 0.0
    added_mass = 0.0
    for element in elements:
        point_z, cd, ca, area, volume = _axial_heave_element(element)
        world_z = center_z + heave + point_z
        if world_z > 0:
            continue
        if k * water_depth > np.pi:
            vertical = np.exp(k * world_z)
        else:
            kh = k * water_depth
            vertical = np.sinh(k * world_z + kh) / np.cosh(kh)
        phase = omega * time
        fluid_speed = (-amplitude * vertical * np.sin(phase)
                       * g * k / omega)
        fluid_acceleration = (-amplitude * vertical * np.cos(phase)
                              * g * k)
        relative_speed = fluid_speed - speed
        force += (0.5 * rho * cd * area * relative_speed * abs(relative_speed)
                  + rho * volume * (1 + ca) * fluid_acceleration)
        added_mass += rho * volume * ca
    return force, added_mass


def regular_wave_axial_morison_terms(
    elements: Sequence[MorisonElement], *, position: Sequence[float],
    velocity: Sequence[float], time: float, wave_height: float,
    wave_period: float, ramp_time: float, water_depth: float,
    rho: float, g: float = 9.81,
) -> tuple[np.ndarray, np.ndarray]:
    """Physical axial-element wrench and added mass for surge/heave/pitch.

    The point and its axial direction rotate with body pitch. The returned
    wrench excludes the body's translational and angular acceleration; the
    six-by-six matrix places that acceleration dependence on the mass side.
    This uses a proper pitch rotation, independently of the pinned MATLAB
    source function's nonorthogonal general rotation diagnostic.
    """
    pose = np.asarray(position, dtype=float)
    speed = np.asarray(velocity, dtype=float)
    if (pose.shape != (6,) or speed.shape != (6,)
            or not np.isfinite(pose).all() or not np.isfinite(speed).all()
            or not np.allclose(pose[[1, 3, 5]], 0, rtol=0, atol=1e-12)
            or not np.allclose(speed[[1, 3, 5]], 0, rtol=0, atol=1e-12)):
        raise ValueError("axial moving Morison needs finite surge/heave/pitch state")
    if (not np.isfinite([time, wave_height, wave_period, ramp_time, rho, g]).all()
            or time < 0 or wave_height < 0 or wave_period <= 0
            or ramp_time < 0 or rho <= 0 or g <= 0
            or not (np.isfinite(water_depth) or np.isposinf(water_depth))
            or water_depth <= 0):
        raise ValueError("regular-wave Morison settings are invalid")
    omega = 2 * np.pi / wave_period
    k = (omega * omega / g if np.isposinf(water_depth) else
         finite_depth_wavenumber(
             np.array([omega]), water_depth=water_depth, gravity=g,
         )[0])
    amplitude = wave_height / 2
    if ramp_time and time < ramp_time:
        amplitude *= (1 - np.cos(np.pi * time / ramp_time)) / 2
    pitch = pose[4]
    axis = np.array([np.sin(pitch), 0.0, np.cos(pitch)])
    angular_velocity = np.array([0.0, speed[4], 0.0])
    wrench = np.zeros(6)
    added_mass = np.zeros((6, 6))
    for element in elements:
        point_z, cd, ca, area, volume = _axial_heave_element(element)
        arm = point_z * axis
        world = pose[:3] + arm
        if world[2] > 0:
            continue
        kh = k * water_depth
        if kh > np.pi:
            horizontal = vertical = np.exp(k * world[2])
        else:
            kz = k * world[2]
            horizontal = np.cosh(kz + kh) / np.cosh(kh)
            vertical = np.sinh(kz + kh) / np.cosh(kh)
        phase = omega * time - k * world[0]
        wave_speed = amplitude * g * k / omega
        wave_acceleration = amplitude * g * k
        fluid_velocity = np.array([
            wave_speed * horizontal * np.cos(phase), 0.0,
            -wave_speed * vertical * np.sin(phase),
        ])
        fluid_acceleration = np.array([
            -wave_acceleration * horizontal * np.sin(phase), 0.0,
            -wave_acceleration * vertical * np.cos(phase),
        ])
        point_velocity = speed[:3] + np.cross(angular_velocity, arm)
        bias_acceleration = np.cross(
            angular_velocity, np.cross(angular_velocity, arm),
        )
        axial_speed = axis @ (fluid_velocity - point_velocity)
        force = (0.5 * rho * cd * area * axial_speed * abs(axial_speed) * axis
                 + rho * volume * fluid_acceleration
                 + rho * volume * ca * axis
                 * (axis @ (fluid_acceleration - bias_acceleration)))
        wrench[:3] += force
        wrench[3:] += np.cross(arm, force)
        added_mass[:3, :3] += rho * volume * ca * np.outer(axis, axis)
    return wrench, added_mass


def solve_fixed_morison_irregular(
    components: IrregularComponents,
    elements: Sequence[MorisonElement],
    *,
    center_gravity: Sequence[float],
    water_depth: float,
    dt: float,
    end_time: float,
    ramp_time: float,
    rho: float = 1025.0,
    g: float = 9.81,
) -> FixedMorisonResponse:
    """Evaluate a stationary body's Cartesian Morison force at every sample."""
    center = np.asarray(center_gravity, dtype=float)
    if center.shape != (3,) or not np.isfinite(center).all():
        raise ValueError("center_gravity needs three finite coordinates")
    if (not np.isfinite([dt, end_time, ramp_time, rho, g]).all()
            or dt <= 0 or end_time < 0 or ramp_time < 0 or rho <= 0 or g <= 0):
        raise ValueError("fixed Morison time and fluid settings are invalid")
    steps = round(end_time / dt)
    if not np.isclose(steps * dt, end_time, rtol=0, atol=1e-10):
        raise ValueError("end_time must be a multiple of dt")
    if not elements:
        raise ValueError("a fixed Morison body needs at least one element")
    prepared = []
    for element in elements:
        if not isinstance(element, MorisonElement):
            raise TypeError("elements must be MorisonElement values")
        point = np.asarray(element.point, dtype=float)
        drag = np.asarray(element.drag_coefficient, dtype=float)
        added = np.asarray(element.added_mass_coefficient, dtype=float)
        area = np.asarray(element.area, dtype=float)
        if (any(value.shape != (3,) or not np.isfinite(value).all()
                for value in (point, drag, added, area))
                or np.any(drag < 0) or np.any(added < 0) or np.any(area < 0)
                or not np.isfinite(element.volume) or element.volume <= 0
                or element.phase_mode not in ("directional", "matlab_shared")):
            raise ValueError("Morison point, coefficients, area, and volume are invalid")
        prepared.append((point, drag, added, area, float(element.volume),
                         element.phase_mode))

    omega = np.asarray(components.omega, dtype=float)
    spectrum = np.asarray(components.spectral_amplitude, dtype=float)
    width = np.asarray(components.d_omega, dtype=float)
    directions = np.asarray(components.directions, dtype=float)
    spreading = np.asarray(components.spreading, dtype=float)
    phase = np.asarray(components.phase, dtype=float)
    if (omega.ndim != 1 or len(omega) < 2
            or spectrum.shape != omega.shape or width.shape != omega.shape
            or directions.ndim != 1 or spreading.shape != directions.shape
            or phase.shape != (len(omega), len(directions))
            or any(not np.isfinite(value).all() for value in
                   (omega, spectrum, width, directions, spreading, phase))
            or np.any(omega <= 0) or np.any(spectrum < 0)
            or np.any(width <= 0) or np.any(spreading < 0)
            or not np.isclose(spreading.sum(), 1, atol=1e-12)):
        raise ValueError("invalid directional irregular-wave components")
    k = finite_depth_wavenumber(omega, water_depth=water_depth, gravity=g)
    time = np.arange(steps + 1) * dt
    ramp = np.ones_like(time)
    if ramp_time:
        before = time < ramp_time
        ramp[before] = (1 - np.cos(np.pi * time[before] / ramp_time)) / 2
    elevation = np.zeros_like(time)
    force = np.zeros((len(time), 6))
    headings = np.deg2rad(directions)

    for start in range(0, len(time), 512):
        stop = min(start + 512, len(time))
        t = time[start:stop, None]
        local_ramp = ramp[start:stop]
        for heading_index, heading in enumerate(headings):
            direction = np.array([np.cos(heading), np.sin(heading)])
            amplitude = np.sqrt(spreading[heading_index] * spectrum * width)
            origin_phase = t * omega + phase[:, heading_index]
            elevation[start:stop] += local_ramp * (
                np.cos(origin_phase) @ amplitude
            )
            for point, drag, added, area, volume, phase_mode in prepared:
                world = center + point
                if world[2] > 0:
                    continue
                kh = k * water_depth
                kz = k * world[2]
                horizontal = np.where(
                    kh > np.pi, np.exp(kz),
                    np.cosh(kz + kh) / np.cosh(kh),
                )
                vertical = np.where(
                    kh > np.pi, np.exp(kz),
                    np.sinh(kz + kh) / np.cosh(kh),
                )
                # The pinned irregWaveMorison.m indexes randPhase(jj,1)
                # inside its heading loop. Keep that source behavior opt-in.
                force_phase = (phase[:, 0] if phase_mode == "matlab_shared"
                               else phase[:, heading_index])
                argument = (t * omega + force_phase
                            - k * (world[:2] @ direction))
                cosine = np.cos(argument)
                sine = np.sin(argument)
                horizontal_velocity = local_ramp * (
                    cosine @ (amplitude * horizontal * g * k / omega)
                )
                vertical_velocity = -local_ramp * (
                    sine @ (amplitude * vertical * g * k / omega)
                )
                horizontal_acceleration = -local_ramp * (
                    sine @ (amplitude * horizontal * g * k)
                )
                vertical_acceleration = -local_ramp * (
                    cosine @ (amplitude * vertical * g * k)
                )
                velocity = np.column_stack((
                    horizontal_velocity * direction[0],
                    horizontal_velocity * direction[1],
                    vertical_velocity,
                ))
                acceleration = np.column_stack((
                    horizontal_acceleration * direction[0],
                    horizontal_acceleration * direction[1],
                    vertical_acceleration,
                ))
                contribution = (rho * volume * (1 + added) * acceleration
                                + 0.5 * rho * drag * area * velocity
                                * np.abs(velocity))
                force[start:stop, :3] += contribution
                force[start:stop, 3:] += np.cross(point, contribution)
    return FixedMorisonResponse(time, elevation, force)
