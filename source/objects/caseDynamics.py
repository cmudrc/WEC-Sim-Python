"""Case-driven entry point for the currently supported device layouts.

The case describes wave, body, constraint, and PTO properties rather than
naming RM3, OSWEC, or Sphere. Each supported layout feeds the common
``GeneralizedDynamics`` engine through a validated adapter. Unsupported
physics fails explicitly instead of silently omitting a force or DOF.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import numpy as np

from .hingePitch import solve_hinged_pitch_from_excitation
from .irregularWave import pm_equal_energy_components, synthesize_irregular_response
from .linearHeave import solve_heave_free_decay
from .rm3Regular import solve_rm3_regular


AXES = ("surge", "sway", "heave", "roll", "pitch", "yaw")


@dataclass(frozen=True)
class CaseResponse:
    time: np.ndarray
    body_position: np.ndarray  # (time, body, six WEC-Sim coordinates)
    body_velocity: np.ndarray
    hydro_files: tuple[Path, ...]
    pto_force: np.ndarray | None = None
    pto_label: str | None = None
    wave_elevation: np.ndarray | None = None
    total_heave_force: np.ndarray | None = None
    auxiliary_files: tuple[Path, ...] = ()


def _section(value, name, required, allowed):
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be an object")
    missing = required - value.keys()
    extra = value.keys() - allowed
    if missing or extra:
        raise ValueError(f"{name}: missing {sorted(missing)}, unsupported {sorted(extra)}")
    return value


def _number(value, name, *, positive=False, nonnegative=False):
    if isinstance(value, bool):
        raise ValueError(f"{name} must be numeric")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if (not np.isfinite(number) or (positive and number <= 0)
            or (nonnegative and number < 0)):
        raise ValueError(f"{name} is outside its supported range")
    return number


def _hydro_file(body, base_dir):
    _section(body, "body", {"hydro_file"},
             {"hydro_file", "hydro_body", "mass", "pitch_inertia"})
    raw = body["hydro_file"]
    if not isinstance(raw, str) or not raw:
        raise ValueError("body.hydro_file must be a file path")
    path = (base_dir / raw).expanduser().resolve(strict=True)
    if not path.is_file():
        raise ValueError(f"body.hydro_file is not a file: {path}")
    return path


def _body_number(body, expected):
    number = body.get("hydro_body", expected)
    if isinstance(number, bool) or number != expected:
        raise ValueError(f"body.hydro_body must be {expected} for this layout")


def _location(constraint):
    location = np.asarray(constraint.get("location", [0, 0, 0]), dtype=float)
    if location.shape != (3,) or not np.isfinite(location).all():
        raise ValueError("constraint.location must have three finite coordinates")
    if not np.isclose(location[:2], 0, atol=1e-10).all():
        raise ValueError("supported joints must lie on the body x=y=0 axis")
    return location


def _pto(pto, kind):
    _section(pto, "pto", {"kind", "damping"},
             {"kind", "damping", "stiffness"})
    if pto["kind"] != kind:
        raise ValueError(f"this layout requires a {kind} PTO")
    return (
        _number(pto["damping"], "pto.damping", nonnegative=True),
        _number(pto.get("stiffness", 0), "pto.stiffness", nonnegative=True),
    )


def run_case(case: Mapping, *, base_dir: str | Path = ".") -> CaseResponse:
    """Run one explicitly supported wave/device configuration.

    ``base_dir`` resolves relative hydro and phase file paths. The accepted
    layouts are one-body heave free decay, one-body fixed-hinge pitch, and a
    two-body floating joint with independent heaves and shared surge/pitch.
    Every one of these layouts uses the generalized dynamics engine.
    """
    case = _section(
        case, "case", {"simulation", "wave", "bodies", "constraint"},
        {"name", "simulation", "wave", "bodies", "constraint", "pto", "body_to_body"},
    )
    sim = _section(case["simulation"], "simulation", {"dt", "end_time"},
                   {"dt", "end_time", "ramp_time", "rho", "g", "radiation_memory"})
    dt = _number(sim["dt"], "simulation.dt", positive=True)
    end_time = _number(sim["end_time"], "simulation.end_time", nonnegative=True)
    ramp_time = _number(sim.get("ramp_time", 100), "simulation.ramp_time", nonnegative=True)
    rho = _number(sim.get("rho", 1000), "simulation.rho", positive=True)
    g = _number(sim.get("g", 9.81), "simulation.g", positive=True)
    wave = _section(case["wave"], "wave", {"type"},
                    {"type", "height", "period", "direction", "directions",
                     "spreading", "seed", "phase_file", "frequency_count"})
    constraint = _section(case["constraint"], "constraint", {"kind"},
                          {"kind", "location", "initial_displacement"})
    bodies = case["bodies"]
    if not isinstance(bodies, list) or not bodies:
        raise ValueError("bodies must be a nonempty list")
    base = Path(base_dir).expanduser().resolve()
    hydro = tuple(_hydro_file(body, base) for body in bodies)
    b2b = case.get("body_to_body", False)
    if not isinstance(b2b, bool):
        raise ValueError("body_to_body must be a boolean")
    kind = constraint["kind"]

    if kind == "heave":
        if len(bodies) != 1 or wave["type"] != "none" or b2b or "pto" in case:
            raise ValueError("heave free decay needs one body, no waves, and no PTO")
        if set(wave) != {"type"} or set(constraint) - {"kind", "initial_displacement"}:
            raise ValueError("heave free decay has no wave or joint-location settings")
        if set(bodies[0]) - {"hydro_file", "hydro_body", "mass"}:
            raise ValueError("heave free decay does not use body inertia settings")
        _body_number(bodies[0], 1)
        if bodies[0].get("mass", "equilibrium") != "equilibrium":
            raise ValueError("heave free decay currently requires equilibrium mass")
        if "ramp_time" in sim:
            raise ValueError("ramp_time is inapplicable to no-wave free decay")
        displacement = _number(
            constraint.get("initial_displacement", 0),
            "constraint.initial_displacement",
        )
        solved = solve_heave_free_decay(
            hydro[0], displacement, dt=dt, end_time=end_time,
            cic_end_time=_number(
                sim.get("radiation_memory", 15), "simulation.radiation_memory",
                positive=True,
            ),
            rho=rho, g=g,
        )
        position = np.zeros((len(solved.time), 1, 6))
        velocity = np.zeros_like(position)
        position[:, 0, 2] = solved.position
        velocity[:, 0, 2] = solved.velocity
        return CaseResponse(
            solved.time, position, velocity, hydro,
            total_heave_force=solved.force_total,
        )

    if kind == "fixed_hinge":
        if len(bodies) != 1 or wave["type"] != "pm" or b2b:
            raise ValueError("fixed-hinge pitch needs one body and PM waves")
        _body_number(bodies[0], 1)
        mass = _number(bodies[0].get("mass"), "body.mass", positive=True)
        inertia = _number(bodies[0].get("pitch_inertia"),
                          "body.pitch_inertia", positive=True)
        if "initial_displacement" in constraint:
            raise ValueError("fixed-hinge initial displacement is not yet supported")
        location = _location(constraint)
        damping, stiffness = _pto(case.get("pto"), "pitch")
        height = _number(wave.get("height"), "wave.height", positive=True)
        period = _number(wave.get("period"), "wave.period", positive=True)
        if "direction" in wave:
            raise ValueError("PM waves use directions and spreading arrays")
        if "seed" in wave and "phase_file" in wave:
            raise ValueError("supply either wave.seed or wave.phase_file")
        if "phase_file" in wave:
            if not isinstance(wave["phase_file"], str) or not wave["phase_file"]:
                raise ValueError("wave.phase_file must be a file path")
            phase_path = (base / wave["phase_file"]).resolve(strict=True)
            phase = np.loadtxt(phase_path, delimiter=",")
            seed = None
            auxiliary_files = (phase_path,)
        else:
            phase = None
            seed = wave.get("seed", 7)
            if not isinstance(seed, int) or isinstance(seed, bool):
                raise ValueError("wave.seed must be an integer")
            auxiliary_files = ()
        count = wave.get("frequency_count", 500)
        components = pm_equal_energy_components(
            hydro[0], significant_height=height, peak_period=period,
            directions=wave.get("directions", [0, 30, 90]),
            spreading=wave.get("spreading", [0.1, 0.2, 0.7]),
            count=count, seed=seed, phase=phase,
        )
        incident = synthesize_irregular_response(
            hydro[0], components, dt=dt, end_time=end_time,
            ramp_time=ramp_time, rho=rho, g=g,
        )
        solved = solve_hinged_pitch_from_excitation(
            hydro[0], incident.excitation_force,
            hinge_z=location[2], body_mass=mass, pitch_inertia=inertia,
            pto_damping=damping, pto_stiffness=stiffness, dt=dt,
            memory_time=_number(
                sim.get("radiation_memory", 30), "simulation.radiation_memory",
                positive=True,
            ),
            rho=rho, g=g,
        )
        position = np.zeros((len(solved.time), 1, 6))
        velocity = np.zeros_like(position)
        position[:, 0, :3] = solved.center_position
        position[:, 0, 4] = solved.angle
        velocity[:, 0, :3] = solved.center_velocity
        velocity[:, 0, 4] = solved.angular_velocity
        return CaseResponse(
            solved.time, position, velocity, hydro,
            pto_force=solved.pto_torque, pto_label="pto_pitch_torque",
            wave_elevation=incident.elevation,
            auxiliary_files=auxiliary_files,
        )

    if kind == "floating_joint":
        if len(bodies) != 2 or wave["type"] != "regular":
            raise ValueError("floating joint needs two bodies and regular waves")
        if hydro[0] != hydro[1]:
            raise ValueError("the current floating-joint layout needs one shared HDF5")
        for number, body in enumerate(bodies, start=1):
            _body_number(body, number)
            if body.get("mass", "equilibrium") != "equilibrium":
                raise ValueError("floating-joint bodies currently require equilibrium mass")
        if "initial_displacement" in constraint:
            raise ValueError("floating-joint initial displacement is not yet supported")
        location = _location(constraint)
        damping, stiffness = _pto(case.get("pto"), "relative_heave")
        height = _number(wave.get("height"), "wave.height", nonnegative=True)
        period = _number(wave.get("period"), "wave.period", positive=True)
        if set(wave) - {"type", "height", "period", "direction"}:
            raise ValueError("regular waves use height, period, and direction")
        direction = _number(wave.get("direction", 0), "wave.direction")
        if direction != 0:
            raise ValueError("floating-joint dynamics currently support 0-degree waves")
        if "radiation_memory" in sim:
            raise ValueError("regular-wave floating-joint dynamics use constant radiation")
        inertias = tuple(
            _number(body.get("pitch_inertia"), "body.pitch_inertia", positive=True)
            for body in bodies
        )
        solved = solve_rm3_regular(
            hydro[0], wave_height=height, wave_period=period,
            pitch_inertias=inertias, pto_damping=damping,
            pto_stiffness=stiffness, b2b=b2b, joint_z=location[2],
            dt=dt, end_time=end_time, ramp_time=ramp_time, rho=rho, g=g,
        )
        ramp = np.ones(len(solved.time))
        if ramp_time > 0:
            early = solved.time < ramp_time
            ramp[early] = (1 - np.cos(np.pi * solved.time[early] / ramp_time)) / 2
        elevation = height / 2 * ramp * np.cos(2 * np.pi * solved.time / period)
        return CaseResponse(
            solved.time, solved.body_position, solved.body_velocity, hydro,
            pto_force=solved.pto_force, pto_label="pto_relative_heave_force",
            wave_elevation=elevation,
        )

    raise ValueError(f"unsupported constraint layout: {kind}")
