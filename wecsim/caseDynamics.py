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

from .bodyClass import BodyClass
from .generalDynamics import BodyMotion, DynamicBody, GeneralizedDynamics
from .hingePitch import (
    solve_hinged_pitch_from_excitation, solve_hinged_pitch_regular,
)
from .irregularWave import pm_equal_energy_components, synthesize_irregular_response
from .linearCoordinates import build_coordinate_maps, initial_coordinate
from .linearHeave import solve_heave_free_decay
from .ptoConnections import build_linear_ptos
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
    pto_generalized_force: np.ndarray | None = None
    extra_outputs: tuple[tuple[str, np.ndarray], ...] = ()


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
    if isinstance(body, Mapping) and body.get("nonhydro", False):
        _section(body, "fixed nonhydrodynamic body",
                 {"nonhydro", "fixed", "center_gravity"},
                 {"nonhydro", "fixed", "center_gravity", "name", "mass", "inertia"})
        if body["nonhydro"] is not True or body["fixed"] is not True:
            raise ValueError("nonhydrodynamic body must be fixed in this layout")
        center = np.asarray(body["center_gravity"], dtype=float)
        if center.shape != (3,) or not np.isfinite(center).all():
            raise ValueError("fixed body center_gravity must be three finite coordinates")
        if "mass" in body:
            _number(body["mass"], "fixed body mass", positive=True)
        if "inertia" in body:
            inertia = np.asarray(body["inertia"], dtype=float)
            if inertia.shape != (3,) or not np.isfinite(inertia).all() or np.any(inertia <= 0):
                raise ValueError("fixed body inertia must have three positive values")
        return None
    _section(body, "body", {"hydro_file"},
             {"hydro_file", "hydro_body", "mass", "pitch_inertia",
              "inertia", "coordinate_map", "name"})
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


def _location(constraint, name="constraint.location"):
    location = np.asarray(constraint.get("location", [0, 0, 0]), dtype=float)
    if location.shape != (3,) or not np.isfinite(location).all():
        raise ValueError(f"{name} must have three finite coordinates")
    if not np.isclose(location[:2], 0, atol=1e-10).all():
        raise ValueError("supported joints must lie on the body x=y=0 axis")
    return location


def _pto(pto, kind, *, allow_location=False):
    _section(pto, "pto", {"kind", "damping"},
             {"kind", "damping", "stiffness", "equilibrium_position",
              "pretension"} | ({"location"} if allow_location else set()))
    if pto["kind"] != kind:
        raise ValueError(f"this layout requires a {kind} PTO")
    if "equilibrium_position" in pto and "pretension" in pto:
        raise ValueError("specify either PTO equilibrium_position or pretension")
    damping = _number(pto["damping"], "pto.damping", nonnegative=True)
    stiffness = _number(pto.get("stiffness", 0), "pto.stiffness", nonnegative=True)
    if "pretension" in pto:
        pretension = _number(pto["pretension"], "pto.pretension")
        if pretension and not stiffness:
            raise ValueError("nonzero PTO pretension needs positive stiffness")
        equilibrium = -pretension / stiffness if stiffness else 0.0
    else:
        equilibrium = _number(pto.get("equilibrium_position", 0),
                              "pto.equilibrium_position")
    if not np.isfinite(equilibrium):
        raise ValueError("PTO equilibrium position is outside its supported range")
    if equilibrium and not stiffness:
        raise ValueError("nonzero PTO equilibrium_position needs positive stiffness")
    return damping, stiffness, equilibrium


def run_case(case: Mapping, *, base_dir: str | Path = ".") -> CaseResponse:
    """Run one explicitly supported wave/device configuration.

    ``base_dir`` resolves relative hydro and phase file paths. The accepted
    layouts are one-body heave free decay, one-body fixed-hinge pitch, a
    two-body floating joint, and mapped linear coordinates for compatible
    bodies and PTO matrices.
    Every one of these layouts uses the generalized dynamics engine.
    """
    case = _section(
        case, "case", {"simulation", "wave", "bodies", "constraint"},
        {"name", "simulation", "wave", "bodies", "constraint", "pto",
         "ptos", "body_to_body"},
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
                          {"kind", "location", "initial_displacement",
                           "initial_coordinate", "initial_speed", "coordinates"})
    bodies = case["bodies"]
    if not isinstance(bodies, list) or not bodies:
        raise ValueError("bodies must be a nonempty list")
    base = Path(base_dir).expanduser().resolve()
    hydro = tuple(_hydro_file(body, base) for body in bodies)
    b2b = case.get("body_to_body", False)
    if not isinstance(b2b, bool):
        raise ValueError("body_to_body must be a boolean")
    kind = constraint["kind"]
    if any(path is None for path in hydro) and not (
        kind == "fixed_hinge" and len(bodies) == 2
        and hydro[0] is not None and hydro[1] is None
    ):
        raise ValueError("a fixed nonhydrodynamic body requires a two-body fixed_hinge")
    if "ptos" in case and kind != "linear_subspace":
        raise ValueError("configurable PTO connections require linear_subspace")

    if kind == "linear_subspace":
        return _run_linear_subspace(
            case, sim, wave, constraint, bodies, hydro, b2b,
            dt, end_time, ramp_time, rho, g,
        )

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
        if (len(bodies) not in (1, 2) or wave["type"] not in ("pm", "regular")
                or (len(bodies) == 2 and hydro[1] is not None) or b2b):
            raise ValueError("fixed-hinge pitch needs one flap, optional fixed base, and PM or regular waves")
        if len(bodies) == 2 and wave["type"] != "regular":
            raise ValueError("the fixed nonhydrodynamic base currently needs regular waves")
        _body_number(bodies[0], 1)
        if set(bodies[0]) - {"hydro_file", "hydro_body", "mass", "pitch_inertia", "name"}:
            raise ValueError("fixed-hinge pitch uses mass and pitch_inertia")
        mass = _number(bodies[0].get("mass"), "body.mass", positive=True)
        inertia = _number(bodies[0].get("pitch_inertia"),
                          "body.pitch_inertia", positive=True)
        if set(constraint) - {"kind", "location"}:
            raise ValueError("fixed-hinge initial conditions are not yet supported")
        location = _location(constraint)
        pto_data = case.get("pto")
        damping, stiffness, equilibrium = _pto(
            pto_data, "pitch", allow_location=True,
        )
        if len(bodies) == 2 and "location" not in pto_data:
            raise ValueError("fixed nonhydrodynamic base needs pto.location")
        hinge = (_location({"location": pto_data["location"]}, "pto.location")
                 if "location" in pto_data else location)
        height = _number(wave.get("height"), "wave.height", positive=True)
        period = _number(wave.get("period"), "wave.period", positive=True)
        if wave["type"] == "regular":
            if set(wave) - {"type", "height", "period", "direction"}:
                raise ValueError("regular waves use height, period, and direction")
            direction = _number(wave.get("direction", 0), "wave.direction")
            if direction != 0:
                raise ValueError("fixed-hinge regular dynamics currently support 0-degree waves")
            if "radiation_memory" in sim:
                raise ValueError("regular fixed-hinge dynamics use constant-frequency radiation")
            solved = solve_hinged_pitch_regular(
                hydro[0], wave_height=height, wave_period=period,
                hinge_z=hinge[2], body_mass=mass, pitch_inertia=inertia,
                pto_damping=damping, pto_stiffness=stiffness,
                pto_equilibrium=equilibrium, dt=dt, end_time=end_time,
                ramp_time=ramp_time, rho=rho, g=g,
            )
            ramp = np.ones(len(solved.time))
            if ramp_time > 0:
                early = solved.time < ramp_time
                ramp[early] = (1 - np.cos(np.pi * solved.time[early] / ramp_time)) / 2
            elevation = height / 2 * ramp * np.cos(2 * np.pi * solved.time / period)
            position = np.zeros((len(solved.time), len(bodies), 6))
            velocity = np.zeros_like(position)
            position[:, 0, :3] = solved.center_position
            position[:, 0, 4] = solved.angle
            velocity[:, 0, :3] = solved.center_velocity
            velocity[:, 0, 4] = solved.angular_velocity
            if len(bodies) == 2:
                position[:, 1, :3] = np.asarray(bodies[1]["center_gravity"], dtype=float)
            return CaseResponse(
                solved.time, position, velocity, (hydro[0],),
                pto_force=solved.pto_torque, pto_label="pto_pitch_torque",
                wave_elevation=elevation,
            )
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
            hinge_z=hinge[2], body_mass=mass, pitch_inertia=inertia,
            pto_damping=damping, pto_stiffness=stiffness,
            pto_equilibrium=equilibrium, dt=dt,
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
        if len(bodies) != 2 or wave["type"] not in ("regular", "regularCIC"):
            raise ValueError("floating joint needs two bodies and regular waves")
        if hydro[0] != hydro[1]:
            raise ValueError("the current floating-joint layout needs one shared HDF5")
        for number, body in enumerate(bodies, start=1):
            _body_number(body, number)
            if set(body) - {"hydro_file", "hydro_body", "mass", "pitch_inertia"}:
                raise ValueError("floating-joint bodies use mass and pitch_inertia")
            if body.get("mass", "equilibrium") != "equilibrium":
                raise ValueError("floating-joint bodies currently require equilibrium mass")
        if set(constraint) - {"kind", "location"}:
            raise ValueError("floating-joint initial conditions are not yet supported")
        location = _location(constraint)
        damping, stiffness, equilibrium = _pto(case.get("pto"), "relative_heave")
        height = _number(wave.get("height"), "wave.height", nonnegative=True)
        period = _number(wave.get("period"), "wave.period", positive=True)
        if set(wave) - {"type", "height", "period", "direction"}:
            raise ValueError("regular waves use height, period, and direction")
        direction = _number(wave.get("direction", 0), "wave.direction")
        if direction != 0:
            raise ValueError("floating-joint dynamics currently support 0-degree waves")
        if wave["type"] == "regular" and "radiation_memory" in sim:
            raise ValueError("regular-wave floating-joint dynamics use constant radiation")
        if wave["type"] == "regularCIC":
            radiation_memory = _number(
                sim.get("radiation_memory", 60), "simulation.radiation_memory",
                positive=True,
            )
        else:
            radiation_memory = None
        inertias = tuple(
            _number(body.get("pitch_inertia"), "body.pitch_inertia", positive=True)
            for body in bodies
        )
        solved = solve_rm3_regular(
            hydro[0], wave_height=height, wave_period=period,
            pitch_inertias=inertias, pto_damping=damping,
            pto_stiffness=stiffness, pto_equilibrium=equilibrium,
            b2b=b2b, radiation_memory=radiation_memory,
            joint_z=location[2],
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


def _run_linear_subspace(case, sim, wave, constraint, bodies, hydro,
                         b2b, dt, end_time, ramp_time, rho, g):
    """Run constant linear coordinate maps with regular or no incident waves."""
    if set(constraint) - {"kind", "initial_coordinate", "initial_speed",
                           "coordinates"}:
        raise ValueError("linear_subspace uses coordinate maps, not joint locations")
    if wave["type"] not in ("regular", "none"):
        raise ValueError("linear_subspace supports regular or no incident waves")
    if b2b and len(set(hydro)) != 1:
        raise ValueError("body-to-body hydrodynamics need one shared HDF5 file")
    if wave["type"] == "none" and set(wave) != {"type"}:
        raise ValueError("no-wave cases have no wave height or period")
    if wave["type"] == "regular":
        if set(wave) - {"type", "height", "period", "direction"}:
            raise ValueError("regular waves use height, period, and direction")
        height = _number(wave.get("height"), "wave.height", nonnegative=True)
        period = _number(wave.get("period"), "wave.period", positive=True)
        direction = _number(wave.get("direction", 0), "wave.direction")
        frequency = 2 * np.pi / period
        if "radiation_memory" in sim:
            raise ValueError("regular-wave linear dynamics use constant radiation")
    else:
        if "ramp_time" in sim:
            raise ValueError("ramp_time is inapplicable to no-wave dynamics")
        memory_time = _number(
            sim.get("radiation_memory", 15), "simulation.radiation_memory",
            positive=True,
        )
        memory_steps = round(memory_time / dt)
        if not np.isclose(memory_steps * dt, memory_time, atol=1e-10):
            raise ValueError("radiation_memory must be a multiple of dt")
        convolution_time = np.arange(memory_steps + 1) * dt
    time = np.arange(round(end_time / dt) + 1) * dt
    if not np.isclose(time[-1], end_time, atol=1e-10):
        raise ValueError("end_time must be an integer multiple of dt")
    body_names = []
    for index, body_spec in enumerate(bodies, start=1):
        _body_number(body_spec, index)
        if "pitch_inertia" in body_spec:
            raise ValueError("linear_subspace uses body.inertia, not pitch_inertia")
        name = body_spec.get("name", f"body{index}")
        if not isinstance(name, str) or not name or name in body_names:
            raise ValueError("body names must be nonempty and unique")
        body_names.append(name)
    loaded_bodies = []
    centers = []
    for index, path in enumerate(hydro, start=1):
        body = BodyClass(str(path))
        body.bodyNumber = index
        body.bodyTotal = len(bodies)
        body.readH5file()
        if int(np.asarray(body.dof).item()) != 6:
            raise ValueError("linear_subspace needs six-DOF hydrodynamic bodies")
        center = np.asarray(body.cg, dtype=float).ravel()
        if center.shape != (3,):
            raise ValueError("body center must have three coordinates")
        loaded_bodies.append(body)
        centers.append(center)
    maps, coordinate_names = build_coordinate_maps(
        constraint, bodies, body_names, centers,
    )
    n = maps[0].shape[1]
    initial_q = initial_coordinate(
        constraint.get("initial_coordinate", [0] * n), coordinate_names,
        "constraint.initial_coordinate",
    )
    initial_v = initial_coordinate(
        constraint.get("initial_speed", [0] * n), coordinate_names,
        "constraint.initial_speed",
    )
    if "pto" in case and "ptos" in case:
        raise ValueError("use either a PTO matrix or PTO connections")
    if "pto" in case:
        pto = _section(case["pto"], "pto", {"kind"},
                       {"kind", "stiffness_matrix", "damping_matrix",
                        "equilibrium_coordinate"})
        if pto["kind"] != "linear":
            raise ValueError("linear_subspace requires a linear PTO matrix")
        stiffness = np.asarray(pto.get("stiffness_matrix", np.zeros((n, n))), dtype=float)
        damping = np.asarray(pto.get("damping_matrix", np.zeros((n, n))), dtype=float)
        if (stiffness.shape != (n, n) or damping.shape != (n, n)
                or not np.isfinite(stiffness).all() or not np.isfinite(damping).all()):
            raise ValueError("PTO matrices must be finite N-by-N matrices")
        equilibrium = np.asarray(pto.get("equilibrium_coordinate", [0] * n),
                                 dtype=float)
        if equilibrium.shape != (n,) or not np.isfinite(equilibrium).all():
            raise ValueError("pto.equilibrium_coordinate must have N finite values")
        if np.any(equilibrium) and not np.any(stiffness @ equilibrium):
            raise ValueError("PTO equilibrium offset must produce a spring force")
    else:
        stiffness = np.zeros((n, n))
        damping = np.zeros((n, n))
        equilibrium = np.zeros(n)
    pto_bias = np.zeros(n)
    connections = ()

    dynamic_bodies = []
    for index, (body_spec, body, mapping) in enumerate(
            zip(bodies, loaded_bodies, maps), start=1):
        mass_setting = body_spec.get("mass", "equilibrium")
        if mass_setting != "equilibrium":
            mass_setting = _number(mass_setting, "body.mass", positive=True)
        body.mass = mass_setting
        inertia = np.asarray(body_spec.get("inertia", [0, 0, 0]), dtype=float)
        if (inertia.shape != (3,) or not np.isfinite(inertia).all()
                or np.any(inertia < 0)):
            raise ValueError("body.inertia must have three nonnegative entries")
        body.hydroStiffness = np.zeros((6, 6))
        body.viscDrag = {
            "Drag": np.zeros((6, 6)), "cd": np.zeros(6),
            "characteristicArea": np.zeros(6),
        }
        body.linearDamping = np.zeros((6, 6))
        wave_amp = np.vstack((time, np.zeros_like(time)))
        if wave["type"] == "regular":
            body.hydroForcePre(
                frequency, [direction], 1, np.array([0.0]), [], dt, rho, g,
                "regular", wave_amp, index, len(bodies), 0, 0, int(b2b),
            )
        else:
            irf_time = body.hydroData["hydro_coeffs"]["radiation_damping"][
                "impulse_response_fun"]["t"]
            if memory_time > np.max(irf_time) + 1e-10:
                raise ValueError("radiation_memory exceeds the HDF5 kernel")
            body.hydroForcePre(
                [], [0], len(convolution_time), convolution_time, [],
                dt, rho, g, "noWaveCIC", wave_amp,
                index, len(bodies), 0, 0, int(b2b),
            )
        physical_mass = float(np.asarray(body.mass).item())
        rigid_mass = np.diag([physical_mass] * 3 + inertia.tolist())
        hydro_force = body.hydroForce
        if b2b:
            added_mass = tuple(
                np.asarray(hydro_force["fAddedMass"])[:, 6*j:6*(j+1)]
                for j in range(len(bodies))
            )
            if wave["type"] == "regular":
                radiation_damping = tuple(
                    np.asarray(hydro_force["fDamping"])[:, 6*j:6*(j+1)]
                    for j in range(len(bodies))
                )
        else:
            added_mass = [np.zeros((6, 6)) for _ in bodies]
            added_mass[index - 1] = np.asarray(hydro_force["fAddedMass"])
            if wave["type"] == "regular":
                radiation_damping = [np.zeros((6, 6)) for _ in bodies]
                radiation_damping[index - 1] = np.asarray(hydro_force["fDamping"])
        if wave["type"] == "none":
            radiation_damping = tuple(np.zeros((6, 6)) for _ in bodies)
            kernel = np.asarray(hydro_force["irkb"])
            if not b2b and len(bodies) > 1:
                independent_kernel = np.zeros((len(kernel), 6, 6 * len(bodies)))
                independent_kernel[:, :, 6 * (index - 1):6 * index] = kernel
                kernel = independent_kernel
        else:
            kernel = None
        center = centers[index - 1]

        def motion(q, v, *, mapping=mapping):
            return BodyMotion(mapping @ q, mapping, np.zeros(6))

        if wave["type"] == "regular":
            real = np.asarray(hydro_force["fExt"]["re"])
            imaginary = np.asarray(hydro_force["fExt"]["im"])

            def excitation(at_time, *, real=real, imaginary=imaginary):
                ramp = (1.0 if ramp_time == 0 or at_time >= ramp_time
                        else (1 - np.cos(np.pi * at_time / ramp_time)) / 2)
                return height / 2 * ramp * (
                    real * np.cos(frequency * at_time)
                    - imaginary * np.sin(frequency * at_time)
                )
        else:
            def excitation(at_time):
                return np.zeros(6)

        dynamic_bodies.append(DynamicBody(
            rigid_mass=rigid_mass,
            added_mass=tuple(added_mass),
            damping=tuple(radiation_damping),
            restoring=np.asarray(hydro_force["linearHydroRestCoef"]),
            static_force=np.array([
                0, 0, (rho * float(np.asarray(body.dispVol).item()) - physical_mass) * g,
                0, 0, 0,
            ]),
            reference_position=np.r_[center, np.zeros(3)],
            motion=motion,
            excitation=excitation,
            radiation_kernel=kernel,
        ))
    if "ptos" in case:
        connections, stiffness, damping, pto_bias = build_linear_ptos(
            case["ptos"], maps, centers, body_names,
        )
    system = GeneralizedDynamics(
        tuple(dynamic_bodies), n, pto_stiffness=stiffness,
        pto_damping=damping, pto_equilibrium=equilibrium,
        pto_bias=pto_bias,
    )
    response = system.integrate(
        dt=dt, end_time=end_time,
        initial_coordinate=initial_q, initial_speed=initial_v,
    )
    if wave["type"] == "regular":
        ramp = np.ones(len(response.time))
        if ramp_time > 0:
            early = response.time < ramp_time
            ramp[early] = (1 - np.cos(np.pi * response.time[early] / ramp_time)) / 2
        elevation = height / 2 * ramp * np.cos(frequency * response.time)
    else:
        elevation = None
    coordinate_outputs = tuple(
        output for index, name in enumerate(coordinate_names)
        for output in (
            (f"coordinate_{name}_position", response.coordinate[:, index]),
            (f"coordinate_{name}_velocity", response.speed[:, index]),
        )
    ) if "coordinates" in constraint else ()
    pto_outputs = tuple(
        output for connection in connections
        for output in connection.outputs(response.coordinate, response.speed)
    )
    return CaseResponse(
        response.time, response.body_position, response.body_velocity,
        hydro, wave_elevation=elevation,
        pto_generalized_force=(
            -(response.coordinate - equilibrium) @ stiffness.T
            - response.speed @ damping.T + pto_bias
            if "pto" in case or "ptos" in case else None
        ),
        extra_outputs=coordinate_outputs + pto_outputs,
    )
