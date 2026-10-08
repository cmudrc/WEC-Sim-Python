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
from scipy.io import loadmat

from .bodyClass import BodyClass
from .generalDynamics import BodyMotion, DynamicBody, GeneralizedDynamics
from .hardStops import LinearHardStops
from .hingePitch import (
    solve_hinged_pitch_from_excitation, solve_hinged_pitch_regular,
)
from .irregularWave import (
    imported_full_directional_components, pm_equal_energy_components,
    synthesize_full_directional_response, synthesize_irregular_response,
    synthesize_multiple_irregular_response,
)
from .linearCoordinates import build_coordinate_maps, initial_coordinate
from .linearHeave import solve_heave_free_decay
from .passiveYaw import PassiveYawExcitation
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
              "inertia", "coordinate_map", "name", "mean_drift", "fixed",
              "passive_yaw"})
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


def _hard_stops(spec):
    required = {"lower_bound", "upper_bound", "lower_stiffness", "upper_stiffness"}
    optional = {"lower_damping", "upper_damping", "lower_transition_width",
                "upper_transition_width"}
    _section(spec, "pto.hard_stops", required, required | optional)
    values = {
        name: _number(value, f"pto.hard_stops.{name}",
                      positive=name.endswith("stiffness") or name.endswith("width"),
                      nonnegative=name.endswith("damping"))
        for name, value in spec.items()
    }
    return LinearHardStops(**values)


def _imported_elevation(wave, base, hydro_file, time, dt, ramp_time, rho, g):
    """Build RM3 body forces from one sampled MATLAB elevation record."""
    _section(wave, "imported wave", {"type", "file"},
             {"type", "file", "variable", "direction", "reapply_force_ramp"})
    raw = wave["file"]
    if not isinstance(raw, (str, Path)) or not str(raw):
        raise ValueError("wave.file must be a MAT file path")
    path = (base / raw).resolve(strict=True)
    if not path.is_file():
        raise ValueError(f"wave.file is not a file: {path}")
    variable = wave.get("variable", "etaData")
    if not isinstance(variable, str) or not variable:
        raise ValueError("wave.variable must be a nonempty MAT variable name")
    direction = _number(wave.get("direction", 0), "wave.direction")
    if direction != 0:
        raise ValueError("floating-joint imported elevation currently supports 0-degree waves")
    second_ramp = wave.get("reapply_force_ramp", False)
    if not isinstance(second_ramp, bool):
        raise ValueError("wave.reapply_force_ramp must be a boolean")
    mat = loadmat(path, variable_names=[variable])
    if variable not in mat:
        raise ValueError(f"wave.variable {variable!r} is absent from {path}")
    samples = np.asarray(mat[variable])
    if not np.issubdtype(samples.dtype, np.number) or not np.isrealobj(samples):
        raise ValueError("imported elevation must contain real numeric samples")
    samples = samples.astype(float, copy=False)
    if (samples.ndim != 2 or samples.shape[1] != 2 or samples.shape[0] < 2
            or not np.isfinite(samples).all()
            or not np.all(np.diff(samples[:, 0]) > 0)):
        raise ValueError("imported elevation must be a finite, increasing N-by-2 time/elevation array")
    if samples[0, 0] > time[0] + 1e-10 or samples[-1, 0] < time[-1] - 1e-10:
        raise ValueError("imported elevation must cover the full simulation time")
    ramp = np.ones_like(time)
    if ramp_time > 0:
        early = time < ramp_time
        ramp[early] = (1 - np.cos(np.pi * time[early] / ramp_time)) / 2
    elevation = np.interp(time, samples[:, 0], samples[:, 1]) * ramp
    force = np.zeros((len(time), 2, 6))
    for number in (1, 2):
        body = BodyClass(str(hydro_file))
        body.bodyNumber = number
        body.bodyTotal = 2
        body.readH5file()
        body.hydroForce["userDefinedFe"] = np.zeros((len(time), 6))
        body.userDefinedExcitation(np.vstack((time, elevation)), dt, [direction], rho, g)
        force[:, number - 1] = body.hydroForce["userDefinedFe"]
    if second_ramp:
        # The pinned MATLAB body block ramps force after waveClass ramped elevation.
        force *= ramp[:, None, None]
    return elevation, force, path


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
         "ptos", "body_to_body", "mooring"},
    )
    sim = _section(case["simulation"], "simulation", {"dt", "end_time"},
                   {"dt", "end_time", "ramp_time", "rho", "g",
                    "radiation_memory", "radiation_method", "added_mass_scheme"})
    dt = _number(sim["dt"], "simulation.dt", positive=True)
    end_time = _number(sim["end_time"], "simulation.end_time", nonnegative=True)
    ramp_time = _number(sim.get("ramp_time", 100), "simulation.ramp_time", nonnegative=True)
    rho = _number(sim.get("rho", 1000), "simulation.rho", positive=True)
    g = _number(sim.get("g", 9.81), "simulation.g", positive=True)
    wave = _section(case["wave"], "wave", {"type"},
                    {"type", "height", "period", "direction", "directions",
                     "spreading", "seed", "phase_file", "frequency_count",
                     "file", "variable", "reapply_force_ramp", "seas",
                     "excitation_interpolation", "force_quadrature"})
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
    for index, body in enumerate(bodies):
        if "fixed" in body and not (kind == "fixed_hinge" and len(bodies) == 2
                                     and index == 1 and body["fixed"] is True):
            raise ValueError("a fixed hydrodynamic body currently requires a two-body fixed_hinge")
    if "radiation_method" in sim and kind != "floating_joint":
        raise ValueError("radiation_method currently applies to floating_joint")
    if "added_mass_scheme" in sim and kind not in ("floating_joint", "fixed_hinge"):
        raise ValueError("added_mass_scheme currently applies to floating_joint or fixed_hinge")
    if any(path is None for path in hydro) and not (
        kind == "fixed_hinge" and len(bodies) == 2
        and hydro[0] is not None and hydro[1] is None
    ):
        raise ValueError("a fixed nonhydrodynamic body requires a two-body fixed_hinge")
    if "ptos" in case and kind != "linear_subspace":
        raise ValueError("configurable PTO connections require linear_subspace")
    if kind != "linear_subspace" and any("mean_drift" in body for body in bodies):
        raise ValueError("body.mean_drift currently requires linear_subspace")
    if "mooring" in case and kind != "floating_joint":
        raise ValueError("the joint surge mooring requires a floating_joint")

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
        if (len(bodies) not in (1, 2)
                or wave["type"] not in ("pm", "pm_multi", "regular",
                                        "spectrumImportFullDir")
                or b2b):
            raise ValueError("fixed-hinge pitch needs one flap, optional fixed base, and PM or regular waves")
        if len(bodies) == 2 and hydro[1] is not None:
            _body_number(bodies[1], 2)
            fixed_body = BodyClass(str(hydro[1]))
            fixed_body.bodyNumber = 2
            fixed_body.readH5file()
            fixed_center = np.asarray(fixed_body.cg, dtype=float).ravel()
        elif len(bodies) == 2:
            fixed_center = np.asarray(bodies[1]["center_gravity"], dtype=float)
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
        if wave["type"] == "regular":
            if sim.get("added_mass_scheme", "implicit") != "implicit":
                raise ValueError("regular fixed-hinge dynamics need implicit added mass")
            height = _number(wave.get("height"), "wave.height", positive=True)
            period = _number(wave.get("period"), "wave.period", positive=True)
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
                position[:, 1, :3] = fixed_center
            return CaseResponse(
                solved.time, position, velocity,
                tuple(path for path in hydro if path is not None),
                pto_force=solved.pto_torque, pto_label="pto_pitch_torque",
                wave_elevation=elevation,
            )
        if wave["type"] == "spectrumImportFullDir":
            if (set(wave) - {"type", "file", "phase_file", "seed",
                             "excitation_interpolation", "force_quadrature"}
                    or "file" not in wave):
                raise ValueError("spectrumImportFullDir needs a MAT spectrum file")
            if "phase_file" in wave and "seed" in wave:
                raise ValueError("supply either wave.phase_file or wave.seed")
            if not isinstance(wave["file"], str) or not wave["file"]:
                raise ValueError("wave.file must be a MAT file path")
            spectrum_path = (base / wave["file"]).resolve(strict=True)
            auxiliary_files = [spectrum_path]
            if "phase_file" in wave:
                if not isinstance(wave["phase_file"], str) or not wave["phase_file"]:
                    raise ValueError("wave.phase_file must be a file path")
                phase_path = (base / wave["phase_file"]).resolve(strict=True)
                phase = np.loadtxt(phase_path, delimiter=",", ndmin=2)
                auxiliary_files.append(phase_path)
                seed = None
            else:
                phase = None
                seed = wave.get("seed", 7)
                if not isinstance(seed, int) or isinstance(seed, bool):
                    raise ValueError("wave.seed must be an integer")
            components = imported_full_directional_components(
                hydro[0], spectrum_path, phase=phase, seed=seed,
            )
            incident = synthesize_full_directional_response(
                hydro[0], components, dt=dt, end_time=end_time,
                ramp_time=ramp_time, rho=rho, g=g,
                excitation_interpolation=wave.get("excitation_interpolation", "linear"),
                force_quadrature=wave.get("force_quadrature", "integrated"),
            )
        elif wave["type"] == "pm_multi":
            if set(wave) - {"type", "seas", "excitation_interpolation"} or "seas" not in wave:
                raise ValueError("pm_multi waves use a seas list")
            seas = wave["seas"]
            if not isinstance(seas, list) or len(seas) < 2:
                raise ValueError("pm_multi needs at least two sea spectra")
            components = []
            auxiliary_files = []
            for index, sea in enumerate(seas):
                _section(sea, f"wave.seas[{index}]", {"height", "period"},
                         {"height", "period", "direction", "directions",
                          "spreading", "seed", "phase_file", "frequency_count"})
                if "direction" in sea:
                    if "directions" in sea or "spreading" in sea:
                        raise ValueError("a sea uses direction or directions and spreading")
                    directions = [_number(sea["direction"],
                                          f"wave.seas[{index}].direction")]
                    spreading = [1.0]
                else:
                    if "directions" not in sea or "spreading" not in sea:
                        raise ValueError("a sea needs direction or directions and spreading")
                    directions, spreading = sea["directions"], sea["spreading"]
                if "phase_file" in sea and "seed" in sea:
                    raise ValueError("a sea uses either phase_file or seed")
                if "phase_file" in sea:
                    if not isinstance(sea["phase_file"], str) or not sea["phase_file"]:
                        raise ValueError("sea.phase_file must be a file path")
                    phase_path = (base / sea["phase_file"]).resolve(strict=True)
                    phase = np.loadtxt(phase_path, delimiter=",", ndmin=2)
                    seed = None
                    auxiliary_files.append(phase_path)
                else:
                    phase = None
                    seed = sea.get("seed", index + 7)
                    if not isinstance(seed, int) or isinstance(seed, bool):
                        raise ValueError("sea.seed must be an integer")
                components.append(pm_equal_energy_components(
                    hydro[0],
                    significant_height=_number(
                        sea["height"], f"wave.seas[{index}].height", positive=True,
                    ),
                    peak_period=_number(
                        sea["period"], f"wave.seas[{index}].period", positive=True,
                    ),
                    directions=directions, spreading=spreading,
                    count=sea.get("frequency_count", 500), seed=seed, phase=phase,
                ))
            incident = synthesize_multiple_irregular_response(
                hydro[0], components, dt=dt, end_time=end_time,
                ramp_time=ramp_time, rho=rho, g=g,
                excitation_interpolation=wave.get("excitation_interpolation", "linear"),
            )
        else:
            if set(wave) - {"type", "height", "period", "directions",
                             "spreading", "seed", "phase_file",
                             "frequency_count", "excitation_interpolation"}:
                raise ValueError("PM waves use height, period, directions, and phase settings")
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
                phase = np.loadtxt(phase_path, delimiter=",", ndmin=2)
                seed = None
                auxiliary_files = [phase_path]
            else:
                phase = None
                seed = wave.get("seed", 7)
                if not isinstance(seed, int) or isinstance(seed, bool):
                    raise ValueError("wave.seed must be an integer")
                auxiliary_files = []
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
                excitation_interpolation=wave.get("excitation_interpolation", "linear"),
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
            added_mass_scheme=sim.get("added_mass_scheme", "implicit"),
            rho=rho, g=g,
        )
        position = np.zeros((len(solved.time), len(bodies), 6))
        velocity = np.zeros_like(position)
        position[:, 0, :3] = solved.center_position
        position[:, 0, 4] = solved.angle
        velocity[:, 0, :3] = solved.center_velocity
        velocity[:, 0, 4] = solved.angular_velocity
        if len(bodies) == 2:
            position[:, 1, :3] = fixed_center
        return CaseResponse(
            solved.time, position, velocity,
            tuple(path for path in hydro if path is not None),
            pto_force=solved.pto_torque, pto_label="pto_pitch_torque",
            wave_elevation=incident.elevation,
            auxiliary_files=tuple(auxiliary_files),
        )

    if kind == "floating_joint":
        if len(bodies) != 2 or wave["type"] not in (
                "regular", "regularCIC", "none", "elevationImport"):
            raise ValueError("floating joint needs two bodies and regular waves, imported elevation, or no waves")
        if hydro[0] != hydro[1]:
            raise ValueError("the current floating-joint layout needs one shared HDF5")
        for number, body in enumerate(bodies, start=1):
            _body_number(body, number)
            if set(body) - {"hydro_file", "hydro_body", "mass", "pitch_inertia"}:
                raise ValueError("floating-joint bodies use mass and pitch_inertia")
            if body.get("mass", "equilibrium") != "equilibrium":
                raise ValueError("floating-joint bodies currently require equilibrium mass")
        if set(constraint) - {"kind", "location", "initial_coordinate", "initial_speed"}:
            raise ValueError("unsupported floating-joint constraint setting")
        location = _location(constraint)
        coordinate_names = ("surge", "float_heave", "spar_heave", "pitch")
        initial_q = initial_coordinate(
            constraint.get("initial_coordinate", [0] * 4),
            coordinate_names, "constraint.initial_coordinate",
        )
        initial_v = initial_coordinate(
            constraint.get("initial_speed", [0] * 4),
            coordinate_names, "constraint.initial_speed",
        )
        pto_spec = case.get("pto")
        hard_stops = None
        if isinstance(pto_spec, Mapping) and "hard_stops" in pto_spec:
            hard_stops = _hard_stops(pto_spec["hard_stops"])
            pto_spec = {key: value for key, value in pto_spec.items()
                        if key != "hard_stops"}
        damping, stiffness, equilibrium = _pto(pto_spec, "relative_heave")
        mooring_stiffness = 0.0
        if "mooring" in case:
            mooring = _section(case["mooring"], "mooring", {"kind", "stiffness"},
                               {"kind", "stiffness"})
            if mooring["kind"] != "joint_surge_spring":
                raise ValueError("floating_joint currently supports joint_surge_spring mooring")
            mooring_stiffness = _number(mooring["stiffness"],
                                        "mooring.stiffness", positive=True)
        imported_force = None
        auxiliary_files = ()
        if wave["type"] == "none":
            if set(wave) != {"type"} or "ramp_time" in sim:
                raise ValueError("no-wave floating joint has no wave or ramp settings")
            height, period, direction = 0.0, 8.0, 0.0
        elif wave["type"] == "elevationImport":
            steps = round(end_time / dt)
            if not np.isclose(steps * dt, end_time, rtol=0, atol=1e-10):
                raise ValueError("end_time must be an integer multiple of dt")
            time = np.arange(steps + 1) * dt
            elevation, imported_force, wave_path = _imported_elevation(
                wave, base, hydro[0], time, dt, ramp_time, rho, g,
            )
            auxiliary_files = (wave_path,)
            height, period, direction = 0.0, 8.0, 0.0
        else:
            height = _number(wave.get("height"), "wave.height", nonnegative=True)
            period = _number(wave.get("period"), "wave.period", positive=True)
            if set(wave) - {"type", "height", "period", "direction"}:
                raise ValueError("regular waves use height, period, and direction")
            direction = _number(wave.get("direction", 0), "wave.direction")
            if direction != 0:
                raise ValueError("floating-joint dynamics currently support 0-degree waves")
        if wave["type"] == "regular" and "radiation_memory" in sim:
            raise ValueError("regular-wave floating-joint dynamics use constant radiation")
        radiation_method = sim.get(
            "radiation_method",
            "constant" if wave["type"] == "regular" else "convolution",
        )
        if wave["type"] == "regular" and radiation_method != "constant":
            raise ValueError("regular waves need constant radiation")
        if (wave["type"] in ("regularCIC", "none", "elevationImport")
                and radiation_method not in ("convolution", "fir")):
            raise ValueError("radiation memory needs convolution or FIR radiation")
        if wave["type"] == "elevationImport" and radiation_method != "convolution":
            raise ValueError("imported elevation currently requires convolution radiation")
        if wave["type"] in ("regularCIC", "none", "elevationImport"):
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
            pto_hard_stops=hard_stops,
            mooring_surge_stiffness=mooring_stiffness,
            b2b=b2b, radiation_memory=radiation_memory,
            radiation_method=radiation_method,
            added_mass_scheme=sim.get("added_mass_scheme", "implicit"),
            no_wave=wave["type"] == "none",
            excitation_force=imported_force,
            initial_coordinate=initial_q, initial_speed=initial_v,
            joint_z=location[2],
            dt=dt, end_time=end_time, ramp_time=ramp_time, rho=rho, g=g,
        )
        if wave["type"] == "none":
            elevation = None
        elif wave["type"] != "elevationImport":
            ramp = np.ones(len(solved.time))
            if ramp_time > 0:
                early = solved.time < ramp_time
                ramp[early] = (1 - np.cos(np.pi * solved.time[early] / ramp_time)) / 2
            elevation = height / 2 * ramp * np.cos(2 * np.pi * solved.time / period)
        extra_outputs = []
        if hard_stops is not None:
            extra_outputs.append(("pto_stop_force", solved.pto_stop_force))
        if mooring_stiffness:
            extra_outputs.extend((
                ("mooring_surge_position", solved.mooring_surge_position),
                ("mooring_surge_force", solved.mooring_surge_force),
            ))
        return CaseResponse(
            solved.time, solved.body_position, solved.body_velocity, hydro,
            pto_force=solved.pto_force, pto_label="pto_relative_heave_force",
            wave_elevation=elevation,
            auxiliary_files=auxiliary_files,
            extra_outputs=tuple(extra_outputs),
        )

    raise ValueError(f"unsupported constraint layout: {kind}")


def _run_linear_subspace(case, sim, wave, constraint, bodies, hydro,
                         b2b, dt, end_time, ramp_time, rho, g):
    """Run constant coordinate maps with regular, regularCIC, or no waves."""
    if set(constraint) - {"kind", "initial_coordinate", "initial_speed",
                           "coordinates"}:
        raise ValueError("linear_subspace uses coordinate maps, not joint locations")
    if wave["type"] not in ("regular", "regularCIC", "none"):
        raise ValueError("linear_subspace supports regular, regularCIC, or no waves")
    if b2b and len(set(hydro)) != 1:
        raise ValueError("body-to-body hydrodynamics need one shared HDF5 file")
    if wave["type"] == "none" and set(wave) != {"type"}:
        raise ValueError("no-wave cases have no wave height or period")
    if wave["type"] in ("regular", "regularCIC"):
        if set(wave) - {"type", "height", "period", "direction"}:
            raise ValueError("regular waves use height, period, and direction")
        height = _number(wave.get("height"), "wave.height", nonnegative=True)
        period = _number(wave.get("period"), "wave.period", positive=True)
        direction = _number(wave.get("direction", 0), "wave.direction")
        frequency = 2 * np.pi / period
        if wave["type"] == "regular" and "radiation_memory" in sim:
            raise ValueError("regular-wave linear dynamics use constant radiation")
    else:
        if "ramp_time" in sim:
            raise ValueError("ramp_time is inapplicable to no-wave dynamics")
    if wave["type"] != "regular":
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
        drift_option = bodies[index - 1].get("mean_drift", "none")
        drift_flags = {"none": 0, "control_surface": 1,
                       "momentum_conservation": 2}
        if not isinstance(drift_option, str) or drift_option not in drift_flags:
            raise ValueError(
                "body.mean_drift must be none, control_surface, or momentum_conservation"
            )
        if drift_option != "none" and wave["type"] == "none":
            raise ValueError("mean drift requires regular incident waves")
        body.meanDriftForce = drift_flags[drift_option]
        body.readH5file()
        if drift_option != "none":
            drift_data = body.hydroData["hydro_coeffs"]["mean_drift"]
            if (drift_data.ndim != 3 or drift_data.shape[0] != 6
                    or not np.isfinite(drift_data).all()):
                raise ValueError(
                    f"body{index} HDF5 lacks finite {drift_option} mean-drift coefficients"
                )
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
    passive_indices = [
        index for index, spec in enumerate(bodies)
        if spec.get("passive_yaw", False) is True
    ]
    if any(not isinstance(spec.get("passive_yaw", False), bool)
           for spec in bodies):
        raise ValueError("body.passive_yaw must be a boolean")
    if passive_indices:
        yaw_map = np.zeros((6, 1))
        yaw_map[5, 0] = 1
        if (len(passive_indices) != 1 or n != 1 or b2b
                or wave["type"] != "regular"
                or not np.array_equal(maps[passive_indices[0]], yaw_map)
                or any(np.any(mapping) for index, mapping in enumerate(maps)
                       if index != passive_indices[0])
                or bodies[passive_indices[0]].get("mean_drift", "none") != "none"):
            raise ValueError(
                "passive yaw currently needs one pure-yaw body, stationary others, "
                "regular waves, and independent radiation"
            )
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
    passive_model = None
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
            regular_memory = wave["type"] == "regularCIC"
            body.hydroForcePre(
                frequency if regular_memory else [],
                [direction if regular_memory else 0],
                len(convolution_time), convolution_time, [],
                dt, rho, g, "regularCIC" if regular_memory else "noWaveCIC",
                wave_amp,
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
        if wave["type"] != "regular":
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

        if wave["type"] in ("regular", "regularCIC"):
            real = np.asarray(hydro_force["fExt"]["re"])
            imaginary = np.asarray(hydro_force["fExt"]["im"])
            drift = np.asarray(hydro_force["fExt"]["md"])

            def excitation(at_time, *, real=real, imaginary=imaginary,
                           drift=drift):
                ramp = (1.0 if ramp_time == 0 or at_time >= ramp_time
                        else (1 - np.cos(np.pi * at_time / ramp_time)) / 2)
                return (height / 2 * ramp * (
                    real * np.cos(frequency * at_time)
                    - imaginary * np.sin(frequency * at_time)
                ) + (height / 2)**2 * ramp * drift)
        else:
            def excitation(at_time):
                return np.zeros(6)

        state_excitation = None
        if body_spec.get("passive_yaw", False):
            passive_model = PassiveYawExcitation.from_hydro_data(
                body.hydroData, omega=frequency,
                incident_direction=direction, amplitude=height / 2,
                ramp_time=ramp_time, rho=rho, g=g,
            )

            def state_excitation(at_time, coordinate, speed, *,
                                 model=passive_model):
                return model.force(at_time, coordinate[0])

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
            state_excitation=state_excitation,
        ))
    if "ptos" in case:
        connections, stiffness, damping, pto_bias = build_linear_ptos(
            case["ptos"], maps, centers, body_names, coordinate_names,
        )
    controlled_connections = tuple(
        connection for connection in connections if connection.control is not None
    )
    system = GeneralizedDynamics(
        tuple(dynamic_bodies), n, pto_stiffness=stiffness,
        pto_damping=damping, pto_equilibrium=equilibrium,
        pto_bias=pto_bias,
        controlled_ptos=controlled_connections,
    )
    response = system.integrate(
        dt=dt, end_time=end_time,
        initial_coordinate=initial_q, initial_speed=initial_v,
    )
    if wave["type"] in ("regular", "regularCIC"):
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
    controlled_forces = (
        {connection.name: response.controlled_pto_force[:, index]
         for index, connection in enumerate(controlled_connections)}
        if controlled_connections else {}
    )
    pto_outputs = tuple(
        output for connection in connections
        for output in connection.outputs(
            response.coordinate, response.speed,
            force_override=controlled_forces.get(connection.name),
        )
    )
    generalized_pto = (-(response.coordinate - equilibrium) @ stiffness.T
                       - response.speed @ damping.T + pto_bias)
    for connection in controlled_connections:
        generalized_pto += np.outer(
            controlled_forces[connection.name], connection.stroke_jacobian,
        )
    drift_outputs = tuple(
        output
        for index, (spec, body, dynamic) in enumerate(
            zip(bodies, loaded_bodies, dynamic_bodies), start=1,
        )
        if spec.get("mean_drift", "none") != "none"
        for output in (
            (f"body{index}_mean_drift_force",
             (height / 2)**2 * ramp[:, None]
             * np.asarray(body.hydroForce["fExt"]["md"])[None, :]),
            (f"body{index}_excitation_force",
             np.stack([dynamic.excitation(t) for t in response.time])),
        )
    ) if wave["type"] in ("regular", "regularCIC") else ()
    passive_outputs = ()
    if passive_model is not None:
        body_index = passive_indices[0] + 1
        passive_outputs = ((
            f"body{body_index}_excitation_force",
            np.stack([
                passive_model.force(t, angle)
                for t, angle in zip(response.time, response.coordinate[:, 0])
            ]),
        ),)
    return CaseResponse(
        response.time, response.body_position, response.body_velocity,
        hydro, wave_elevation=elevation,
        pto_generalized_force=(
            generalized_pto if "pto" in case or "ptos" in case else None
        ),
        extra_outputs=(coordinate_outputs + pto_outputs + drift_outputs
                       + passive_outputs),
    )
