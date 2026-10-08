"""Python interface for configuring supported linearized WEC devices.

The objects in this module assemble the same validated case used by the JSON
runner. Users can define bodies, named motions, attachment points, PTOs, and
waves directly in Python and receive NumPy arrays without writing JSON or CSV.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from numbers import Real
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from .caseDynamics import CaseResponse, run_case
from .controls import DeclutchingControl, LatchingControl
from .directLinearGenerator import DirectLinearGenerator
from .hardStops import LinearHardStops
from .morison import MorisonElement


@dataclass(frozen=True)
class WorldPoint:
    """A fixed point in world coordinates, also usable as a rotation pivot."""

    x: float
    y: float
    z: float

    def coordinates(self) -> list[float]:
        return [self.x, self.y, self.z]


@dataclass(frozen=True)
class BodyPoint:
    """A point in a body's reference frame, relative to its center of gravity."""

    body: Body
    x: float
    y: float
    z: float

    def coordinates(self) -> list[float]:
        return [self.x, self.y, self.z]


@dataclass(frozen=True)
class Motion:
    body: Body
    dof: str
    scale: float = 1.0
    pivot: WorldPoint | BodyPoint | None = None


@dataclass(frozen=True)
class HydroState:
    """One BEM dataset and rigid properties for a variable-draft body."""

    hydro_file: str | Path
    mass: str | float = "equilibrium"
    inertia: tuple[float, float, float] = (0, 0, 0)


@dataclass(frozen=True)
class VariableHydro:
    """Ordered states; each switch time activates the next state."""

    states: tuple[HydroState, ...]
    switch_times: tuple[float, ...]


@dataclass(frozen=True)
class Body:
    name: str
    hydro_file: str | Path | None
    mass: str | float = "equilibrium"
    inertia: tuple[float, float, float] = (0, 0, 0)
    hydro_body: int | None = None
    mean_drift: str = "none"
    passive_yaw: bool = False
    geometry_file: str | Path | None = None
    nonlinear_hydro: str | None = None
    drag_coefficient: float = 0.0
    drag_area: float = 0.0
    variable_hydro: VariableHydro | None = None
    passive_yaw_threshold: float = 0.0
    yaw_heading_bank: tuple[float, ...] | None = None
    fixed: bool = False
    center_gravity: tuple[float, float, float] | None = None
    volume: float = 0.0

    def at(self, x: float, y: float, z: float) -> BodyPoint:
        """Locate an attachment or Morison point relative to this body's CG."""
        return BodyPoint(self, x, y, z)

    def move(self, dof: str, *, scale: float = 1.0,
             pivot: WorldPoint | BodyPoint | None = None) -> Motion:
        """Describe this body's contribution to one device coordinate."""
        return Motion(self, dof, scale, pivot)


@dataclass(frozen=True)
class _FloatingJoint:
    float_body: Body
    spar_body: Body
    location: WorldPoint
    pto_name: str
    damping: float
    stiffness: float
    equilibrium_position: float
    mooring_surge_stiffness: float
    hard_stops: LinearHardStops | None
    radiation_method: str | None
    added_mass_scheme: str


@dataclass(frozen=True)
class Coordinate:
    name: str
    motions: tuple[Motion, ...]


@dataclass(frozen=True)
class LinearPTO:
    name: str
    from_point: BodyPoint | WorldPoint
    to_point: BodyPoint | WorldPoint
    axis: tuple[float, float, float] | None = None
    damping: float = 0.0
    stiffness: float = 0.0
    equilibrium_position: float | None = None
    pretension: float | None = None
    control: DeclutchingControl | LatchingControl | None = None
    direct_drive: SimpleDirectDrive | None = None
    linear_generator: DirectLinearGenerator | None = None


@dataclass(frozen=True)
class SimpleDirectDrive:
    """Reactive PI control and simple generator/drivetrain parameters."""

    kp: float
    ki: float
    torque_constant: float
    gear_ratio: float
    drivetrain_inertia: float
    drivetrain_friction: float
    winding_resistance: float
    winding_inductance: float


@dataclass(frozen=True)
class RotationalPTO:
    """A torsional spring and damper acting on a rotation coordinate."""

    name: str
    coordinate: Coordinate
    damping: float = 0.0
    stiffness: float = 0.0
    equilibrium_angle: float = 0.0


@dataclass(frozen=True)
class RegularWave:
    height: float
    period: float
    direction: float = 0.0

    def as_case(self) -> dict:
        return {"type": "regular", "height": self.height,
                "period": self.period, "direction": self.direction}


@dataclass(frozen=True)
class RegularCICWave:
    """Regular incident waves with convolution-integral radiation."""

    height: float
    period: float
    direction: float = 0.0

    def as_case(self) -> dict:
        return {"type": "regularCIC", "height": self.height,
                "period": self.period, "direction": self.direction}


@dataclass(frozen=True)
class PMWave:
    """Pierson–Moskowitz sea for supported hydrodynamic or Morison bodies.

    ``height`` is significant wave height in metres, ``period`` is peak
    period in seconds, and ``direction`` is the incident heading in degrees.
    A saved phase CSV replays a MATLAB realization. Without one, ``seed``
    selects a reproducible Python realization rather than MATLAB's RNG stream.
    ``frequency_range`` narrows the BEM frequency interval in rad/s. A
    ``water_depth`` override is currently used only for fixed Morison bodies.
    """

    height: float
    period: float
    direction: float = 0.0
    seed: int | None = None
    phase_file: str | Path | None = None
    frequency_count: int = 500
    directions: tuple[float, ...] | None = None
    spreading: tuple[float, ...] | None = None
    frequency_range: tuple[float, float] | None = None
    water_depth: float | None = None

    def as_case(self) -> dict:
        if self.seed is not None and self.phase_file is not None:
            raise ValueError("PMWave uses either seed or phase_file")
        if (self.directions is None) != (self.spreading is None):
            raise ValueError("PMWave directions and spreading must be supplied together")
        wave = {"type": "pm", "height": self.height,
                "period": self.period,
                "directions": (list(self.directions) if self.directions is not None
                               else [self.direction]),
                "spreading": (list(self.spreading) if self.spreading is not None
                              else [1.0]),
                "frequency_count": self.frequency_count}
        if self.frequency_range is not None:
            wave["frequency_range"] = list(self.frequency_range)
        if self.water_depth is not None:
            wave["water_depth"] = self.water_depth
        if self.seed is not None:
            wave["seed"] = self.seed
        if self.phase_file is not None:
            wave["phase_file"] = str(self.phase_file)
        return wave


@dataclass(frozen=True)
class JONSWAPWave(PMWave):
    """JONSWAP sea using WEC-Sim's default or an explicit peak factor."""

    gamma: float | None = None

    def as_case(self) -> dict:
        wave = super().as_case()
        wave["type"] = "jonswap"
        if self.gamma is not None:
            wave["gamma"] = self.gamma
        return wave


@dataclass(frozen=True)
class ImportedSpectrumWave:
    """WEC-Sim three-column MAT spectrum with its saved phase realization.

    ``file`` resolves from ``WEC.run(base_dir=...)``. The MAT file contains
    frequency in Hz, density in m²/Hz, and phase in radians.
    """

    file: str | Path

    def as_case(self) -> dict:
        return {"type": "spectrumImport", "file": str(self.file)}


@dataclass(frozen=True)
class ImportedElevationWave:
    """Sampled time/elevation MAT record in WEC-Sim's ``elevationImport`` mode.

    ``file`` resolves from ``WEC.run(base_dir=...)``. The named MAT variable
    contains increasing time and elevation columns in seconds and metres.
    """

    file: str | Path
    variable: str = "etaData"
    reapply_force_ramp: bool = False

    def as_case(self) -> dict:
        return {"type": "elevationImport", "file": str(self.file),
                "variable": self.variable,
                "reapply_force_ramp": self.reapply_force_ramp}


@dataclass(frozen=True)
class NoWave:
    def as_case(self) -> dict:
        return {"type": "none"}


@dataclass(frozen=True)
class MotionHistory:
    position: np.ndarray
    velocity: np.ndarray


@dataclass(frozen=True)
class FlexibleModeHistory:
    position: np.ndarray
    velocity: np.ndarray
    acceleration: np.ndarray


@dataclass(frozen=True)
class DirectDriveHistory:
    shaft_velocity: np.ndarray
    shaft_torque: np.ndarray
    inertia_torque: np.ndarray
    friction_torque: np.ndarray
    generator_torque: np.ndarray
    current: np.ndarray
    voltage: np.ndarray
    resistance_loss: np.ndarray
    electrical_power: np.ndarray
    mechanical_power: np.ndarray


@dataclass(frozen=True)
class LinearGeneratorHistory:
    flux_d: np.ndarray
    flux_q: np.ndarray
    electrical_angle: np.ndarray
    friction_force: np.ndarray
    electrical_power: np.ndarray
    phase_current: np.ndarray
    phase_voltage: np.ndarray


@dataclass(frozen=True)
class PTOHistory:
    """PTO stroke in metres, or angle in radians for a rotational PTO."""

    stroke: np.ndarray
    velocity: np.ndarray
    force: np.ndarray
    absorbed_power: np.ndarray
    direct_drive: DirectDriveHistory | None = None
    linear_generator: LinearGeneratorHistory | None = None


@dataclass(frozen=True)
class WECResult:
    time: np.ndarray
    bodies: dict[str, MotionHistory]
    coordinates: dict[str, MotionHistory]
    ptos: dict[str, PTOHistory]
    wave_elevation: np.ndarray | None
    case: dict
    raw: CaseResponse
    flexible_modes: dict[str, FlexibleModeHistory] = field(default_factory=dict)
    body_forces: dict[str, np.ndarray] = field(default_factory=dict)


class WEC:
    """Build supported WEC motion/PTO layouts or fixed Morison bodies.

    Rotations and PTO stroke are linearized about the reference pose. The
    PTO axis stays fixed in world coordinates during a run. Hydrodynamic body
    order must match the body order in the supplied HDF5 file.
    """

    def __init__(self, name: str = "WEC", *, body_to_body: bool = False):
        if not isinstance(name, str) or not name:
            raise ValueError("WEC name must be a nonempty string")
        if not isinstance(body_to_body, bool):
            raise ValueError("body_to_body must be a boolean")
        self.name = name
        self.body_to_body = body_to_body
        self.bodies: list[Body] = []
        self.coordinates: list[Coordinate] = []
        self.ptos: list[LinearPTO] = []
        self.rotational_ptos: list[RotationalPTO] = []
        self._floating_gbm_body: Body | None = None
        self._floating_joint: _FloatingJoint | None = None
        self._morison_elements: list[tuple[Body, MorisonElement]] = []

    def body(self, name: str, hydro_file: str | Path, *,
             mass: str | float = "equilibrium",
             inertia: Sequence[float] = (0, 0, 0),
             hydro_body: int | None = None,
             mean_drift: str = "none",
             passive_yaw: bool = False,
             passive_yaw_threshold: float = 0.0,
             yaw_heading_bank: Sequence[float] | None = None,
             geometry_file: str | Path | None = None,
             nonlinear_hydro: str | None = None,
             drag_coefficient: float = 0.0,
             drag_area: float = 0.0) -> Body:
        if any(existing.name == name for existing in self.bodies):
            raise ValueError(f"body name already exists: {name}")
        if (isinstance(passive_yaw_threshold, bool)
                or not isinstance(passive_yaw_threshold, Real)
                or not np.isfinite(passive_yaw_threshold)
                or passive_yaw_threshold < 0
                or (passive_yaw_threshold > 0 and passive_yaw is not True)):
            raise ValueError("passive_yaw_threshold needs a nonnegative degree value and passive_yaw=True")
        if yaw_heading_bank is not None:
            headings = np.asarray(yaw_heading_bank, dtype=float)
            if (not passive_yaw or passive_yaw_threshold
                    or headings.ndim != 1 or len(headings) < 2
                    or not np.isfinite(headings).all()
                    or not np.all(np.diff(headings) > 0)
                    or headings[0] < -180 or headings[-1] > 180):
                raise ValueError("yaw_heading_bank needs passive_yaw, no threshold, and ordered directions in [-180, 180]")
            yaw_heading_bank = tuple(float(value) for value in headings)
        body = Body(name, hydro_file, mass, tuple(inertia), hydro_body,
                    mean_drift, passive_yaw, geometry_file, nonlinear_hydro,
                    drag_coefficient, drag_area,
                    passive_yaw_threshold=passive_yaw_threshold,
                    yaw_heading_bank=yaw_heading_bank)
        self.bodies.append(body)
        return body

    def variable_body(self, name: str, states: Sequence[HydroState], *,
                      switch_times: Sequence[float]) -> Body:
        """Add a body whose draft, mass, and BEM data switch together.

        Current dynamics support one heave-only body in regular waves.
        ``switch_times`` must contain one grid-aligned time per transition.
        """
        ordered = tuple(states)
        if (len(ordered) < 2 or not all(isinstance(s, HydroState) for s in ordered)
                or len(switch_times) != len(ordered) - 1):
            raise ValueError("variable body needs states and one fewer switch times")
        if any(existing.name == name for existing in self.bodies):
            raise ValueError(f"body name already exists: {name}")
        first = ordered[0]
        body = Body(name, first.hydro_file, first.mass, first.inertia,
                    variable_hydro=VariableHydro(ordered, tuple(switch_times)))
        self.bodies.append(body)
        return body

    def fixed_body(self, name: str, *, center_gravity: Sequence[float],
                   mass: str | float = "equilibrium",
                   inertia: Sequence[float] = (0, 0, 0),
                   volume: float = 0.0) -> Body:
        """Add a stationary body without HDF5 hydrodynamics."""
        if any(existing.name == name for existing in self.bodies):
            raise ValueError(f"body name already exists: {name}")
        center = tuple(center_gravity)
        if len(center) != 3 or not np.isfinite(center).all():
            raise ValueError("center_gravity needs three finite coordinates")
        body = Body(name, None, mass, tuple(inertia), fixed=True,
                    center_gravity=center, volume=volume)
        self.bodies.append(body)
        return body

    def morison_element(self, body: Body, *, point: BodyPoint,
                        drag_coefficient: Sequence[float],
                        added_mass_coefficient: Sequence[float],
                        area: Sequence[float], volume: float,
                        phase_mode: str = "directional") -> MorisonElement:
        """Attach a Cartesian Morison element at a body-local point.

        Moving hydrodynamic bodies support axial elements in pure heave or
        regular-wave surge/heave/pitch motion; the case runner validates
        those layouts.
        """
        if not any(body is item for item in self.bodies):
            raise ValueError("Morison body must belong to this WEC")
        if not isinstance(point, BodyPoint) or point.body is not body:
            raise ValueError("Morison point must belong to its body")
        element = MorisonElement(
            tuple(point.coordinates()), tuple(drag_coefficient),
            tuple(added_mass_coefficient), tuple(area), volume, phase_mode,
        )
        self._morison_elements.append((body, element))
        return element

    def floating_gbm(self, body: Body) -> None:
        """Select a floating surge/heave/pitch joint with HDF5 flexible modes.

        The joint and body center of gravity must coincide at the origin.
        Current support is one body in zero-heading regular waves, without a
        PTO or mooring.
        """
        if not any(body is item for item in self.bodies):
            raise ValueError("floating GBM body must belong to this WEC")
        if self._floating_gbm_body is not None:
            raise ValueError("a floating GBM body is already selected")
        self._floating_gbm_body = body

    def floating_joint(self, float_body: Body, spar_body: Body, *,
                       location: WorldPoint = WorldPoint(0, 0, 0),
                       pto_name: str = "relative_heave",
                       damping: float = 0.0, stiffness: float = 0.0,
                       equilibrium_position: float = 0.0,
                       mooring_surge_stiffness: float = 0.0,
                       hard_stops: LinearHardStops | None = None,
                       radiation_method: str | None = None,
                       added_mass_scheme: str = "implicit") -> None:
        """Select the paired two-body surge/heave/pitch slider joint.

        The PTO acts on float heave minus spar heave. The reduced joint does
        not model off-axis PTO endpoints or arbitrary Simscape constraints.
        """
        if (len(self.bodies) != 2 or self.bodies[0] is not float_body
                or self.bodies[1] is not spar_body or self._floating_joint is not None
                or self._floating_gbm_body is not None or self.coordinates
                or self.ptos or self.rotational_ptos or self._morison_elements):
            raise ValueError("floating_joint needs exactly two ordered bodies and no other layout or PTO")
        if not isinstance(location, WorldPoint) or location.x != 0 or location.y != 0:
            raise ValueError("floating_joint location needs a world point on the z axis")
        if not isinstance(pto_name, str) or not pto_name:
            raise ValueError("floating_joint PTO name must be nonempty")
        if hard_stops is not None and not isinstance(hard_stops, LinearHardStops):
            raise TypeError("hard_stops must be LinearHardStops")
        if radiation_method not in (None, "constant", "convolution", "fir"):
            raise ValueError("unsupported floating_joint radiation method")
        if added_mass_scheme not in ("implicit", "simulink_delay"):
            raise ValueError("unsupported floating_joint added-mass scheme")
        self._floating_joint = _FloatingJoint(
            float_body, spar_body, location, pto_name, damping, stiffness,
            equilibrium_position, mooring_surge_stiffness, hard_stops,
            radiation_method, added_mass_scheme,
        )

    def coordinate(self, name: str, *motions: Motion) -> Coordinate:
        if any(existing.name == name for existing in self.coordinates):
            raise ValueError(f"coordinate name already exists: {name}")
        if not motions:
            raise ValueError("a coordinate needs body motions")
        if any(not isinstance(motion, Motion)
               or not any(motion.body is body for body in self.bodies)
               for motion in motions):
            raise ValueError("coordinate motions must use bodies in this WEC")
        coordinate = Coordinate(name, tuple(motions))
        self.coordinates.append(coordinate)
        return coordinate

    def pto(self, name: str, from_point: BodyPoint | WorldPoint,
            to_point: BodyPoint | WorldPoint, *,
            axis: Sequence[float] | None = None,
            damping: float = 0.0, stiffness: float = 0.0,
            equilibrium_position: float | None = None,
            pretension: float | None = None,
            control: DeclutchingControl | LatchingControl | None = None,
            direct_drive: SimpleDirectDrive | None = None,
            linear_generator: DirectLinearGenerator | None = None) -> LinearPTO:
        if any(existing.name == name for existing in
               (*self.ptos, *self.rotational_ptos)):
            raise ValueError(f"PTO name already exists: {name}")
        for point in (from_point, to_point):
            if not isinstance(point, (BodyPoint, WorldPoint)):
                raise TypeError("PTO endpoints must be body points or world points")
            if (isinstance(point, BodyPoint)
                    and not any(point.body is body for body in self.bodies)):
                raise ValueError("PTO body point must belong to this WEC")
        if control is not None:
            if not isinstance(control, (DeclutchingControl, LatchingControl)):
                raise TypeError("PTO control must be a supported control law")
            if damping or stiffness or equilibrium_position is not None or pretension is not None:
                raise ValueError("sampled PTO control sets its own gain and has no spring")
        if direct_drive is not None:
            if not isinstance(direct_drive, SimpleDirectDrive):
                raise TypeError("direct_drive must be SimpleDirectDrive")
            if (control is not None or damping or stiffness
                    or equilibrium_position is not None or pretension is not None
                    or linear_generator is not None):
                raise ValueError("direct drive supplies its own PTO force")
        if linear_generator is not None:
            if not isinstance(linear_generator, DirectLinearGenerator):
                raise TypeError("linear_generator must be DirectLinearGenerator")
            if (control is not None or damping or stiffness
                    or equilibrium_position is not None or pretension is not None):
                raise ValueError("linear generator supplies its own PTO force")
        pto = LinearPTO(
            name, from_point, to_point,
            tuple(axis) if axis is not None else None,
            damping, stiffness, equilibrium_position, pretension, control,
            direct_drive, linear_generator,
        )
        self.ptos.append(pto)
        return pto

    def rotational_pto(self, name: str, coordinate: Coordinate, *,
                       damping: float = 0.0, stiffness: float = 0.0,
                       equilibrium_angle: float = 0.0) -> RotationalPTO:
        """Attach a torsional PTO to a named rotation coordinate.

        Damping is in N m s/rad, stiffness in N m/rad, and equilibrium angle
        in radians. The returned PTO history reports angle as ``stroke``.
        """
        if not any(coordinate is item for item in self.coordinates):
            raise ValueError("rotational PTO coordinate must belong to this WEC")
        if any(existing.name == name for existing in
               (*self.ptos, *self.rotational_ptos)):
            raise ValueError(f"PTO name already exists: {name}")
        pto = RotationalPTO(name, coordinate, damping, stiffness,
                            equilibrium_angle)
        self.rotational_ptos.append(pto)
        return pto

    def to_case(
        self, wave: RegularWave | RegularCICWave | PMWave | JONSWAPWave | ImportedSpectrumWave | ImportedElevationWave | NoWave, *,
        dt: float, end_time: float,
        ramp_time: float | None = None,
        radiation_memory: float | None = None,
        rho: float | None = None, g: float | None = None,
        initial_coordinate: Mapping[str, float] | Sequence[float] | None = None,
        initial_speed: Mapping[str, float] | Sequence[float] | None = None,
    ) -> dict:
        """Return the case mapping used by the validated dynamics runner."""
        if not isinstance(wave, (RegularWave, RegularCICWave, PMWave,
                                 JONSWAPWave, ImportedSpectrumWave,
                                 ImportedElevationWave, NoWave)):
            raise TypeError("wave must be RegularWave, RegularCICWave, PMWave, JONSWAPWave, ImportedSpectrumWave, ImportedElevationWave, or NoWave")
        simulation = {"dt": dt, "end_time": end_time}
        for key, value in (
            ("ramp_time", ramp_time), ("radiation_memory", radiation_memory),
            ("rho", rho), ("g", g),
        ):
            if value is not None:
                simulation[key] = value
        if self._floating_joint is not None:
            return self._floating_joint_case(
                wave, simulation, initial_coordinate, initial_speed,
            )
        bodies = []
        for index, body in enumerate(self.bodies, start=1):
            if body.fixed:
                body_case = {
                    "name": body.name, "nonhydro": True, "fixed": True,
                    "center_gravity": list(body.center_gravity),
                    "mass": body.mass, "inertia": list(body.inertia),
                    "volume": body.volume,
                    "morison_elements": [
                        {"point": list(element.point),
                         "drag_coefficient": list(element.drag_coefficient),
                         "added_mass_coefficient": list(element.added_mass_coefficient),
                         "area": list(element.area), "volume": element.volume,
                         "phase_mode": element.phase_mode}
                        for attached_body, element in self._morison_elements
                        if attached_body is body
                    ],
                }
                bodies.append(body_case)
                continue
            body_case = {
                "name": body.name,
                "hydro_file": str(body.hydro_file),
                "hydro_body": body.hydro_body if body.hydro_body is not None else index,
                "mass": body.mass,
                "inertia": list(body.inertia),
            }
            if body.mean_drift != "none":
                body_case["mean_drift"] = body.mean_drift
            if body.passive_yaw:
                body_case["passive_yaw"] = True
            if body.passive_yaw_threshold:
                body_case["passive_yaw_threshold"] = body.passive_yaw_threshold
            if body.yaw_heading_bank is not None:
                body_case["yaw_heading_bank"] = list(body.yaw_heading_bank)
            if body.geometry_file is not None:
                body_case["geometry_file"] = str(body.geometry_file)
            if body.nonlinear_hydro is not None:
                body_case["nonlinear_hydro"] = body.nonlinear_hydro
            if body.drag_coefficient or body.drag_area:
                body_case["drag_coefficient"] = body.drag_coefficient
                body_case["drag_area"] = body.drag_area
            if body.variable_hydro is not None:
                body_case["variable_hydro"] = {
                    "states": [
                        {"hydro_file": str(state.hydro_file),
                         "mass": state.mass,
                         "inertia": list(state.inertia)}
                        for state in body.variable_hydro.states
                    ],
                    "switch_times": list(body.variable_hydro.switch_times),
                }
            elements = [element for attached, element in self._morison_elements
                        if attached is body]
            if elements:
                body_case["morison_elements"] = [
                    {"point": list(element.point),
                     "drag_coefficient": list(element.drag_coefficient),
                     "added_mass_coefficient": list(element.added_mass_coefficient),
                     "area": list(element.area), "volume": element.volume,
                     "phase_mode": element.phase_mode}
                    for element in elements
                ]
            bodies.append(body_case)
        constraint = {"kind": "linear_subspace", "coordinates": []}
        if any(body.fixed for body in self.bodies):
            if (not all(body.fixed for body in self.bodies)
                    or self.coordinates or self.ptos or self.rotational_ptos
                    or self.body_to_body or radiation_memory is not None
                    or initial_coordinate is not None or initial_speed is not None):
                raise ValueError("fixed Morison bodies need no moving coordinates, PTOs, or radiation")
            constraint = {"kind": "fixed_morison"}
        if self._floating_gbm_body is not None:
            if (len(self.bodies) != 1 or self.bodies[0] is not self._floating_gbm_body
                    or self.coordinates or self.ptos or self.rotational_ptos
                    or self.body_to_body or not isinstance(wave, RegularWave)
                    or initial_coordinate is not None or initial_speed is not None
                    or radiation_memory is not None):
                raise ValueError("floating GBM needs one body, regular waves, and no PTO or custom coordinates")
            constraint = {"kind": "floating_gbm", "location": [0, 0, 0]}
        for coordinate in self.coordinates:
            motions = []
            for motion in coordinate.motions:
                entry = {"body": motion.body.name, "dof": motion.dof,
                         "scale": motion.scale}
                if isinstance(motion.pivot, WorldPoint):
                    entry["pivot"] = {"world": motion.pivot.coordinates()}
                elif isinstance(motion.pivot, BodyPoint):
                    if motion.pivot.body is not motion.body:
                        raise ValueError("body-local pivot must belong to its motion body")
                    entry["pivot"] = {"point": motion.pivot.coordinates()}
                elif motion.pivot is not None:
                    raise TypeError("pivot must be a world or body point")
                motions.append(entry)
            constraint["coordinates"].append({"name": coordinate.name,
                                               "motions": motions})
        if initial_coordinate is not None:
            constraint["initial_coordinate"] = _state(initial_coordinate)
        if initial_speed is not None:
            constraint["initial_speed"] = _state(initial_speed)
        case = {
            "name": self.name,
            "simulation": simulation,
            "wave": wave.as_case(),
            "bodies": bodies,
            "constraint": constraint,
        }
        if self.body_to_body:
            case["body_to_body"] = True
        if self.ptos or self.rotational_ptos:
            case["ptos"] = ([self._pto_case(pto) for pto in self.ptos]
                            + [self._rotational_pto_case(pto)
                               for pto in self.rotational_ptos])
        return case

    def _floating_joint_case(self, wave, simulation,
                             initial_coordinate, initial_speed) -> dict:
        joint = self._floating_joint
        if (len(self.bodies) != 2 or self.bodies[0] is not joint.float_body
                or self.bodies[1] is not joint.spar_body or self.coordinates
                or self.ptos or self.rotational_ptos or self._morison_elements
                or self._floating_gbm_body is not None):
            raise ValueError("floating_joint cannot combine with other bodies, coordinates, or PTOs")
        if not isinstance(wave, (RegularWave, RegularCICWave,
                                 ImportedElevationWave, NoWave)):
            raise ValueError("floating_joint supports regular, regularCIC, imported elevation, or no waves")
        if (not np.isfinite(joint.location.coordinates()).all()
                or not np.isfinite([joint.damping, joint.stiffness,
                                    joint.equilibrium_position,
                                    joint.mooring_surge_stiffness]).all()
                or joint.mooring_surge_stiffness < 0):
            raise ValueError("floating_joint location, PTO, and mooring settings must be finite")
        bodies = []
        for index, body in enumerate((joint.float_body, joint.spar_body), start=1):
            if (body.fixed or body.hydro_file is None or body.mass != "equilibrium"
                    or body.mean_drift != "none" or body.passive_yaw
                    or body.geometry_file is not None or body.nonlinear_hydro is not None
                    or body.variable_hydro is not None or body.drag_coefficient
                    or body.drag_area or body.passive_yaw_threshold
                    or body.yaw_heading_bank is not None):
                raise ValueError("floating_joint needs equilibrium-mass hydrodynamic bodies without extra force models")
            inertia = np.asarray(body.inertia, dtype=float)
            if inertia.shape != (3,) or not np.isfinite(inertia).all() or inertia[1] <= 0:
                raise ValueError("floating_joint needs a positive pitch inertia on each body")
            bodies.append({
                "hydro_file": str(body.hydro_file),
                "hydro_body": body.hydro_body if body.hydro_body is not None else index,
                "mass": "equilibrium", "pitch_inertia": float(inertia[1]),
            })
        simulation["added_mass_scheme"] = joint.added_mass_scheme
        if joint.radiation_method is not None:
            simulation["radiation_method"] = joint.radiation_method
        constraint = {"kind": "floating_joint",
                      "location": joint.location.coordinates()}
        if initial_coordinate is not None:
            constraint["initial_coordinate"] = _state(initial_coordinate)
        if initial_speed is not None:
            constraint["initial_speed"] = _state(initial_speed)
        pto = {"kind": "relative_heave", "damping": joint.damping,
               "stiffness": joint.stiffness,
               "equilibrium_position": joint.equilibrium_position}
        if joint.hard_stops is not None:
            pto["hard_stops"] = asdict(joint.hard_stops)
        case = {"name": self.name, "simulation": simulation,
                "wave": wave.as_case(), "bodies": bodies,
                "constraint": constraint, "pto": pto}
        if self.body_to_body:
            case["body_to_body"] = True
        if joint.mooring_surge_stiffness:
            case["mooring"] = {"kind": "joint_surge_spring",
                               "stiffness": joint.mooring_surge_stiffness}
        return case

    def run(
        self, wave: RegularWave | RegularCICWave | PMWave | JONSWAPWave | ImportedSpectrumWave | ImportedElevationWave | NoWave, *,
        dt: float, end_time: float,
        ramp_time: float | None = None,
        radiation_memory: float | None = None,
        rho: float | None = None, g: float | None = None,
        initial_coordinate: Mapping[str, float] | Sequence[float] | None = None,
        initial_speed: Mapping[str, float] | Sequence[float] | None = None,
        base_dir: str | Path = ".",
    ) -> WECResult:
        """Simulate the WEC and return named NumPy motion and PTO histories."""
        case = self.to_case(
            wave, dt=dt, end_time=end_time, ramp_time=ramp_time,
            radiation_memory=radiation_memory, rho=rho, g=g,
            initial_coordinate=initial_coordinate, initial_speed=initial_speed,
        )
        response = run_case(case, base_dir=base_dir)
        extras = dict(response.extra_outputs)
        bodies = {
            body.name: MotionHistory(
                response.body_position[:, index, :],
                response.body_velocity[:, index, :],
            )
            for index, body in enumerate(self.bodies)
        }
        if self._floating_joint is not None:
            coordinates = {
                name: MotionHistory(
                    response.coordinate_position[:, index],
                    response.coordinate_velocity[:, index],
                )
                for index, name in enumerate(
                    ("surge", "float_heave", "spar_heave", "pitch")
                )
            }
            pto = PTOHistory(
                response.pto_stroke,
                response.pto_velocity,
                response.pto_force,
                response.pto_absorbed_power,
            )
            return WECResult(
                response.time, bodies, coordinates,
                {self._floating_joint.pto_name: pto},
                response.wave_elevation, case, response,
            )
        coordinates = {
            coordinate.name: MotionHistory(
                extras[f"coordinate_{coordinate.name}_position"],
                extras[f"coordinate_{coordinate.name}_velocity"],
            )
            for coordinate in self.coordinates
        }
        ptos = {
            pto.name: PTOHistory(
                extras[f"pto_{pto.name}_stroke"],
                extras[f"pto_{pto.name}_velocity"],
                extras[f"pto_{pto.name}_force"],
                extras[f"pto_{pto.name}_absorbed_power"],
                (DirectDriveHistory(*(
                    extras[f"pto_{pto.name}_drive_{field}"] for field in (
                        "shaft_velocity", "shaft_torque", "inertia_torque",
                        "friction_torque", "generator_torque", "current",
                        "voltage", "resistance_loss", "electrical_power",
                        "mechanical_power",
                    )
                )) if isinstance(pto, LinearPTO) and pto.direct_drive is not None
                 else None),
                (LinearGeneratorHistory(*(
                    extras[f"pto_{pto.name}_generator_{key}"] for key in (
                        "flux_d", "flux_q", "angle", "friction_force",
                        "electrical_power", "phase_current", "phase_voltage",
                    )
                )) if isinstance(pto, LinearPTO) and pto.linear_generator is not None
                 else None),
            )
            for pto in (*self.ptos, *self.rotational_ptos)
        }
        flexible_modes = ({
            self._floating_gbm_body.name: FlexibleModeHistory(
                extras["flex_position"], extras["flex_velocity"],
                extras["flex_acceleration"],
            )
        } if self._floating_gbm_body is not None else {})
        body_forces = {
            body.name: extras[f"morison_force_{body.name}"]
            for body in self.bodies
            if body.fixed or any(attached is body
                                 for attached, _ in self._morison_elements)
        }
        return WECResult(response.time, bodies, coordinates, ptos,
                         response.wave_elevation, case, response,
                         flexible_modes, body_forces)

    @staticmethod
    def _pto_case(pto: LinearPTO) -> dict:
        def endpoint(point):
            if isinstance(point, BodyPoint):
                return {"body": point.body.name, "point": point.coordinates()}
            return {"ground": point.coordinates()}

        item = {
            "name": pto.name,
            "kind": "linear_actuator",
            "from": endpoint(pto.from_point),
            "to": endpoint(pto.to_point),
            "damping": (pto.control.gain if pto.control is not None
                        else pto.damping),
            "stiffness": pto.stiffness,
        }
        if pto.control is not None:
            if isinstance(pto.control, DeclutchingControl):
                item["control"] = {
                    "kind": "declutching",
                    "declutch_time": pto.control.declutch_time,
                    "minimum_on_time": pto.control.minimum_on_time,
                }
            else:
                item["control"] = {
                    "kind": "latching",
                    "latch_time": pto.control.latch_time,
                    "latch_damping": pto.control.latch_damping,
                    "minimum_normal_time": pto.control.minimum_normal_time,
                }
        if pto.axis is not None:
            item["axis"] = list(pto.axis)
        if pto.equilibrium_position is not None:
            item["equilibrium_position"] = pto.equilibrium_position
        if pto.pretension is not None:
            item["pretension"] = pto.pretension
        if pto.direct_drive is not None:
            item["direct_drive"] = vars(pto.direct_drive).copy()
        if pto.linear_generator is not None:
            item["linear_generator"] = vars(pto.linear_generator).copy()
        return item

    @staticmethod
    def _rotational_pto_case(pto: RotationalPTO) -> dict:
        return {
            "name": pto.name,
            "kind": "coordinate_torque",
            "coordinate": pto.coordinate.name,
            "damping": pto.damping,
            "stiffness": pto.stiffness,
            "equilibrium_position": pto.equilibrium_angle,
        }


def _state(value: Mapping[str, float] | Sequence[float]) -> dict | list:
    if isinstance(value, Mapping):
        return dict(value)
    return list(value)
