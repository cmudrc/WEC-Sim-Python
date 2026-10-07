"""Linearized PTO connections between body-local points and fixed anchors.

Each connection measures stroke along a fixed world-space axis from its
``from`` endpoint to its ``to`` endpoint. At the reference pose stroke is
zero. Body attachment points are measured from the body's HDF5 center of
gravity. Small body rotations contribute through angular displacement cross
the attachment offset, which is consistent with ``linear_subspace`` maps.
"""

from dataclasses import dataclass
import re
from typing import Mapping

import numpy as np

from .controls import DeclutchingControl


@dataclass(frozen=True)
class LinearPTOConnection:
    name: str
    stroke_jacobian: np.ndarray
    damping: float
    stiffness: float
    equilibrium_position: float
    control: DeclutchingControl | None = None

    def force(self, coordinate: np.ndarray, speed: np.ndarray) -> np.ndarray:
        if self.control is not None:
            raise ValueError("controlled PTO force requires controller memory")
        stroke = coordinate @ self.stroke_jacobian
        stroke_speed = speed @ self.stroke_jacobian
        return (-self.stiffness * (stroke - self.equilibrium_position)
                - self.damping * stroke_speed)

    def outputs(self, coordinate: np.ndarray,
                speed: np.ndarray, *,
                force_override: np.ndarray | None = None) -> tuple[tuple[str, np.ndarray], ...]:
        stroke = coordinate @ self.stroke_jacobian
        stroke_speed = speed @ self.stroke_jacobian
        prefix = f"pto_{self.name}"
        if self.control is not None and force_override is None:
            raise ValueError("controlled PTO output needs the simulated force")
        force = (self.force(coordinate, speed) if force_override is None
                 else np.asarray(force_override, dtype=float))
        if force.shape != stroke_speed.shape:
            raise ValueError("controlled PTO force must match the time grid")
        absorbed_power = (self.damping * stroke_speed**2
                          if force_override is None else -force * stroke_speed)
        return (
            (f"{prefix}_stroke", stroke),
            (f"{prefix}_velocity", stroke_speed),
            (f"{prefix}_force", force),
            (f"{prefix}_absorbed_power", absorbed_power),
        )


def _object(value, name, required, allowed):
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be an object")
    missing = required - value.keys()
    extra = value.keys() - allowed
    if missing or extra:
        raise ValueError(f"{name}: missing {sorted(missing)}, unsupported {sorted(extra)}")
    return value


def _vector(value, name):
    vector = np.asarray(value, dtype=float)
    if vector.shape != (3,) or not np.isfinite(vector).all():
        raise ValueError(f"{name} must be a finite three-vector")
    return vector


def _number(value, name, *, nonnegative=False):
    if isinstance(value, bool):
        raise ValueError(f"{name} must be numeric")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not np.isfinite(number) or (nonnegative and number < 0):
        raise ValueError(f"{name} is outside its supported range")
    return number


def _endpoint(spec, name, maps, centers, names):
    _object(spec, name, set(), {"body", "point", "ground"})
    if "ground" in spec:
        if set(spec) != {"ground"}:
            raise ValueError(f"{name} ground endpoint only uses ground")
        point = _vector(spec["ground"], f"{name}.ground")
        return point, np.zeros((3, maps[0].shape[1]))
    if set(spec) != {"body", "point"}:
        raise ValueError(f"{name} body endpoint needs body and point")
    body = spec["body"]
    if isinstance(body, bool):
        raise ValueError(f"{name}.body must identify a body")
    if isinstance(body, int):
        index = body - 1
    elif isinstance(body, str) and body in names:
        index = names.index(body)
    else:
        raise ValueError(f"{name}.body must identify a body")
    if index < 0 or index >= len(maps):
        raise ValueError(f"{name}.body must identify a body")
    offset = _vector(spec["point"], f"{name}.point")
    mapping = maps[index]
    point_jacobian = mapping[:3] + np.cross(mapping[3:].T, offset).T
    return centers[index] + offset, point_jacobian


def build_linear_ptos(specs, maps, centers, names):
    """Return connections and their total linear generalized force law.

    ``K``, ``C``, and ``bias`` satisfy ``F = -K q - C q_dot + bias``.
    """
    if not isinstance(specs, list) or not specs:
        raise ValueError("ptos must be a nonempty list")
    if len(maps) != len(centers) or len(maps) != len(names):
        raise ValueError("PTO body maps, centers, and names must align")
    n = maps[0].shape[1]
    stiffness = np.zeros((n, n))
    damping = np.zeros((n, n))
    bias = np.zeros(n)
    connections = []
    seen = set()
    for position, raw in enumerate(specs, start=1):
        spec = _object(raw, f"ptos[{position}]", {"name", "kind", "from", "to"},
                       {"name", "kind", "from", "to", "axis", "damping",
                        "stiffness", "equilibrium_position", "pretension",
                        "control"})
        name = spec["name"]
        if (not isinstance(name, str)
                or re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name) is None
                or re.fullmatch(r"coordinate[1-9][0-9]*", name) is not None
                or name in seen):
            raise ValueError("PTO names must be unique identifiers")
        seen.add(name)
        if spec["kind"] != "linear_actuator":
            raise ValueError("ptos currently support linear_actuator only")
        first, first_j = _endpoint(spec["from"], f"ptos[{position}].from",
                                   maps, centers, names)
        second, second_j = _endpoint(spec["to"], f"ptos[{position}].to",
                                     maps, centers, names)
        direction = (_vector(spec["axis"], f"ptos[{position}].axis")
                     if "axis" in spec else second - first)
        norm = np.linalg.norm(direction)
        if norm == 0 or not np.isfinite(norm):
            raise ValueError("PTO axis needs noncoincident endpoints or an explicit axis")
        direction = direction / norm
        jacobian = direction @ (second_j - first_j)
        if not np.any(jacobian):
            raise ValueError("PTO attachment has no motion in its force axis")
        coefficient = _number(spec.get("damping", 0),
                              f"ptos[{position}].damping", nonnegative=True)
        spring = _number(spec.get("stiffness", 0),
                         f"ptos[{position}].stiffness", nonnegative=True)
        if not coefficient and not spring:
            raise ValueError("PTO actuator needs damping or stiffness")
        if "equilibrium_position" in spec and "pretension" in spec:
            raise ValueError("specify either PTO equilibrium_position or pretension")
        if "pretension" in spec:
            pretension = _number(spec["pretension"], f"ptos[{position}].pretension")
            if pretension and not spring:
                raise ValueError("nonzero PTO pretension needs positive stiffness")
            equilibrium = -pretension / spring if spring else 0.0
        else:
            equilibrium = _number(spec.get("equilibrium_position", 0),
                                  f"ptos[{position}].equilibrium_position")
        if not np.isfinite(equilibrium) or (equilibrium and not spring):
            raise ValueError("nonzero PTO equilibrium needs positive stiffness")
        control = None
        if "control" in spec:
            control_spec = _object(
                spec["control"], f"ptos[{position}].control",
                {"kind", "declutch_time"},
                {"kind", "declutch_time", "minimum_on_time"},
            )
            if control_spec["kind"] != "declutching":
                raise ValueError("only declutching PTO control is supported")
            if not coefficient or spring or equilibrium:
                raise ValueError("declutching PTO needs damping and no spring")
            control = DeclutchingControl(
                coefficient,
                _number(control_spec["declutch_time"],
                        f"ptos[{position}].control.declutch_time"),
                _number(control_spec.get("minimum_on_time", 0.2),
                        f"ptos[{position}].control.minimum_on_time"),
            )
        connection = LinearPTOConnection(name, jacobian, coefficient,
                                         spring, equilibrium, control)
        connections.append(connection)
        stiffness += spring * np.outer(jacobian, jacobian)
        if control is None:
            damping += coefficient * np.outer(jacobian, jacobian)
        bias += spring * equilibrium * jacobian
    return tuple(connections), stiffness, damping, bias
