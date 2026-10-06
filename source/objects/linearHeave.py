"""Linear, single-body heave free decay using WEC-Sim hydrodynamic data.

This solver covers the published Sphere ``noWaveCIC`` free-decay cases. It
uses the production BodyClass preprocessing and a fixed-step trapezoidal
integration of the Cummins radiation convolution. Other WEC-Sim degrees of
freedom, excitation, constraints, PTOs, and nonlinear forces are outside its
scope.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .bodyClass import BodyClass


@dataclass(frozen=True)
class HeaveResponse:
    time: np.ndarray
    position: np.ndarray
    displacement: np.ndarray
    velocity: np.ndarray
    force_total: np.ndarray


def _sample_count(duration, dt):
    count = round(duration / dt)
    if count < 0 or not np.isclose(count * dt, duration, rtol=0, atol=1e-10):
        raise ValueError("durations must be nonnegative multiples of dt")
    return count


def solve_heave_free_decay(
    h5_file: str | Path,
    initial_displacement: float,
    *,
    dt: float = 0.01,
    end_time: float = 40.0,
    cic_end_time: float = 15.0,
    rho: float = 1000.0,
    g: float = 9.81,
) -> HeaveResponse:
    """Return a linear heave response relative to the HDF5 equilibrium center.

    The model assumes one unconstrained heave DOF, zero incident waves, zero
    initial velocity, and no PTO, mooring, viscous drag, or body interaction.
    ``position`` is the absolute vertical center of gravity in meters;
    ``displacement`` is relative to its equilibrium value.
    """
    if dt <= 0 or cic_end_time <= 0 or rho <= 0 or g <= 0:
        raise ValueError("dt, cic_end_time, rho, and g must be positive")
    steps = _sample_count(end_time, dt)
    memory_steps = _sample_count(cic_end_time, dt)
    time = np.arange(steps + 1) * dt
    memory_time = np.arange(memory_steps + 1) * dt

    body = BodyClass(str(h5_file))
    body.bodyNumber = 1
    body.readH5file()
    if int(np.asarray(body.dof).item()) != 6:
        raise ValueError("the heave solver requires one six-DOF hydrodynamic body")
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {
        "Drag": np.zeros((6, 6)),
        "cd": np.zeros(6),
        "characteristicArea": np.zeros(6),
    }
    body.linearDamping = np.zeros((6, 6))
    body.mass = "equilibrium"
    body.hydroForcePre(
        [], [0], len(memory_time), memory_time, [], dt, rho, g,
        "noWaveCIC", np.vstack((time, np.zeros_like(time))),
        1, [], 0, 0, 0,
    )

    stiffness = float(body.hydroForce["linearHydroRestCoef"][2, 2])
    mass = float(np.asarray(body.mass).item() + body.hydroForce["fAddedMass"][2, 2])
    kernel = np.asarray(body.hydroForce["irkb"][:, 2, 2])
    equilibrium_z = float(body.hydroData["properties"]["cg"][0, 2])
    if mass <= 0:
        raise ValueError("body mass plus infinite-frequency added mass must be positive")

    displacement = np.zeros(steps + 1)
    velocity = np.zeros(steps + 1)
    acceleration = np.zeros(steps + 1)
    displacement[0] = initial_displacement
    acceleration[0] = -stiffness * initial_displacement / mass
    denominator = 1 + dt / (2 * mass) * (stiffness * dt / 2 + dt * kernel[0] / 2)

    for step in range(1, steps + 1):
        memory = min(step, memory_steps)
        known_radiation = dt * np.dot(
            kernel[1:memory + 1], velocity[step - memory:step][::-1],
        )
        velocity[step] = (
            velocity[step - 1] + dt * acceleration[step - 1] / 2
            - dt / (2 * mass) * (
                stiffness * displacement[step - 1]
                + stiffness * dt * velocity[step - 1] / 2
                + known_radiation
            )
        ) / denominator
        displacement[step] = (
            displacement[step - 1]
            + dt * (velocity[step - 1] + velocity[step]) / 2
        )
        acceleration[step] = (
            -stiffness * displacement[step]
            - known_radiation
            - dt * kernel[0] * velocity[step] / 2
        ) / mass

    return HeaveResponse(
        time=time,
        position=equilibrium_z + displacement,
        displacement=displacement,
        velocity=velocity,
        force_total=mass * acceleration,
    )
