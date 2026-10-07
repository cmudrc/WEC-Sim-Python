"""Reduced four-coordinate model of the published RM3 regular-wave case.

The two bodies share surge and pitch at the floating joint and have separate
heave coordinates. Hydrodynamic coefficients come from production ``BodyClass``
preprocessing. This is a reference-case model, not the general WEC-Sim runner.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .bodyClass import BodyClass
from .generalDynamics import BodyMotion, DynamicBody, GeneralizedDynamics


@dataclass(frozen=True)
class RM3RegularResponse:
    time: np.ndarray
    body_position: np.ndarray
    body_velocity: np.ndarray
    pto_force: np.ndarray


def solve_rm3_regular(
    h5_file: str | Path,
    *,
    wave_height: float = 2.5,
    wave_period: float = 8.0,
    pitch_inertias: tuple[float, float] = (21_306_090.66, 94_407_091.24),
    pto_damping: float = 1_200_000.0,
    pto_stiffness: float = 0.0,
    pto_equilibrium: float = 0.0,
    b2b: bool = False,
    radiation_memory: float | None = None,
    radiation_method: str | None = None,
    joint_z: float = 0.0,
    dt: float = 0.1,
    end_time: float = 400.0,
    ramp_time: float = 100.0,
    rho: float = 1000.0,
    g: float = 9.81,
) -> RM3RegularResponse:
    """Solve RM3 surge, two heaves, and shared pitch with a relative heave PTO.

    The bodies use equilibrium displaced-volume masses, fixed-frequency
    added mass and radiation damping, hydrostatic restoring, and regular-wave
    excitation. Rigid-body rotation changes their surge/heave Jacobians at
    each step. Without ``radiation_memory``, fixed-frequency damping uses
    classical RK4. With radiation memory, infinite-frequency added mass and
    the radiation impulse-response kernel use a trapezoidal history step.
    ``radiation_method="fir"`` instead samples the same kernel as a discrete
    FIR filter and holds its force through each RK4 step.
    ``b2b=True`` includes cross-body radiation blocks. Sway, roll, yaw, and
    full Simscape joint forces are outside this reduced model.
    """
    inputs = [wave_height, wave_period, pto_damping, pto_stiffness,
              pto_equilibrium,
              joint_z, dt, end_time, ramp_time, rho, g, *pitch_inertias]
    if radiation_memory is not None:
        inputs.append(radiation_memory)
    if not np.isfinite(inputs).all():
        raise ValueError("solver inputs must be finite")
    if (len(pitch_inertias) != 2 or wave_height < 0 or wave_period <= 0
            or any(inertia <= 0 for inertia in pitch_inertias)
            or pto_damping < 0 or pto_stiffness < 0 or dt <= 0
            or end_time < 0 or ramp_time < 0 or rho <= 0 or g <= 0):
        raise ValueError("invalid RM3 wave, body, PTO, or time parameters")
    if radiation_memory is not None and radiation_memory <= 0:
        raise ValueError("radiation_memory must be positive")
    if radiation_method is None:
        radiation_method = ("convolution" if radiation_memory is not None
                            else "constant")
    if radiation_method not in ("constant", "convolution", "fir"):
        raise ValueError("unsupported RM3 radiation method")
    if ((radiation_method == "constant" and radiation_memory is not None)
            or (radiation_method != "constant" and radiation_memory is None)):
        raise ValueError("constant radiation has no memory; convolution and FIR need it")
    if not isinstance(b2b, bool):
        raise ValueError("b2b must be a boolean")
    steps = round(end_time / dt)
    if not np.isclose(steps * dt, end_time, rtol=0, atol=1e-10):
        raise ValueError("end_time must be an integer multiple of dt")
    time = np.arange(steps + 1) * dt
    omega = 2 * np.pi / wave_period
    if radiation_memory is not None:
        memory_steps = round(radiation_memory / dt)
        if not np.isclose(memory_steps * dt, radiation_memory, rtol=0, atol=1e-10):
            raise ValueError("radiation_memory must be an integer multiple of dt")
        convolution_time = np.arange(memory_steps + 1) * dt
    data = []
    for index, pitch_inertia in enumerate(pitch_inertias, start=1):
        body = BodyClass(str(h5_file))
        body.bodyNumber = index
        body.bodyTotal = 2
        body.readH5file()
        if int(np.asarray(body.dof).item()) != 6:
            raise ValueError("each RM3 hydrodynamic body must have six DOFs")
        body.mass = "equilibrium"
        body.hydroStiffness = np.zeros((6, 6))
        body.viscDrag = {
            "Drag": np.zeros((6, 6)),
            "cd": np.zeros(6),
            "characteristicArea": np.zeros(6),
        }
        body.linearDamping = np.zeros((6, 6))
        if radiation_memory is None:
            body.hydroForcePre(
                omega, [0], 1, np.array([0.0]), [], dt, rho, g,
                "regular", np.vstack((time, np.zeros_like(time))),
                index, 2, 0, 0, int(b2b),
            )
        else:
            irf_time = body.hydroData["hydro_coeffs"]["radiation_damping"][
                "impulse_response_fun"]["t"]
            if radiation_memory > np.max(irf_time) + 1e-10:
                raise ValueError("radiation_memory exceeds the HDF5 kernel")
            body.hydroForcePre(
                omega, [0], len(convolution_time), convolution_time, [],
                dt, rho, g, "regularCIC",
                np.vstack((time, np.zeros_like(time))),
                index, 2, 0, 0, int(b2b),
            )
        mass = float(np.asarray(body.mass).item())
        center = np.asarray(body.cg).ravel()
        if not np.isclose(center[0:2], 0.0, atol=1e-10).all():
            raise ValueError("the RM3 model requires body centers on the joint axis")
        lever = float(center[2] - joint_z)
        rigid_mass = np.diag([mass, mass, mass, 0.0, pitch_inertia, 0.0])
        hydro = body.hydroForce
        if b2b:
            added_mass = tuple(np.asarray(hydro["fAddedMass"])[:, 6*j:6*(j+1)] for j in range(2))
            damping = (tuple(np.asarray(hydro["fDamping"])[:, 6*j:6*(j+1)] for j in range(2))
                       if radiation_memory is None else tuple(np.zeros((6, 6)) for _ in range(2)))
        else:
            added_mass = [np.zeros((6, 6)), np.zeros((6, 6))]
            damping = [np.zeros((6, 6)), np.zeros((6, 6))]
            added_mass[index - 1] = np.asarray(hydro["fAddedMass"])
            if radiation_memory is None:
                damping[index - 1] = np.asarray(hydro["fDamping"])
        if radiation_memory is None:
            kernel = None
        else:
            raw_kernel = np.asarray(hydro["irkb"])
            if b2b:
                kernel = raw_kernel
            else:
                kernel = np.zeros((len(raw_kernel), 6, 12))
                kernel[:, :, 6 * (index - 1):6 * index] = raw_kernel
        data.append({
            "center_z": float(center[2]),
            "lever": lever,
            "rigid_mass": rigid_mass,
            "added_mass": added_mass,
            "damping": damping,
            "kernel": kernel,
            "restoring": np.asarray(hydro["linearHydroRestCoef"]),
            "re": np.asarray(hydro["fExt"]["re"]),
            "im": np.asarray(hydro["fExt"]["im"]),
            "vertical_bias": (rho * float(np.asarray(body.dispVol).item()) - mass) * g,
        })

    pto_coupling = np.zeros((4, 4))
    pto_coupling[1, 1] = pto_coupling[2, 2] = 1
    pto_coupling[1, 2] = pto_coupling[2, 1] = -1

    dynamic_bodies = []
    for index, body in enumerate(data):
        lever = body["lever"]

        def motion(q, v, *, index=index, lever=lever):
            angle = q[3]
            jacobian = np.zeros((6, 4))
            jacobian[0, 0] = 1.0
            jacobian[2, index + 1] = 1.0
            jacobian[0, 3] = lever * np.cos(angle)
            jacobian[2, 3] = -lever * np.sin(angle)
            jacobian[4, 3] = 1.0
            curvature = np.array([
                -lever * np.sin(angle), 0.0, -lever * np.cos(angle),
                0.0, 0.0, 0.0,
            ]) * v[3]**2
            displacement = np.array([
                q[0] + lever * np.sin(angle), 0.0,
                q[index + 1] + lever * (np.cos(angle) - 1.0),
                0.0, angle, 0.0,
            ])
            return BodyMotion(displacement, jacobian, curvature)

        def excitation(at_time, *, re=body["re"], im=body["im"]):
            ramp = (1.0 if ramp_time == 0 or at_time >= ramp_time
                    else (1.0 - np.cos(np.pi * at_time / ramp_time)) / 2)
            return wave_height / 2 * ramp * (
                re * np.cos(omega * at_time) - im * np.sin(omega * at_time)
            )

        dynamic_bodies.append(DynamicBody(
            rigid_mass=body["rigid_mass"],
            added_mass=tuple(body["added_mass"]),
            damping=tuple(body["damping"]),
            restoring=body["restoring"],
            static_force=np.array([0, 0, body["vertical_bias"], 0, 0, 0]),
            reference_position=np.array([0, 0, body["center_z"], 0, 0, 0]),
            motion=motion,
            excitation=excitation,
            radiation_kernel=body["kernel"],
        ))

    system = GeneralizedDynamics(
        tuple(dynamic_bodies), 4,
        pto_stiffness=pto_stiffness * pto_coupling,
        pto_damping=pto_damping * pto_coupling,
        pto_equilibrium=np.array([0, pto_equilibrium / 2,
                                  -pto_equilibrium / 2, 0]),
        radiation_discretization=("fir" if radiation_method == "fir"
                                  else "trapezoid"),
    )
    solved = system.integrate(dt=dt, end_time=end_time)
    q = solved.coordinate
    v = solved.speed
    return RM3RegularResponse(
        time=solved.time, body_position=solved.body_position,
        body_velocity=solved.body_velocity,
        pto_force=(-pto_damping * (v[:, 1] - v[:, 2])
                   - pto_stiffness * (q[:, 1] - q[:, 2] - pto_equilibrium)),
    )
