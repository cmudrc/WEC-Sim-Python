"""Source-convention radiation diagnostics for pinned WEC-Sim HDF5 data."""

from pathlib import Path

import h5py
import numpy as np
from scipy.signal import lsim


def replay_rm3_fitted_radiation(
    hydro_file: str | Path, time: np.ndarray, body_velocity: np.ndarray, *,
    body_to_body: bool = False,
) -> np.ndarray:
    """Replay RM3's MATLAB state-space radiation on prescribed body motion.

    Return force with shape ``(time, body, six DOFs)``. The pinned MATLAB
    block multiplies HDF5 C by water density and sets D to zero, even though
    HDF5 contains a direct term. Each fitted channel starts with zero state.
    Linear interpolation between saved velocities approximates the source's
    internal solver input; this function does not advance a WEC trajectory.
    The published fit has negative damping in active joint coordinates, so
    this replay is separate from the default convolution dynamics.
    """
    t = np.asarray(time, dtype=float)
    velocity = np.asarray(body_velocity, dtype=float)
    if (t.ndim != 1 or len(t) < 2 or t[0] != 0
            or not np.isfinite(t).all() or not np.all(np.diff(t) > 0)
            or not np.allclose(np.diff(t), t[1] - t[0], rtol=1e-8, atol=1e-12)
            or velocity.shape != (len(t), 2, 6)
            or not np.isfinite(velocity).all()):
        raise ValueError("radiation replay needs uniform time from zero and two finite six-DOF velocity histories")
    if not isinstance(body_to_body, bool):
        raise TypeError("body_to_body must be a boolean")

    force = np.zeros_like(velocity)
    with h5py.File(hydro_file) as h5:
        rho = float(np.asarray(h5["simulation_parameters/rho"]).item())
        if not np.isfinite(rho) or rho <= 0:
            raise ValueError("hydrodynamic density must be positive")
        for output_body in range(2):
            fit = h5[f"body{output_body + 1}/hydro_coeffs/radiation_damping/state_space"]
            orders = np.asarray(fit["it"])
            A = np.asarray(fit["A/all"])
            B = np.asarray(fit["B/all"])
            C = np.asarray(fit["C/all"])
            if (orders.shape != (6, 12) or A.shape[:2] != (6, 12)
                    or B.shape[:2] != (6, 12) or C.shape[:2] != (6, 12)
                    or not np.all(orders == np.floor(orders))):
                raise ValueError("expected a two-body six-DOF radiation fit")
            for output_dof in range(6):
                for input_body in range(2):
                    if input_body != output_body and not body_to_body:
                        continue
                    for input_dof in range(6):
                        column = 6 * input_body + input_dof
                        order = int(orders[output_dof, column])
                        if order < 0 or order > A.shape[-1]:
                            raise ValueError("invalid fitted radiation order")
                        if order == 0:
                            continue
                        system = (
                            A[output_dof, column, :order, :order],
                            B[output_dof, column, :order, :1],
                            C[output_dof, column, :1, :order],
                            np.zeros((1, 1)),
                        )
                        _, response, _ = lsim(
                            system, U=velocity[:, input_body, input_dof], T=t,
                        )
                        force[:, output_body, output_dof] += rho * response
    return force
