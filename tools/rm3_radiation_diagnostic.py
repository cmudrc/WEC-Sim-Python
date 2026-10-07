"""Compare source BEM and fitted radiation damping for RM3 common surge.

The fitted response deliberately omits HDF5 feedthrough D because the pinned
MATLAB WEC-Sim body class uses zero D for the state-space radiation block.
This diagnoses the model used by Cases 5–6; it does not certify passivity.
"""

import argparse
from pathlib import Path

import h5py
import numpy as np


def common_surge_curves(hydro_file):
    """Return frequencies and effective common-surge damping in N s/m.

    Row zero omits cross-body radiation, as in Case 5. Row one includes it,
    as in Case 6. The value is the total surge force divided by a shared
    surge velocity, so both body outputs contribute.
    """
    with h5py.File(hydro_file) as h5:
        rho = float(np.asarray(h5["simulation_parameters/rho"]).item())
        if not np.isfinite(rho) or rho <= 0:
            raise ValueError("hydrodynamic density must be positive")
        source_frequency = np.asarray(h5["simulation_parameters/w"]).ravel()
        if (len(source_frequency) < 2
                or not np.all(np.diff(source_frequency) > 0)):
            raise ValueError("hydrodynamic frequencies must increase")
        fitted_frequency = np.r_[0.0, source_frequency]
        source = np.zeros((2, len(source_frequency)))
        fitted = np.zeros((2, len(fitted_frequency)))

        for output_body in range(2):
            prefix = f"body{output_body + 1}/hydro_coeffs/radiation_damping"
            source_blocks = np.asarray(h5[f"{prefix}/all"])
            fit = h5[f"{prefix}/state_space"]
            orders = np.asarray(fit["it"], dtype=int)
            A = np.asarray(fit["A/all"])
            B = np.asarray(fit["B/all"])
            C = np.asarray(fit["C/all"])
            if source_blocks.shape != (6, 12, len(source_frequency)):
                raise ValueError("expected two six-DOF RM3 hydrodynamic bodies")

            for input_body in range(2):
                column = 6 * input_body
                variants = ((0, 1) if input_body == output_body else (1,))
                for coupled in variants:
                    source[coupled] += rho * source_blocks[0, column]
                    order = orders[0, column]
                    if order < 0 or order > A.shape[-1]:
                        raise ValueError("invalid fitted radiation order")
                    if order == 0:
                        continue
                    state_A = A[0, column, :order, :order]
                    state_B = B[0, column, :order, 0]
                    state_C = rho * C[0, column, 0, :order]
                    identity = np.eye(order)
                    for index, frequency in enumerate(fitted_frequency):
                        fitted[coupled, index] += np.real(
                            state_C @ np.linalg.solve(
                                1j * frequency * identity - state_A,
                                state_B,
                            )
                        )

    if not (np.isfinite(source).all() and np.isfinite(fitted).all()):
        raise ValueError("radiation damping contains nonfinite values")
    return source_frequency, source, fitted_frequency, fitted


def plot_curves(source_frequency, source, fitted_frequency, fitted, output):
    """Save a two-panel diagnostic plot for the published RM3 Cases 5–6."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(
        1, 2, figsize=(11, 4.5), sharey=True, layout="constrained",
    )
    for coupled, axis in enumerate(axes):
        axis.axhline(0, color="#303030", linewidth=0.8)
        axis.axvline(2 * np.pi / 8, color="#777777", linestyle=":",
                     linewidth=1.5)
        axis.plot(source_frequency, source[coupled] / 1000,
                  color="#236b82", linewidth=2, label="Source BEM damping")
        axis.plot(fitted_frequency, fitted[coupled] / 1000,
                  color="#c74558", linewidth=2, label="MATLAB state-space fit")
        axis.fill_between(fitted_frequency, fitted[coupled] / 1000, 0,
                          where=fitted[coupled] < 0,
                          color="#c74558", alpha=0.12)
        axis.set_xlim(0, 1.5)
        axis.set_ylim(-30, 220)
        axis.set_xlabel("Angular frequency (rad/s)")
        axis.set_title("Cross-body radiation " + ("on" if coupled else "off"))
        axis.grid(alpha=0.2)
    axes[0].set_ylabel("Effective common-surge damping (kN s/m)")
    axes[0].legend(loc="upper left", frameon=False)
    fig.suptitle("RM3 Cases 5–6: fitted damping changes sign near zero frequency")
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hydro_file", type=Path)
    parser.add_argument("--plot", type=Path)
    args = parser.parse_args()
    source_frequency, source, fitted_frequency, fitted = common_surge_curves(
        args.hydro_file,
    )
    for coupled in (0, 1):
        print(f"cross-body {'on' if coupled else 'off'}: "
              f"fit DC {fitted[coupled, 0] / 1000:.4f} kN s/m; "
              f"source minimum {source[coupled].min() / 1000:.4g} kN s/m")
    if args.plot:
        plot_curves(source_frequency, source, fitted_frequency, fitted,
                    args.plot)
        print(args.plot)


if __name__ == "__main__":
    main()
