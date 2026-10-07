"""Configure a two-body WEC and PTO directly in Python.

Run from the repository root:

    python -m examples.python.configurable_rm3_pto

This is a linearized configuration example, not a MATLAB parity case.
"""

from pathlib import Path

from wecsim_python import RegularWave, WEC, WorldPoint


HYDRO = Path(__file__).resolve().parents[2] / "source/objects/rm3.h5"


def build_wec() -> WEC:
    wec = WEC("Configurable two-body WEC with offset PTO points")
    float_body = wec.body(
        "float", HYDRO, inertia=(0, 21_306_090.66, 0),
    )
    spar = wec.body(
        "spar", HYDRO, inertia=(0, 94_407_091.24, 0),
    )
    joint = WorldPoint(0, 0, 0)
    wec.coordinate("shared_surge", float_body.move("surge"), spar.move("surge"))
    wec.coordinate("float_heave", float_body.move("heave"))
    wec.coordinate("spar_heave", spar.move("heave"))
    wec.coordinate(
        "shared_pitch", float_body.move("pitch", pivot=joint),
        spar.move("pitch", pivot=joint),
    )
    wec.pto(
        "main", float_body.at(2, 0, 0.72), spar.at(-1, 0, 21.29),
        axis=(0, 0, 1), damping=1_200_000,
    )
    return wec


if __name__ == "__main__":
    result = build_wec().run(
        RegularWave(height=2.5, period=8),
        dt=0.1, end_time=20, ramp_time=5,
    )
    print(f"Simulated {len(result.time)} time steps")
    print(f"Maximum PTO force: {abs(result.ptos['main'].force).max():.0f} N")
    print(f"Maximum float pitch: {abs(result.coordinates['shared_pitch'].position).max():.4f} rad")
