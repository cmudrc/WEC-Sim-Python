"""Pair the Python-configured, fully coupled RM3 MoorDyn run with MATLAB."""

import os
from pathlib import Path
import shutil

import numpy as np
import pytest

from wecsim import ImportedElevationWave, MoorDyn, WEC


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
LIBRARY = os.environ.get("WEC_SIM_MOORDYN_LIBRARY")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE and LIBRARY),
    reason="pinned RM3 MoorDyn input, output, and native library not provided",
)


def _bound(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit:.6g}"


def test_public_rm3_moordyn_coupled_400_seconds(tmp_path):
    apps = Path(APPLICATIONS)
    reference = Path(REFERENCE)
    source_lines = apps / "Mooring/MoorDyn/Mooring/lines.txt"
    input_dir = tmp_path / "Mooring"
    input_dir.mkdir()
    lines = input_dir / "lines.txt"
    shutil.copyfile(source_lines, lines)

    wec = WEC("RM3 MoorDyn")
    float_body = wec.body(
        "float", "_Common_Input_Files/RM3/hydroData/rm3.h5",
        inertia=(0, 21_306_090.66, 0),
    )
    spar = wec.body(
        "spar", "_Common_Input_Files/RM3/hydroData/rm3.h5",
        inertia=(0, 94_407_091.24, 0),
    )
    wec.floating_joint(
        float_body, spar, damping=1_200_000,
        moordyn=MoorDyn(LIBRARY, lines),
        moordyn_point=spar.at(0, 0, 21.5),
    )
    result = wec.run(
        ImportedElevationWave(
            "Mooring/MoorDyn/etaData.mat", reapply_force_ramp=True,
        ),
        dt=0.01, end_time=400, ramp_time=40, radiation_memory=60,
        initial_coordinate={"spar_heave": -0.21}, base_dir=apps,
    )
    assert result.time.shape == (40_001,)
    wave = np.loadtxt(reference / "RM3_MOORDYN_wave.csv", delimiter=",")
    _bound(result.wave_elevation, wave[:, 1], 1e-10, "wave elevation")

    for number, name in ((1, "float"), (2, "spar")):
        saved = np.loadtxt(reference / f"RM3_MOORDYN_MoorDyn_body{number}.csv",
                           delimiter=",")
        _bound(result.time, saved[:, 0], 1e-8, "time")
        for axis, position_limit, speed_limit in (
            (0, 0.003, 0.001),
            (2, 0.0005, 0.0005),
            (4, 0.00005, 0.00005),
        ):
            _bound(result.bodies[name].position[:, axis],
                   saved[:, 1 + axis], position_limit,
                   f"{name} position axis {axis}")
            _bound(result.bodies[name].velocity[:, axis],
                   saved[:, 7 + axis], speed_limit,
                   f"{name} velocity axis {axis}")

    pto = np.loadtxt(reference / "RM3_MOORDYN_MoorDyn_pto1.csv",
                     delimiter=",")
    _bound(result.ptos["relative_heave"].force, pto[:, 15],
           400, "PTO force")
    outputs = dict(result.raw.extra_outputs)
    coupling = np.loadtxt(reference / "RM3_MOORDYN_coupling.csv",
                          delimiter=",")
    for kind, column, limit in (
        ("position", 1, (0.003, 1e-6, 0.0005, 1e-6, 0.00005, 1e-6)),
        ("velocity", 7, (0.001, 1e-6, 0.0005, 1e-6, 0.00005, 1e-6)),
        ("force", 13, (500, 1e-3, 2_500, 1e-3, 1_200, 1e-3)),
    ):
        actual = outputs[f"moordyn_connection_{kind}"]
        for axis, bound in enumerate(limit):
            _bound(actual[:, axis], coupling[:, column + axis], bound,
                   f"connection {kind} axis {axis}")
    tensions = np.loadtxt(input_dir / "lines.out", skiprows=1)
    saved_tensions = np.loadtxt(
        reference / "RM3_MOORDYN_fairlead_tension.csv", delimiter=",",
    )
    _bound(tensions[:, 1:], saved_tensions[:, 1:], 2_000,
           "three fairlead tensions")
