"""Pair the published RM3 MooringMatrix motion and PTO with MATLAB."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import run_case


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="paired MATLAB mooring output and Applications checkout not provided",
)


def _max_error(actual, expected, limit, label):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = np.max(np.abs(actual - expected))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit:.6g}"


def test_rm3_mooring_matrix_full_duration_against_matlab():
    apps = Path(APPLICATIONS)
    reference = Path(REFERENCE)
    prefix = "RM3_MOORING_MATRIX"
    hydro = apps / "_Common_Input_Files/RM3/hydroData/rm3.h5"
    dt, end_time, ramp_time = 0.01, 400.0, 40.0
    result = run_case(
        {
            "simulation": {"dt": dt, "end_time": end_time,
                           "ramp_time": ramp_time, "radiation_memory": 60},
            "wave": {"type": "elevationImport",
                     "file": "Mooring/MooringMatrix/etaData.mat",
                     "variable": "etaData", "direction": 0,
                     "reapply_force_ramp": True},
            "bodies": [
                {"hydro_file": str(hydro), "hydro_body": number,
                 "mass": "equilibrium", "pitch_inertia": inertia}
                for number, inertia in ((1, 21_306_090.66), (2, 94_407_091.24))
            ],
            "constraint": {"kind": "floating_joint", "location": [0, 0, 0],
                           "initial_coordinate": {"spar_heave": -0.21}},
            "pto": {"kind": "relative_heave", "damping": 1_200_000},
            "mooring": {"kind": "joint_surge_spring", "stiffness": 100_000},
        },
        base_dir=apps,
    )
    assert result.auxiliary_files == (apps / "Mooring/MooringMatrix/etaData.mat",)
    assert result.wave_elevation.shape == result.time.shape
    outputs = dict(result.extra_outputs)
    source_wave = np.loadtxt(reference / f"{prefix}_wave.csv", delimiter=",")
    _max_error(result.wave_elevation, source_wave[:, 1], 1e-10,
               "imported wave elevation")
    bodies = [np.loadtxt(reference / f"{prefix}_MooringMatrix_body{number}.csv",
                         delimiter=",") for number in (1, 2)]
    for index, saved in enumerate(bodies):
        np.testing.assert_allclose(result.time, saved[:, 0], rtol=0, atol=1e-8)
        for dof, name, position_limit, velocity_limit in (
            (0, "surge", 0.015, 0.006),
            (2, "heave", 0.004, 0.0006),
            (4, "pitch", 0.0009, 0.0004),
        ):
            _max_error(result.body_position[:, index, dof], saved[:, 1 + dof],
                       position_limit, f"body{index+1} {name} position")
            _max_error(result.body_velocity[:, index, dof], saved[:, 7 + dof],
                       velocity_limit, f"body{index+1} {name} velocity")

    mooring = np.loadtxt(reference / f"{prefix}_mooring1.csv", delimiter=",")
    _max_error(outputs["mooring_surge_position"], mooring[:, 1], 0.015,
               "mooring surge position")
    _max_error(outputs["mooring_surge_force"], mooring[:, 13], 1_500,
               "mooring surge force")

    pto = np.loadtxt(reference / f"{prefix}_MooringMatrix_pto1.csv",
                     delimiter=",")
    # Source PTO translation is relative to its initialized frame. Recover
    # that datum from the reported body centers and common pitch.
    initial_gap = bodies[0][0, 3] - bodies[1][0, 3]
    pitch = result.body_position[:, 0, 4]
    pitch_speed = result.body_velocity[:, 0, 4]
    stroke = ((result.body_position[:, 0, 2]
               - result.body_position[:, 1, 2]) / np.cos(pitch) - initial_gap)
    speed = ((result.body_velocity[:, 0, 2]
              - result.body_velocity[:, 1, 2]
              + (initial_gap + stroke) * np.sin(pitch) * pitch_speed)
             / np.cos(pitch))
    _max_error(stroke, pto[:, 3], 0.004, "PTO stroke")
    _max_error(speed, pto[:, 9], 0.0005, "PTO speed")
    _max_error(result.pto_force, pto[:, 15], 500, "PTO force")
    _max_error(result.pto_force * speed, pto[:, 21], 500, "PTO power")
