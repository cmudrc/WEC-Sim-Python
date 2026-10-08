"""Pair a sampled-elevation Sphere heave WEC with pinned MATLAB WEC-Sim."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import ImportedElevationWave, WEC


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="derived MATLAB Sphere imported-elevation output not provided",
)


def _max_error(actual, expected, limit, label):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: max error {error:.6g} exceeds {limit}"


def test_public_sphere_imported_elevation_motion_and_force():
    root = Path(APPLICATIONS)
    reference = Path(REFERENCE)
    source = np.loadtxt(
        reference / "SPHERE_ELEVATION_IMPORT_0m_body1.csv", delimiter=",",
    )
    source_wave = np.loadtxt(
        reference / "SPHERE_ELEVATION_IMPORT_wave.csv", delimiter=",",
    )
    wec = WEC("Sphere with imported elevation")
    sphere = wec.body(
        "sphere", "_Common_Input_Files/Sphere/hydroData/sphere.h5",
        inertia=(20_907_301, 21_306_090.66, 37_085_481.11),
    )
    wec.coordinate("heave", sphere.move("heave"))
    result = wec.run(
        ImportedElevationWave(
            "Free_Decay/0m/etaData.mat", reapply_force_ramp=True,
        ),
        dt=0.01, end_time=40, ramp_time=10, radiation_memory=15,
        base_dir=root,
    )
    assert source.shape == (4001, 25)
    _max_error(result.time, source[:, 0], 1e-10, "time")
    _max_error(result.wave_elevation, source_wave[:, 1], 1e-10,
               "wave elevation")
    force = dict(result.raw.extra_outputs)["body1_excitation_force"]
    _max_error(force, source[:, 19:25], 1e-3, "six excitation forces")
    _max_error(result.bodies["sphere"].position[:, 2], source[:, 3], 0.01,
               "heave position")
    _max_error(result.bodies["sphere"].velocity[:, 2], source[:, 9], 0.01,
               "heave velocity")
    assert result.raw.auxiliary_files == (
        (root / "Free_Decay/0m/etaData.mat").resolve(),
    )
