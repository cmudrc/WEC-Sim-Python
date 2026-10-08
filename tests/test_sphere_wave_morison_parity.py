"""Pair a wave-driven axial Morison element on the published Sphere model."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularWave, WEC
from wecsim.morison import MorisonElement, regular_wave_heave_morison_terms


SPHERE_H5 = os.environ.get("WEC_SIM_SPHERE_H5")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (SPHERE_H5 and REFERENCE),
    reason="paired wave-driven Morison Sphere output and HDF5 absent",
)


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    error = float(np.max(np.abs(actual - expected)))
    assert np.isfinite(error) and error < limit, (
        f"{label}: {error:.6g} exceeds {limit}"
    )


def test_wave_driven_sphere_morison_force_and_motion():
    folder = Path(REFERENCE)
    source = np.loadtxt(
        folder / "SPHERE_MOVING_MORISON_WAVE_1m-ME_body1.csv",
        delimiter=",",
    )
    wave = np.loadtxt(folder / "SPHERE_MOVING_MORISON_WAVE_wave.csv",
                      delimiter=",")
    assert source.shape == (4001, 37)
    assert wave.shape == (4001, 2)
    assert np.isfinite(source).all() and np.isfinite(wave).all()
    element = MorisonElement(
        point=(0, 0, -2), drag_coefficient=(0, 0, 1),
        added_mass_coefficient=(0, 0, 1), area=(0, 0, 100), volume=20,
    )
    source_terms = np.array([
        regular_wave_heave_morison_terms(
            [element], center_z=-2, heave=row[3] + 2,
            speed=row[9], time=row[0], wave_height=1, wave_period=8,
            ramp_time=10, water_depth=np.inf, rho=1000,
        )
        for row in source
    ])
    source_physical = -source[:, 27]
    reconstructed = source_terms[:, 0] - source_terms[:, 1] * source[:, 33]
    _max_error(reconstructed[50:], source_physical[50:], 100,
               "wave Morison force on MATLAB's saved state after startup")

    wec = WEC("Sphere with wave-driven axial Morison element")
    sphere = wec.body("sphere", Path(SPHERE_H5).resolve())
    wec.coordinate("heave", sphere.move("heave"))
    wec.morison_element(
        sphere, point=sphere.at(0, 0, -2),
        drag_coefficient=(0, 0, 1),
        added_mass_coefficient=(0, 0, 1),
        area=(0, 0, 100), volume=20,
    )
    result = wec.run(
        RegularWave(1, 8), dt=.01, end_time=40, ramp_time=10,
        initial_coordinate={"heave": 1},
    )
    _max_error(result.time, source[:, 0], 1e-12, "time")
    _max_error(result.wave_elevation, wave[:, 1], 1e-10,
               "regular-wave elevation")
    _max_error(result.bodies["sphere"].position[:, 2], source[:, 3],
               .01, "sphere heave")
    _max_error(result.bodies["sphere"].velocity[:, 2], source[:, 9],
               .02, "sphere speed")
    physical = result.body_forces["sphere"][:, 2]
    _max_error(physical[100:], source_physical[100:], 1000,
               "wave-driven moving Morison force after startup")
