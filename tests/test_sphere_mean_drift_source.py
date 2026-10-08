"""Check the pinned MATLAB Sphere Mean_Drift export before dynamics pairing."""

import os
from pathlib import Path

import h5py
import numpy as np
import pytest


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="paired MATLAB mean-drift output and Applications checkout not provided",
)


def test_published_sphere_mean_drift_has_a_complete_source_response():
    root = Path(REFERENCE)
    body = np.loadtxt(root / "SPHERE_MEAN_DRIFT_Mean_Drift_body1.csv",
                      delimiter=",")
    wave = np.loadtxt(root / "SPHERE_MEAN_DRIFT_wave.csv", delimiter=",")
    coefficients = np.loadtxt(
        root / "SPHERE_MEAN_DRIFT_excitation_coefficients.csv", delimiter=",")
    mass = np.loadtxt(root / "SPHERE_MEAN_DRIFT_mass_properties.csv",
                      delimiter=",")
    assert body.shape == (10_001, 43)
    assert wave.shape == (10_001, 2)
    assert coefficients.shape == (6, 3)
    assert mass.shape == (4,)
    for values in (body, wave, coefficients, mass):
        assert np.isfinite(values).all()
    time = np.arange(10_001) * 0.01
    np.testing.assert_allclose(body[:, 0], time, rtol=0, atol=1e-8)
    np.testing.assert_allclose(wave[:, 0], time, rtol=0, atol=1e-8)
    assert np.max(np.abs(wave[:, 1])) > 0.01
    assert np.max(np.abs(body[:, 19])) > 1
    assert np.max(np.abs(coefficients[:, 2])) > 0

    hydro = (Path(APPLICATIONS) / "Mean_Drift/hydroData/sphere.h5").resolve()
    with h5py.File(hydro, "r") as source:
        mean_drift = source["body1/hydro_coeffs/mean_drift/control_surface/val"][:]
    assert np.isfinite(mean_drift).all()
    assert np.max(np.abs(mean_drift)) > 0
