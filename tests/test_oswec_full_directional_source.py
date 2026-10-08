"""Check the pinned OSWEC full-directional source realization."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE
         and os.environ.get("WEC_SIM_REFERENCE_MODEL") == "OSWEC_FULL_DIR"),
    reason="paired MATLAB full-directional output not provided",
)


def test_imported_directional_spectrum_and_wave_replay():
    reference = Path(REFERENCE)
    spectrum_file = (Path(APPLICATIONS) / "Full_Directional_Waves"
                     / "fullDirSpectrum.mat")
    imported = loadmat(spectrum_file)
    frequency = np.loadtxt(reference / "OSWEC_FULL_DIR_frequency.csv", delimiter=",")
    heading = np.loadtxt(reference / "OSWEC_FULL_DIR_directions.csv", delimiter=",")
    spread = np.loadtxt(reference / "OSWEC_FULL_DIR_spread.csv", delimiter=",")
    amplitude = np.loadtxt(reference / "OSWEC_FULL_DIR_amplitude.csv", delimiter=",")
    phase = np.loadtxt(reference / "OSWEC_FULL_DIR_phase.csv", delimiter=",")
    wave = np.loadtxt(reference / "OSWEC_FULL_DIR_wave.csv", delimiter=",")

    assert frequency.shape == (47, 3)
    assert heading.shape == (180, 2)
    assert spread.shape == amplitude.shape == phase.shape == (47, 180)
    assert wave.shape == (8001, 2)
    assert np.isfinite(phase).all() and np.isfinite(wave).all()
    np.testing.assert_allclose(frequency[:, 0], imported["frequencies"].ravel() * 2 * np.pi,
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(frequency[:, 2], imported["spectrum"].ravel(),
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(heading[:, 0], imported["directions"].ravel(),
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(spread, imported["spread"], rtol=0, atol=1e-12)
    np.testing.assert_allclose(amplitude, spread * frequency[:, 2, None] / np.pi,
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose((spread * heading[None, :, 1]).sum(axis=1), 1,
                               rtol=0, atol=1e-12)

    # Replay sparse time samples independently from MATLAB's saved phase.
    indices = np.linspace(0, len(wave) - 1, 41, dtype=int)
    time = wave[indices, 0]
    height = np.sqrt(amplitude * frequency[:, 1, None] * heading[None, :, 1])
    elevation = (height[None] * np.cos(
        time[:, None, None] * frequency[None, :, 0, None]
        + phase[None],
    )).sum(axis=(1, 2))
    ramp = np.ones(len(time))
    early = time < 100
    ramp[early] = (1 - np.cos(np.pi * time[early] / 100)) / 2
    np.testing.assert_allclose(elevation * ramp, wave[indices, 1],
                               rtol=0, atol=1e-11)

    for number in (1, 2):
        body = np.loadtxt(
            reference / f"OSWEC_FULL_DIR_Full_Directional_Waves_body{number}.csv",
            delimiter=",",
        )
        assert body.shape == (8001, 25)
        assert np.isfinite(body).all()
        if number == 2:
            np.testing.assert_allclose(body[:, 1:13], body[0, 1:13],
                                       rtol=0, atol=1e-10)
    pto = np.loadtxt(
        reference / "OSWEC_FULL_DIR_Full_Directional_Waves_pto1.csv",
        delimiter=",",
    )
    assert pto.shape == (8001, 25)
    assert np.isfinite(pto).all()
