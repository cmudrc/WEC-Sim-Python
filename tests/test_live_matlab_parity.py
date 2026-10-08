"""Compare production Python output with a pinned, executed MATLAB waveClass."""

import os
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]

from wecsim.waveClass import WaveClass  # noqa: E402

REFERENCE = os.environ.get("WEC_SIM_MATLAB_REFERENCE_DIR")
pytestmark = pytest.mark.skipif(not REFERENCE, reason="MATLAB reference output not provided")


def test_regular_wave_against_executed_matlab():
    wave = WaveClass("regular")
    wave.T = 8.0
    wave.H = 2.5
    wave.waveDir = [0]
    wave.wavegauge1loc = [5, 5]
    wave.wavegauge2loc = [10, 0]
    wave.wavegauge3loc = [0, -10]
    wave.waveSetup(
        [5.19999512307279, 0.0199999977946844],
        "infinite", 100.0, 0.1, 2000, 9.81, 1000.0, 200.0,
    )
    reference = Path(REFERENCE)
    origin = np.loadtxt(reference / "regular_origin.csv", delimiter=",")
    markers = np.loadtxt(reference / "regular_markers.csv", delimiter=",")
    np.testing.assert_allclose(np.asarray(wave.waveAmpTime).T, origin, rtol=0, atol=1e-10)
    for index, attribute in enumerate(("waveAmpTime1", "waveAmpTime2", "waveAmpTime3")):
        np.testing.assert_allclose(
            np.asarray(getattr(wave, attribute))[1], markers[:, index + 1],
            rtol=0, atol=1e-10,
        )


def test_finite_depth_against_executed_matlab():
    wave = WaveClass("regular")
    wave.T = 8.0
    wave.H = 2.5
    wave.waveDir = [0]
    wave.wavegauge1loc = [0, 0]
    wave.wavegauge2loc = [0, 0]
    wave.wavegauge3loc = [0, 0]
    wave.waveSetup([0.4, 2.0], 25.0, 0, 0.1, 2000, 9.81, 1000.0, 200.0)
    expected = np.loadtxt(Path(REFERENCE) / "finite_depth.csv", delimiter=",")
    np.testing.assert_allclose([np.asarray(wave.k).item(), np.asarray(wave.Pw).item()],
                               expected, rtol=1e-11)


@pytest.mark.parametrize("spectrum,height", [("PM", 2.5), ("JS", 4.0)])
@pytest.mark.parametrize("discretization", ["Traditional", "EqualEnergy"])
def test_current_irregular_spectrum_against_executed_matlab(
    spectrum, height, discretization,
):
    wave = WaveClass("irregular")
    wave.T = 8
    wave.H = height
    wave.spectrumType = spectrum
    wave.freqDisc = discretization
    wave.numFreq = 64
    wave.phaseSeed = 1
    if spectrum == "PM":
        wave.waveDir = [0, 30, 90]
        wave.waveSpread = [0.1, 0.2, 0.7]
    wave.wavegauge1loc = [5, 5]
    wave.wavegauge2loc = [10, 0]
    wave.wavegauge3loc = [0, -10]
    reference = Path(REFERENCE)
    label = spectrum.lower() + ("_equal" if discretization == "EqualEnergy" else "")
    wave.phaseData = np.loadtxt(reference / f"{label}_phase.csv",
                                delimiter=",")
    wave.waveSetup([0.4, 2.0], "infinite", 1, 0.1, 20, 9.81, 1000, 2)

    if discretization == "EqualEnergy":
        expected = np.loadtxt(reference / f"{label}_bins.csv", delimiter=",")
        actual = np.column_stack((wave.w, wave.dw, wave.S))
    else:
        expected = np.loadtxt(reference / f"{label}_spectrum.csv", delimiter=",")
        actual = np.column_stack((wave.w, wave.S))
    np.testing.assert_allclose(actual, expected,
                               rtol=2e-12, atol=1e-14)
    power = np.loadtxt(reference / f"{label}_power.csv",
                       delimiter=",")
    np.testing.assert_allclose(wave.Pw, power, rtol=2e-12, atol=1e-9)
    if spectrum == "JS":
        gamma = np.loadtxt(reference / f"{label}_gamma.csv", delimiter=",")
        np.testing.assert_allclose(wave.gamma, gamma, rtol=2e-12)
    source_elevation = np.loadtxt(
        reference / f"{label}_elevation.csv", delimiter=",")
    np.testing.assert_allclose(np.asarray(wave.waveAmpTime).T,
                               source_elevation, rtol=0, atol=2e-12)
    source_markers = np.loadtxt(
        reference / f"{label}_markers.csv", delimiter=",")
    for index, attribute in enumerate(("waveAmpTime1", "waveAmpTime2",
                                       "waveAmpTime3")):
        np.testing.assert_allclose(
            np.asarray(getattr(wave, attribute))[1], source_markers[:, index + 1],
            rtol=0, atol=2e-12,
        )
