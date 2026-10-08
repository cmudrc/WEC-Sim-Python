"""Pair the published Sphere MPC waves and prediction matrices with MATLAB."""

import os
from pathlib import Path

import h5py
import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import JONSWAPWave
from wecsim.irregularWave import jonswap_equal_energy_components
from wecsim.mpcPlant import (
    build_sphere_mpc_matrices, replay_sphere_mpc_plant,
)
from wecsim.waveClass import WaveClass


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="published MATLAB Sphere MPC output and HDF5 absent",
)


def _load(name):
    return np.loadtxt(Path(REFERENCE) / f"SPHERE_MPC_{name}.csv",
                      delimiter=",", ndmin=2)


def _match(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


def test_published_sphere_mpc_waves_and_prediction_model():
    apps = Path(APPLICATIONS)
    hydro = apps / "_Common_Input_Files/Sphere/hydroData/sphere.h5"
    coefficients = apps / "Controls/MPC/coeff.mat"
    source = loadmat(Path(REFERENCE) / "SPHERE_MPC_setup.mat",
                     simplify_cells=True)
    matrices = build_sphere_mpc_matrices(hydro, coefficients)
    for name in ("A", "Bu", "Bv", "C"):
        _match(getattr(matrices, name), np.asarray(source["plant"][name]),
               1e-10, f"prediction plant {name}")
    for name in ("Sx", "Su", "Sv", "Q", "H"):
        _match(getattr(matrices, name), np.asarray(source["setup"][name]),
               1e-10, f"MPC {name}")
    assert matrices.Sx.shape == (90, 9)
    assert matrices.Su.shape == (90, 31)
    assert np.linalg.eigvalsh((matrices.H + matrices.H.T) / 2).min() > 0

    components = _load("components")
    assert components.shape == (500, 4)
    phase = components[:, 3, None]
    assert JONSWAPWave(2.5, 8, phase_file="phases.csv").as_case()["type"] == "jonswap"
    sea = jonswap_equal_energy_components(
        hydro, significant_height=2.5, peak_period=8,
        directions=np.array([0]), spreading=np.array([1]), phase=phase,
    )
    _match(sea.omega, components[:, 0], 1e-12, "public JONSWAP frequencies")
    _match(sea.spectral_amplitude, components[:, 1], 1e-12,
           "public JONSWAP spectrum")
    _match(sea.d_omega, components[:, 2], 1e-12, "public JONSWAP widths")
    with h5py.File(hydro) as h5:
        frequencies = np.asarray(h5["/simulation_parameters/w"]).ravel()
    wave = WaveClass("irregular")
    wave.T = 8
    wave.H = 2.5
    wave.spectrumType = "JS"
    wave.freqDisc = "EqualEnergy"
    wave.numFreq = 500
    wave.phaseData = phase
    wave.waveDir = [0]
    wave.waveSpread = [1]
    wave.waveSetup([frequencies.min(), frequencies.max()], "infinite",
                   100, .01, 40000, 9.81, 1000, 400)
    assert wave.gamma == 1.0
    _match(np.asarray(wave.w).ravel(), components[:, 0], 1e-12,
           "JONSWAP frequencies")
    _match(np.asarray(wave.A).ravel(), components[:, 1], 1e-12,
           "JONSWAP amplitudes")
    _match(np.asarray(wave.dw).ravel(), components[:, 2], 1e-12,
           "JONSWAP widths")
    _match(np.asarray(wave.waveAmpTime).T, _load("wave"), 3e-12,
           "JONSWAP elevation")

    controller = _load("controller")
    body = _load("MPC_body1")
    plant = _load("plant_output")
    assert controller.shape == (40001, 13)
    assert body.shape == (40001, 25)
    assert plant.shape == (40001, 4)
    _match(plant[:, 0], body[:, 0], 1e-12, "plant sample time")
    _match(plant[:, 3], controller[:, 3], 1e-8,
           "internal plant PTO force")
    replay = replay_sphere_mpc_plant(
        matrices, controller[:, 3], body[:, 21], dt=.01,
    )
    _match(replay[:, 0], plant[:, 1], 7e-5, "internal plant heave speed")
    _match(replay[:, 1], plant[:, 2], 5e-5, "internal plant heave position")
    _match(controller[:, 9], controller[:, 3] * body[:, 9], 1e-7,
           "logged controller power")
    assert controller[np.flatnonzero(np.abs(controller[:, 3]) > 1e-9)[0], 0] == 205.51
    assert 2e6 < np.max(np.abs(controller[:, 3])) < 2.5e6
    assert np.max(np.abs(np.diff(controller[:, 3]) / .01)) < 1.500001e6
