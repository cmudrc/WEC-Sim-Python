"""Published fixed hydro monopile against pinned MATLAB WEC-Sim."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim.irregularWave import (
    pm_equal_energy_components, synthesize_irregular_response,
)


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
HYDRO = os.environ.get("WEC_SIM_MONOPILE_H5")
pytestmark = pytest.mark.skipif(
    not REFERENCE or not HYDRO,
    reason="pinned hydro monopile MATLAB output and HDF5 absent",
)


def _source(name):
    return np.loadtxt(Path(REFERENCE) / f"MONOPILE_HYDRO_{name}.csv",
                      delimiter=",", ndmin=2)


def test_published_fixed_monopile_wave_and_excitation():
    components, wave = _source("components"), _source("wave")
    body1, body2 = _source("monopile_body1"), _source("monopile_body2")
    assert components.shape[1] == 4 and len(components) >= 100
    assert wave.shape == (40_001, 2)
    assert body1.shape == body2.shape == (40_001, 25)
    assert np.isfinite(components).all()
    assert np.isfinite(wave).all()
    for body in (body1, body2):
        np.testing.assert_allclose(body[:, 0], wave[:, 0],
                                   rtol=0, atol=1e-10)
        assert np.max(np.abs(body[:, 1:7] - body[0, 1:7])) < 1e-6
        assert np.max(np.abs(body[:, 7:13])) < 1e-6

    sea = pm_equal_energy_components(
        HYDRO, significant_height=2, peak_period=5,
        directions=(0,), spreading=(1,), count=len(components),
        phase=components[:, 3:4],
    )
    np.testing.assert_allclose(sea.omega, components[:, 0],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(sea.spectral_amplitude, components[:, 1],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(sea.d_omega, components[:, 2],
                               rtol=0, atol=1e-12)
    result = synthesize_irregular_response(
        HYDRO, sea, dt=.01, end_time=400, ramp_time=100,
        rho=1025, g=9.81,
    )
    np.testing.assert_allclose(result.elevation, wave[:, 1],
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.excitation_force[:, :6],
                               body1[:, 19:25], rtol=0, atol=1e-5)
    assert np.max(np.abs(body2[:, 19:25])) < 1e-7
