"""Check numeric output from the pinned floating OWC MATLAB application."""

import os
from pathlib import Path

import h5py
import numpy as np
import pytest


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
HYDRO = os.environ.get("WEC_SIM_FLOATING_OWC_H5")
pytestmark = pytest.mark.skipif(
    not (REFERENCE and HYDRO),
    reason="pinned floating OWC MATLAB output and BEMIO HDF5 not provided",
)


def _read(name):
    values = np.loadtxt(Path(REFERENCE) / f"FLOATING_OWC_SOURCE_{name}.csv",
                        delimiter=",")
    assert values.ndim == 2 and np.isfinite(values).all()
    assert np.all(np.diff(values[:, 0]) > 0)
    return values


def test_published_floating_owc_source_has_coupled_motion_and_power():
    wave = _read("wave")
    assert wave.shape[1] == 2
    assert abs(wave[-1, 0] - 500) < 1e-8
    assert np.max(np.abs(wave[:, 1])) > 1

    bodies = [_read(f"FloatingOWC_body{number}") for number in (1, 2)]
    for body in bodies:
        assert body.shape[1] == 25
        np.testing.assert_allclose(body[:, 0], wave[:, 0], rtol=0, atol=1e-8)
    assert np.ptp(bodies[0][:, 3] - bodies[1][:, 3]) > 0.01

    coupling = _read("coupling")
    assert coupling.shape[1] == 19
    assert np.max(np.linalg.norm(coupling[:, 13:16], axis=1)) > 100
    fairlead = _read("fairlead_tension")
    assert fairlead.shape[1] >= 4
    assert np.max(fairlead[:, 1:]) > 100

    for name in ("x_vTurb", "u", "P_pneumatic", "P_turb", "eta_turb"):
        trace = _read(name)
        assert trace.shape[1] >= 2
        assert trace[-1, 0] >= 499.9
    with h5py.File(HYDRO) as hydro:
        assert "body1" in hydro and "body2" in hydro
        assert "hydro_coeffs" in hydro["body1"]
