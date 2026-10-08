"""Check the published WaveBot system-identification input and MATLAB export."""

import os
from pathlib import Path

import numpy as np
import pytest


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="pinned WaveBot MATLAB output is absent",
)


def _read(name):
    values = np.loadtxt(
        Path(REFERENCE) / f"WAVEBOT_IMPEDANCE_{name}.csv",
        delimiter=",", ndmin=2,
    )
    assert np.isfinite(values).all(), name
    return values


def test_published_multisine_and_mooring_force_laws():
    input_signal = _read("multisine")
    actuation = _read("actuation")
    parameters = _read("parameters")
    mooring = _read("mooring")
    body = _read("CalcImpedance_body1")
    forces = _read("body1_forces")

    assert input_signal.shape == (149_999, 4)
    assert actuation.shape[1] == 4
    assert mooring.shape == (60_001, 19)
    assert body.shape == (60_001, 25)
    assert forces.shape == (60_001, 37)
    np.testing.assert_allclose(body[:, 0], np.arange(60_001) * .01,
                               rtol=0, atol=1e-8)
    np.testing.assert_allclose(mooring[:, 0], body[:, 0],
                               rtol=0, atol=1e-8)
    np.testing.assert_allclose(forces[:, 0], body[:, 0],
                               rtol=0, atol=1e-8)
    np.testing.assert_allclose(parameters[0, :4], [1156.5, 84, 84, 84])
    np.testing.assert_allclose(parameters[0, 4:], [50, 100, 4])

    commanded = np.column_stack([
        np.interp(actuation[:, 0], input_signal[:, 0], input_signal[:, i + 1])
        for i in range(3)
    ]) * parameters[0, 4:]
    np.testing.assert_allclose(actuation[:, 1:], commanded,
                               rtol=0, atol=1e-8)

    stiffness = np.diag([24_000, 24_000, 5_000, 0, 1_000, 0])
    damping = np.diag([5, 0, 5, 0, 0, 0])
    spring_damper = (-mooring[:, 1:7] @ stiffness.T
                     - mooring[:, 7:13] @ damping.T)
    np.testing.assert_allclose(mooring[:, 13:19], spring_damper,
                               rtol=0, atol=1e-5)
    assert np.max(np.abs(mooring[:, 13:19])) > 100
    assert np.max(np.abs(body[:, 1:7])) > .001


def test_published_hydrodynamic_matrices_and_body_force_channels():
    added_mass = _read("added_mass")
    hydrostatic = _read("hydrostatic")
    forces = _read("body1_forces")
    body = _read("CalcImpedance_body1")
    assert added_mass.shape == hydrostatic.shape == (6, 6)
    assert np.max(np.abs(added_mass)) > 0
    assert np.max(np.abs(hydrostatic)) > 0
    np.testing.assert_allclose(body[:, 19:25], 0, rtol=0, atol=1e-10)
    linear_damping = np.diag([1000, 0, 1000, 0, 100, 0])
    np.testing.assert_allclose(
        forces[:, 25:31], -body[:, 7:13] @ linear_damping.T,
        rtol=0, atol=1e-6,
    )
