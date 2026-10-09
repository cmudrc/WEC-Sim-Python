"""Replay the published RM3 coupling signal through Python's MoorDyn binding."""

import os
from pathlib import Path
import shutil

import numpy as np
import pytest

from wecsim import MoorDyn


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
LIBRARY = os.environ.get("WEC_SIM_MOORDYN_LIBRARY")
pytestmark = pytest.mark.skipif(
    not (REFERENCE and APPLICATIONS and LIBRARY),
    reason="pinned MATLAB MoorDyn output, input, and native library not provided",
)


def test_pinned_moordyn_forces_from_source_kinematics(tmp_path):
    reference = Path(REFERENCE)
    coupling = np.loadtxt(reference / "RM3_MOORDYN_coupling.csv", delimiter=",")
    fairlead = np.loadtxt(reference / "RM3_MOORDYN_fairlead_tension.csv",
                          delimiter=",")
    input_dir = tmp_path / "Mooring"
    input_dir.mkdir()
    input_file = input_dir / "lines.txt"
    shutil.copyfile(Path(APPLICATIONS) / "Mooring/MoorDyn/Mooring/lines.txt",
                    input_file)

    force = np.zeros((len(coupling), 6))
    with MoorDyn(LIBRARY, input_file).start(coupling[0, 1:7],
                                            coupling[0, 7:13]) as mooring:
        for i in range(1, len(coupling)):
            force[i] = mooring.step(
                coupling[i, 1:7], coupling[i, 7:13],
                coupling[i - 1, 0], coupling[i, 0] - coupling[i - 1, 0],
            )

    line_output = np.loadtxt(input_dir / "lines.out", skiprows=1)
    assert line_output.shape == fairlead.shape == (40_001, 4)
    np.testing.assert_allclose(line_output[:, 0], fairlead[:, 0],
                               rtol=0, atol=1e-7)
    tension_error = line_output[:, 1:] - fairlead[:, 1:]
    tension_peak = np.max(np.abs(fairlead[:, 1:]))
    tension_max = np.max(np.abs(tension_error))
    tension_rms = np.sqrt(np.mean(tension_error**2))
    assert tension_max < 0.005 * tension_peak, tension_max
    assert tension_rms < 0.0001 * tension_peak, tension_rms

    expected_force = coupling[1:, 13:19]
    force_error = force[1:] - expected_force
    for axis in (0, 2, 4):
        peak = np.max(np.abs(expected_force[:, axis]))
        maximum = np.max(np.abs(force_error[:, axis]))
        rms = np.sqrt(np.mean(force_error[:, axis]**2))
        assert maximum < 0.01 * peak, (axis, maximum)
        assert rms < 0.001 * peak, (axis, rms)
    for axis in (1, 3, 5):
        assert np.max(np.abs(force_error[:, axis])) < 1e-3, axis
    print(f"MoorDyn replay: tension max/rms {tension_max:.3f}/{tension_rms:.3f} N; "
          f"force max {np.max(np.abs(force_error), axis=0)}")
