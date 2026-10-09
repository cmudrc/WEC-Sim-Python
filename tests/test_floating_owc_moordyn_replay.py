"""Pair the published floating OWC's nine-line MoorDyn response."""

import os
from pathlib import Path
import shutil

import h5py
import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from wecsim import MoorDyn


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
LIBRARY = os.environ.get("WEC_SIM_MOORDYN_LIBRARY")
HYDRO = os.environ.get("WEC_SIM_FLOATING_OWC_H5")
pytestmark = pytest.mark.skipif(
    not (REFERENCE and APPLICATIONS and LIBRARY and HYDRO),
    reason="pinned floating OWC output, application, MoorDyn, and HDF5 absent",
)


def _read(name):
    return np.loadtxt(
        Path(REFERENCE) / f"FLOATING_OWC_SOURCE_{name}.csv",
        delimiter=",", ndmin=2,
    )


def test_published_nine_line_moordyn_from_floater_motion(tmp_path):
    """Advance native MoorDyn from the body state, not its saved force."""
    floater = _read("FloatingOWC_body1")
    coupling = _read("coupling")
    expected_tension = _read("fairlead_tension")
    assert floater.shape == (50_001, 25)
    assert coupling.shape == (50_001, 19)
    np.testing.assert_allclose(floater[:, 0], coupling[:, 0], atol=1e-8)

    with h5py.File(HYDRO) as hydro:
        center = np.asarray(hydro["body1/properties/cg"]).reshape(3)
    # The published MoorDyn coupled body is at the world origin in equilibrium.
    # Its attachment is therefore -CG in the floater's body coordinates.
    offset = Rotation.from_euler("xyz", floater[:, 4:7]).apply(
        np.broadcast_to(-center, (len(floater), 3))
    )
    pose = floater[:, 1:7].copy()
    pose[:, :3] += offset
    velocity = floater[:, 7:13].copy()
    velocity[:, :3] += np.cross(velocity[:, 3:6], offset)
    assert np.max(np.abs(pose - coupling[:, 1:7])) < 1e-10
    assert np.max(np.abs(velocity - coupling[:, 7:13])) < 1e-10

    input_dir = tmp_path / "Mooring"
    input_dir.mkdir()
    input_file = input_dir / "lines.txt"
    shutil.copyfile(
        Path(APPLICATIONS) / "OWC/FloatingOWC/Mooring/lines.txt",
        input_file,
    )
    force = np.zeros((len(coupling), 6))
    with MoorDyn(LIBRARY, input_file).start(pose[0], velocity[0]) as mooring:
        for i in range(1, len(coupling)):
            force[i] = mooring.step(
                pose[i], velocity[i],
                coupling[i - 1, 0], coupling[i, 0] - coupling[i - 1, 0],
            )

    expected_force = coupling[1:, 13:19]
    force_error = force[1:] - expected_force
    force_peak = np.max(np.abs(expected_force), axis=0)
    force_max = np.max(np.abs(force_error), axis=0)
    force_rms = np.sqrt(np.mean(force_error**2, axis=0))
    assert np.all(force_max < 0.001 * force_peak), force_max / force_peak
    assert np.all(force_rms < 0.0001 * force_peak), force_rms / force_peak

    # The input repeats FairTen3. MATLAB's struct has unique field names, so
    # its exported fairlead trace contains one copy of that channel.
    output_file = input_dir / "lines.out"
    with output_file.open() as output:
        output_names = output.readline().split()
    source_names = (
        Path(REFERENCE) / "FLOATING_OWC_SOURCE_fairlead_columns.csv"
    ).read_text().splitlines()
    columns = [output_names.index(name) for name in source_names]
    tension = np.loadtxt(output_file, skiprows=1)[:, columns]
    assert tension.shape == expected_tension.shape == (50_001, 6)
    np.testing.assert_allclose(tension[:, 0], expected_tension[:, 0], atol=1e-8)
    tension_error = tension[:, 1:] - expected_tension[:, 1:]
    tension_peak = np.max(np.abs(expected_tension[:, 1:]), axis=0)
    tension_max = np.max(np.abs(tension_error), axis=0)
    tension_rms = np.sqrt(np.mean(tension_error**2, axis=0))
    assert np.all(tension_max < 0.001 * tension_peak), tension_max / tension_peak
    assert np.all(tension_rms < 0.0001 * tension_peak), tension_rms / tension_peak
    print(f"floating OWC MoorDyn force max/rms: {force_max} / {force_rms}; "
          f"tension max/rms: {tension_max} / {tension_rms}")
