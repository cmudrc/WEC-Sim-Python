"""Independent floating OWC body, mooring, chamber, and turbine comparison."""

import os
from pathlib import Path
import shutil

import numpy as np
import pytest

from wecsim import MoorDyn, solve_floating_owc


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
LIBRARY = os.environ.get("WEC_SIM_MOORDYN_LIBRARY")
HYDRO = os.environ.get("WEC_SIM_FLOATING_OWC_H5")
pytestmark = pytest.mark.skipif(
    not (REFERENCE and APPLICATIONS and LIBRARY and HYDRO),
    reason="pinned MATLAB application and native MoorDyn not provided",
)


def _read(name):
    return np.loadtxt(
        Path(REFERENCE) / f"FLOATING_OWC_SOURCE_{name}.csv",
        delimiter=",", ndmin=2,
    )


def _bound(label, predicted, expected, maximum):
    error = np.abs(predicted - expected)
    peak = float(np.max(error))
    rms = float(np.sqrt(np.mean(error**2)))
    assert peak < maximum, f"{label}: maximum {peak:.6g} exceeds {maximum}"
    print(f"{label}: max {peak:.6g}, RMS {rms:.6g}")


def test_published_floating_owc_independent_coupled_motion(tmp_path):
    input_dir = tmp_path / "Mooring"
    input_dir.mkdir()
    input_file = input_dir / "lines.txt"
    shutil.copyfile(
        Path(APPLICATIONS) / "OWC/FloatingOWC/Mooring/lines.txt",
        input_file,
    )
    response = solve_floating_owc(HYDRO, MoorDyn(LIBRARY, input_file))
    floater = _read("FloatingOWC_body1")
    column = _read("FloatingOWC_body2")
    pto = _read("FloatingOWC_pto1")
    pressure = _read("x_deltaP")
    rotor = _read("x_vTurb")
    turbine_power = _read("P_turb")
    pneumatic_power = _read("P_pneumatic")
    coupling = _read("coupling")

    assert response.time.shape == (50_001,)
    np.testing.assert_allclose(response.time, floater[:, 0], rtol=0, atol=1e-8)
    _bound("floater translation (m)", response.floater_pose[:, :3],
           floater[:, 1:4], 0.01)
    _bound("floater rotation (rad)", response.floater_pose[:, 3:6],
           floater[:, 4:7], 0.0002)
    _bound("floater linear velocity (m/s)", response.floater_velocity[:, :3],
           floater[:, 7:10], 0.005)
    _bound("column translation (m)", response.column_pose[:, :3],
           column[:, 1:4], 0.005)
    _bound("column linear velocity (m/s)", response.column_velocity[:, :3],
           column[:, 7:10], 0.003)
    _bound("slider stroke (m)", response.stroke, pto[:, 3], 0.005)
    _bound("chamber gauge pressure (Pa)", response.chamber_pressure,
           pressure[:, 1], 15)
    _bound("rotor speed (rad/s)", response.turbine_speed, rotor[:, 1], 0.5)
    _bound("turbine load power (W)", response.turbine_power,
           turbine_power[:, 1], 25)
    _bound("pneumatic power (W)", response.pneumatic_power,
           pneumatic_power[:, 1], 75)
    force_peak = np.max(np.abs(coupling[1:, 13:19]), axis=0)
    force_error = np.max(
        np.abs(response.mooring_force[1:] - coupling[1:, 13:19]), axis=0)
    assert np.all(force_error < 0.01 * force_peak), force_error / force_peak
    print("live mooring force max error:", force_error)

    output_file = input_dir / "lines.out"
    with output_file.open() as output:
        output_names = output.readline().split()
    source_names = (
        Path(REFERENCE) / "FLOATING_OWC_SOURCE_fairlead_columns.csv"
    ).read_text().splitlines()
    tension = np.loadtxt(output_file, skiprows=1)[
        :, [output_names.index(name) for name in source_names]
    ]
    source_tension = _read("fairlead_tension")
    assert tension.shape == source_tension.shape == (50_001, 6)
    tension_peak = np.max(np.abs(source_tension[:, 1:]), axis=0)
    tension_error = np.max(
        np.abs(tension[:, 1:] - source_tension[:, 1:]), axis=0)
    assert np.all(tension_error < 0.01 * tension_peak), tension_error / tension_peak
    print("live fairlead tension max error:", tension_error)
