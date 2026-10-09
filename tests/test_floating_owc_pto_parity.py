"""Pair a configured nonzero floating OWC slider PTO with MATLAB."""

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
    reason="derived floating OWC MATLAB output and native MoorDyn absent",
)
PREFIX = "FLOATING_OWC_PTO_"


def _read(name):
    return np.loadtxt(Path(REFERENCE) / f"{PREFIX}{name}.csv",
                      delimiter=",", ndmin=2)


def _bound(name, predicted, source, limit):
    error = predicted - source
    maximum = float(np.max(np.abs(error)))
    rms = float(np.sqrt(np.mean(error**2)))
    assert maximum < limit, f"{name}: {maximum:.6g} exceeds {limit}"
    print(f"{name}: maximum {maximum:.6g}, RMS {rms:.6g}")


def test_nonzero_floating_owc_pto_dynamics(tmp_path):
    floater = _read("FloatingOWC_body1")
    column = _read("FloatingOWC_body2")
    pto = _read("FloatingOWC_pto1")
    pressure = _read("deltaP")
    rotor = _read("vTurb")
    turbine_power = _read("turbine_power")
    coupling = _read("coupling")
    wave = _read("wave")
    assert floater.shape == column.shape == pto.shape == (15_001, 25)
    assert pressure.shape == rotor.shape == turbine_power.shape == wave.shape == (15_001, 2)
    assert coupling.shape == (15_001, 19)
    np.testing.assert_allclose(floater[:, 0], np.arange(15_001) * .01,
                               rtol=0, atol=1e-8)
    source_pto_force = pto[:, 15]
    assert np.max(np.abs(source_pto_force)) > 10_000
    _bound("MATLAB PTO spring and damper law (N)", source_pto_force,
           -15_000 * pto[:, 3] - 60_000 * pto[:, 9], 1e-5)
    _bound("MATLAB PTO mechanical power (W)", pto[:, 21],
           source_pto_force * pto[:, 9], 1e-5)

    input_dir = tmp_path / "Mooring"
    input_dir.mkdir()
    input_file = input_dir / "lines.txt"
    shutil.copyfile(Path(APPLICATIONS) / "OWC/FloatingOWC/Mooring/lines.txt",
                    input_file)
    response = solve_floating_owc(
        HYDRO, MoorDyn(LIBRARY, input_file), end_time=150,
        pto_stiffness=15_000, pto_damping=60_000,
    )
    np.testing.assert_allclose(response.time, floater[:, 0], rtol=0, atol=1e-8)
    _bound("floater position (m)", response.floater_pose[:, :3],
           floater[:, 1:4], .015)
    _bound("floater rotation (rad)", response.floater_pose[:, 3:6],
           floater[:, 4:7], 3e-4)
    _bound("column position (m)", response.column_pose[:, :3],
           column[:, 1:4], .01)
    _bound("PTO stroke (m)", response.stroke, pto[:, 3], .015)
    _bound("PTO force (N)", response.pto_force, source_pto_force, 1_500)
    _bound("PTO mechanical power (W)", -response.pto_mechanical_power,
           pto[:, 21], 3_000)
    _bound("chamber pressure (Pa)", response.chamber_pressure,
           pressure[:, 1], 20)
    _bound("rotor speed (rad/s)", response.turbine_speed, rotor[:, 1], 1)
    _bound("turbine load power (W)", response.turbine_power,
           turbine_power[:, 1], 50)
    force_peak = np.max(np.abs(coupling[1:, 13:19]), axis=0)
    force_error = np.max(np.abs(
        response.mooring_force[1:] - coupling[1:, 13:19]), axis=0)
    assert np.all(force_error < .02 * force_peak), force_error / force_peak
