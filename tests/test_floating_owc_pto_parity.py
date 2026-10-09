"""Pair a configured nonzero floating OWC slider PTO with MATLAB."""

import os
from pathlib import Path
import shutil

import numpy as np
import pytest

from wecsim import MoorDyn, RegularWave, WEC


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
    wec = WEC("Floating OWC")
    float_body = wec.body(
        "floater", HYDRO, hydro_body=1,
        inertia=(1.531e9, 1.531e9, 0.1118e9),
    )
    column_body = wec.body(
        "column", HYDRO, hydro_body=2, mass=4_493_450,
    )
    wec.floating_owc(
        float_body, column_body, moordyn=MoorDyn(LIBRARY, input_file),
        moordyn_point=float_body.at(0, 0, 31.945),
        column_height=50.69, column_diameter=5.89,
        pto_name="turbine_slider", pto_stiffness=15_000,
        pto_damping=60_000,
    )
    response = wec.run(RegularWave(4.5, 11.2), dt=.01,
                       end_time=150, ramp_time=50)
    extras = dict(response.raw.extra_outputs)
    np.testing.assert_allclose(response.time, floater[:, 0], rtol=0, atol=1e-8)
    _bound("wave elevation (m)", response.wave_elevation, wave[:, 1], 1e-10)
    np.testing.assert_allclose(response.coordinates["column_stroke"].position,
                               response.ptos["turbine_slider"].stroke,
                               rtol=0, atol=0)
    _bound("floater position (m)", response.bodies["floater"].position[:, :3],
           floater[:, 1:4], .007)
    _bound("floater rotation (rad)", response.bodies["floater"].position[:, 3:6],
           floater[:, 4:7], 1e-4)
    _bound("column position (m)", response.bodies["column"].position[:, :3],
           column[:, 1:4], .003)
    _bound("PTO stroke (m)", response.ptos["turbine_slider"].stroke,
           pto[:, 3], .008)
    _bound("PTO force (N)", response.ptos["turbine_slider"].force,
           source_pto_force, 400)
    _bound("PTO mechanical power (W)", -extras["pto_mechanical_power"],
           pto[:, 21], 1_500)
    _bound("chamber pressure (Pa)", extras["chamber_pressure"],
           pressure[:, 1], 20)
    _bound("rotor speed (rad/s)", extras["turbine_speed"], rotor[:, 1], .8)
    _bound("turbine load power (W)", extras["turbine_power"],
           turbine_power[:, 1], 50)
    force_peak = np.max(np.abs(coupling[1:, 13:19]), axis=0)
    force_error = np.max(np.abs(
        extras["moordyn_connection_force"][1:] - coupling[1:, 13:19]), axis=0)
    assert np.all(force_error < .015 * force_peak), force_error / force_peak
