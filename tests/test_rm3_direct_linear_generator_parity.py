"""Paired 400 s RM3 PTO-Sim direct linear generator trajectory."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import DirectLinearGenerator, run_rm3_direct_linear_generator


RM3_H5 = os.environ.get("WEC_SIM_RM3_H5")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (RM3_H5 and REFERENCE),
    reason="paired MATLAB RM3 direct-generator output and HDF5 absent",
)


@pytest.fixture(scope="module")
def paired():
    folder = Path(REFERENCE)
    bodies = [
        np.loadtxt(folder / f"RM3_DD_PTO_RM3_DD_PTO_body{body}.csv",
                   delimiter=",")
        for body in (1, 2)
    ]
    pto = np.loadtxt(folder / "RM3_DD_PTO_RM3_DD_PTO_pto1.csv",
                     delimiter=",")
    drive = np.loadtxt(folder / "RM3_DD_PTO_drive.csv", delimiter=",")
    assert all(body.shape == (40_001, 25) for body in bodies)
    assert pto.shape == (40_001, 25)
    assert drive.shape == (40_001, 12)
    generator = DirectLinearGenerator(
        stator_resistance=4.58, friction=-100, pole_pitch=0.072,
        magnet_flux=8, inductance=0.285, load_resistance=-117.6471,
    )
    response = run_rm3_direct_linear_generator(
        RM3_H5, generator=generator,
    )
    return response, bodies, pto, drive


def _max_error(actual, reference, limit, label):
    assert actual.shape == reference.shape, label
    assert np.isfinite(actual).all(), label
    error = float(np.max(np.abs(actual - reference)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


def test_active_body_and_pto_motion(paired):
    response, bodies, pto, drive = paired
    _max_error(response.time, drive[:, 0], 1e-9, "time")
    for index, body in enumerate(bodies):
        _max_error(response.body_heave[:, index], body[:, 3],
                   1e-5, f"body {index + 1} heave")
        _max_error(response.body_heave_velocity[:, index], body[:, 9],
                   1e-5, f"body {index + 1} heave velocity")
    _max_error(response.pto_stroke, pto[:, 3], 1e-5, "PTO stroke")
    _max_error(response.pto_velocity, pto[:, 9], 1e-5, "PTO velocity")
    _max_error(response.pto_force, drive[:, 2], .005, "generator force")
    _max_error(response.absorbed_power, drive[:, 1],
               .005, "absorbed mechanical power")


def test_generator_electrical_output(paired):
    response, _, _, drive = paired
    _max_error(response.friction_force, drive[:, 3],
               .001, "friction force")
    _max_error(response.electrical_power, drive[:, 11],
               .005, "electrical power")
    _max_error(response.phase_current, drive[:, 4:7],
               .0001, "three-phase currents")
    _max_error(response.phase_voltage, drive[:, 7:10],
               .01, "three-phase voltages")
