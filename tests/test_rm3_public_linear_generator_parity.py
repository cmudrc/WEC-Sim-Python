"""Public body-local PTO configuration against the published RM3 PTO-Sim run."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import DirectLinearGenerator, RegularWave, WEC


RM3_H5 = os.environ.get("WEC_SIM_RM3_H5")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (RM3_H5 and REFERENCE),
    reason="paired MATLAB RM3 direct-generator output and HDF5 absent",
)


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    error = float(np.max(np.abs(actual - expected)))
    assert np.isfinite(error) and error < limit, (
        f"{label}: {error:.6g} exceeds {limit}"
    )


def test_configured_two_body_generator_matches_matlab():
    folder = Path(REFERENCE)
    bodies = [
        np.loadtxt(folder / f"RM3_DD_PTO_RM3_DD_PTO_body{index}.csv",
                   delimiter=",")
        for index in (1, 2)
    ]
    pto = np.loadtxt(folder / "RM3_DD_PTO_RM3_DD_PTO_pto1.csv",
                     delimiter=",")
    drive = np.loadtxt(folder / "RM3_DD_PTO_drive.csv", delimiter=",")
    assert all(body.shape == (40_001, 25) for body in bodies)
    assert pto.shape == (40_001, 25)
    assert drive.shape == (40_001, 12)

    wec = WEC("RM3 with linear generator")
    floating = wec.body("float", RM3_H5)
    spar = wec.body("spar", RM3_H5)
    wec.coordinate("float_heave", floating.move("heave"))
    wec.coordinate("spar_heave", spar.move("heave"))
    wec.pto(
        "PTO1", spar.at(0, 0, 0), floating.at(0, 0, 0),
        axis=(0, 0, 1),
        linear_generator=DirectLinearGenerator(
            stator_resistance=4.58, friction=-100, pole_pitch=.072,
            magnet_flux=8, inductance=.285, load_resistance=-117.6471,
        ),
    )
    result = wec.run(RegularWave(height=2.5, period=8), dt=.0005,
                     end_time=400, ramp_time=100)
    sample = slice(None, None, 20)
    history = result.ptos["PTO1"]
    generator = history.linear_generator
    assert generator is not None
    _max_error(result.time[sample], drive[:, 0], 1e-9, "time")
    for index, name in enumerate(("float", "spar")):
        motion = result.bodies[name]
        _max_error(motion.position[sample, 2], bodies[index][:, 3],
                   1e-5, f"{name} heave")
        _max_error(motion.velocity[sample, 2], bodies[index][:, 9],
                   1e-5, f"{name} speed")
    _max_error(history.stroke[sample], pto[:, 3], 1e-5, "PTO stroke")
    _max_error(history.velocity[sample], pto[:, 9], 1e-5, "PTO speed")
    _max_error(history.force[sample], drive[:, 2], .005, "generator force")
    _max_error(history.absorbed_power[sample], drive[:, 1],
               .005, "absorbed power")
    _max_error(generator.friction_force[sample], drive[:, 3],
               .001, "friction force")
    _max_error(generator.electrical_power[sample], drive[:, 11],
               .005, "electrical power")
    _max_error(generator.phase_current[sample], drive[:, 4:7],
               .0001, "phase current")
    _max_error(generator.phase_voltage[sample], drive[:, 7:10],
               .01, "phase voltage")
