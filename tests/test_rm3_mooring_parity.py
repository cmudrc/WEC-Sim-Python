"""Pair the published RM3 MooringMatrix motion and PTO with MATLAB."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim.bodyClass import BodyClass
from wecsim.rm3Regular import solve_rm3_regular


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="paired MATLAB mooring output and Applications checkout not provided",
)


def _max_error(actual, expected, limit, label):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = np.max(np.abs(actual - expected))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit:.6g}"


def test_rm3_mooring_matrix_full_duration_against_matlab():
    apps = Path(APPLICATIONS)
    reference = Path(REFERENCE)
    prefix = "RM3_MOORING_MATRIX"
    hydro = apps / "_Common_Input_Files/RM3/hydroData/rm3.h5"
    dt, end_time, ramp_time = 0.01, 400.0, 40.0
    time = np.arange(40_001) * dt
    input_wave = loadmat(apps / "Mooring/MooringMatrix/etaData.mat")["etaData"]
    ramp = np.ones_like(time)
    early = time < ramp_time
    ramp[early] = (1 - np.cos(np.pi * time[early] / ramp_time)) / 2
    elevation = np.interp(time, input_wave[:, 0], input_wave[:, 1]) * ramp
    force = np.zeros((len(time), 2, 6))
    for number in (1, 2):
        body = BodyClass(str(hydro))
        body.bodyNumber = number
        body.bodyTotal = 2
        body.readH5file()
        body.hydroForce["userDefinedFe"] = np.zeros((len(time), 6))
        body.userDefinedExcitation(np.vstack((time, elevation)),
                                   dt, [0], 1000, 9.81)
        # The source body block applies the wave ramp a second time.
        force[:, number - 1] = body.hydroForce["userDefinedFe"] * ramp[:, None]

    result = solve_rm3_regular(
        hydro, wave_height=0, radiation_memory=60,
        excitation_force=force, mooring_surge_stiffness=100_000,
        initial_coordinate=np.array([0.0, 0.0, -0.21, 0.0]),
        dt=dt, end_time=end_time, ramp_time=ramp_time,
    )
    bodies = [np.loadtxt(reference / f"{prefix}_MooringMatrix_body{number}.csv",
                         delimiter=",") for number in (1, 2)]
    for index, saved in enumerate(bodies):
        np.testing.assert_allclose(result.time, saved[:, 0], rtol=0, atol=1e-8)
        for dof, name, position_limit, velocity_limit in (
            (0, "surge", 0.015, 0.006),
            (2, "heave", 0.004, 0.0006),
            (4, "pitch", 0.0009, 0.0004),
        ):
            _max_error(result.body_position[:, index, dof], saved[:, 1 + dof],
                       position_limit, f"body{index+1} {name} position")
            _max_error(result.body_velocity[:, index, dof], saved[:, 7 + dof],
                       velocity_limit, f"body{index+1} {name} velocity")

    mooring = np.loadtxt(reference / f"{prefix}_mooring1.csv", delimiter=",")
    _max_error(result.mooring_surge_position, mooring[:, 1], 0.015,
               "mooring surge position")
    _max_error(result.mooring_surge_force, mooring[:, 13], 1_500,
               "mooring surge force")

    pto = np.loadtxt(reference / f"{prefix}_MooringMatrix_pto1.csv",
                     delimiter=",")
    # Source PTO translation is relative to its initialized frame. Recover
    # that datum from the reported body centers and common pitch.
    initial_gap = bodies[0][0, 3] - bodies[1][0, 3]
    pitch = result.body_position[:, 0, 4]
    pitch_speed = result.body_velocity[:, 0, 4]
    stroke = ((result.body_position[:, 0, 2]
               - result.body_position[:, 1, 2]) / np.cos(pitch) - initial_gap)
    speed = ((result.body_velocity[:, 0, 2]
              - result.body_velocity[:, 1, 2]
              + (initial_gap + stroke) * np.sin(pitch) * pitch_speed)
             / np.cos(pitch))
    _max_error(stroke, pto[:, 3], 0.004, "PTO stroke")
    _max_error(speed, pto[:, 9], 0.0005, "PTO speed")
    _max_error(result.pto_force, pto[:, 15], 500, "PTO force")
    _max_error(result.pto_force * speed, pto[:, 21], 500, "PTO power")
