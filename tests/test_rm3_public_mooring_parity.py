"""Pair a configured Python RM3 joint with the published MooringMatrix case."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import ImportedElevationWave, WEC


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="paired MATLAB RM3 MooringMatrix output not provided",
)


def _max_error(actual, expected, limit, label):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


def test_public_joint_with_imported_elevation_and_surge_mooring():
    root = Path(APPLICATIONS)
    reference = Path(REFERENCE)
    hydro = "_Common_Input_Files/RM3/hydroData/rm3.h5"
    wec = WEC("RM3 MooringMatrix")
    float_body = wec.body("float", hydro,
                          inertia=(0, 21_306_090.66, 0))
    spar = wec.body("spar", hydro,
                    inertia=(0, 94_407_091.24, 0))
    wec.floating_joint(
        float_body, spar, damping=1_200_000,
        mooring_surge_stiffness=100_000,
    )
    result = wec.run(
        ImportedElevationWave(
            "Mooring/MooringMatrix/etaData.mat",
            reapply_force_ramp=True,
        ),
        dt=0.01, end_time=400, ramp_time=40,
        radiation_memory=60,
        initial_coordinate={"spar_heave": -0.21},
        base_dir=root,
    )
    wave = np.loadtxt(reference / "RM3_MOORING_MATRIX_wave.csv",
                      delimiter=",")
    _max_error(result.wave_elevation, wave[:, 1], 1e-10, "wave elevation")
    for index, name in enumerate(("float", "spar"), start=1):
        saved = np.loadtxt(
            reference / f"RM3_MOORING_MATRIX_MooringMatrix_body{index}.csv",
            delimiter=",")
        _max_error(result.time, saved[:, 0], 1e-8, "time")
        for dof, position_limit, speed_limit in (
            (0, 0.015, 0.006),
            (2, 0.004, 0.0006),
            (4, 0.0009, 0.0004),
        ):
            _max_error(result.bodies[name].position[:, dof],
                       saved[:, 1 + dof], position_limit,
                       f"{name} position DOF {dof}")
            _max_error(result.bodies[name].velocity[:, dof],
                       saved[:, 7 + dof], speed_limit,
                       f"{name} velocity DOF {dof}")
    source_mooring = np.loadtxt(
        reference / "RM3_MOORING_MATRIX_mooring1.csv", delimiter=",")
    outputs = dict(result.raw.extra_outputs)
    _max_error(outputs["mooring_surge_force"], source_mooring[:, 13],
               1_500, "mooring surge force")
    source_pto = np.loadtxt(
        reference / "RM3_MOORING_MATRIX_MooringMatrix_pto1.csv",
        delimiter=",")
    pto = result.ptos["relative_heave"]
    _max_error(pto.stroke - pto.stroke[0], source_pto[:, 3],
               0.004, "PTO stroke from initial datum")
    _max_error(pto.velocity, source_pto[:, 9], 0.0005, "PTO velocity")
    _max_error(pto.force, source_pto[:, 15], 500, "PTO force")
    _max_error(pto.absorbed_power, -source_pto[:, 21], 500,
               "PTO absorbed power")
