"""Run JSON dynamics cases through the public CLI against paired MATLAB data."""

import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
MODEL = os.environ.get("WEC_SIM_REFERENCE_MODEL")
CORE = os.environ.get("WEC_SIM_MATLAB_CORE_DIR")
APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
SPHERE_H5 = os.environ.get("WEC_SIM_SPHERE_H5")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")


def _run(tmp_path, case):
    case_file = tmp_path / "case.json"
    output = tmp_path / "motion.csv"
    case_file.write_text(json.dumps(case), encoding="utf-8")
    subprocess.run(
        [sys.executable, "-m", "source.objects.wecSimPython",
         str(case_file), "--output", str(output)],
        cwd=ROOT, check=True, capture_output=True, text=True,
    )
    values = np.loadtxt(output, delimiter=",", skiprows=1)
    columns = output.read_text(encoding="utf-8").splitlines()[0].split(",")
    metadata = json.loads(output.with_suffix(".json").read_text(encoding="utf-8"))
    assert np.isfinite(values).all()
    assert metadata["csv_columns"] == columns
    assert metadata["time_steps"] == len(values)
    assert all(len(item["sha256"]) == 64 for item in metadata["hydro_files"])
    return values, {name: values[:, i] for i, name in enumerate(columns)}, metadata


def _rm3_case(hydro, *, b2b=False):
    return {
        "simulation": {"dt": 0.1, "end_time": 400, "ramp_time": 100},
        "wave": {"type": "regular", "height": 2.5, "period": 8},
        "bodies": [
            {"hydro_file": str(hydro), "hydro_body": 1,
             "mass": "equilibrium", "pitch_inertia": 21_306_090.66},
            {"hydro_file": str(hydro), "hydro_body": 2,
             "mass": "equilibrium", "pitch_inertia": 94_407_091.24},
        ],
        "constraint": {"kind": "floating_joint", "location": [0, 0, 0]},
        "pto": {"kind": "relative_heave", "damping": 1_200_000},
        "body_to_body": b2b,
    }


@pytest.mark.skipif(MODEL != "RM3" or not (CORE and REFERENCE),
                    reason="paired MATLAB RM3 baseline not provided")
def test_rm3_general_runner(tmp_path):
    hydro = Path(CORE) / "examples/RM3/hydroData/rm3.h5"
    _, columns, _ = _run(tmp_path, _rm3_case(hydro))
    for body in (1, 2):
        expected = np.loadtxt(
            Path(REFERENCE) / f"RM3_RM3_body{body}.csv", delimiter=",",
        )
        for axis, dof, limit in (("surge", 0, 0.04), ("heave", 2, 0.005),
                                 ("pitch", 4, 0.001)):
            assert np.max(np.abs(
                columns[f"body{body}_{axis}_position"] - expected[:, 1 + dof]
            )) < limit
    expected_pto = np.loadtxt(Path(REFERENCE) / "RM3_RM3_pto1.csv", delimiter=",")
    assert np.max(np.abs(
        columns["pto_relative_heave_force"] - expected_pto[:, 15]
    )) < 5_000


@pytest.mark.skipif(MODEL != "RM3_B2B" or not (APPLICATIONS and REFERENCE),
                    reason="paired MATLAB RM3 body-to-body baseline not provided")
@pytest.mark.parametrize("case,b2b", [("B2B_Case1", False), ("B2B_Case2", True)])
def test_rm3_body_to_body_general_runner(tmp_path, case, b2b):
    hydro = Path(APPLICATIONS) / "_Common_Input_Files/RM3/hydroData/rm3.h5"
    _, columns, _ = _run(tmp_path, _rm3_case(hydro, b2b=b2b))
    for body in (1, 2):
        expected = np.loadtxt(
            Path(REFERENCE) / f"RM3_B2B_{case}_body{body}.csv", delimiter=",",
        )
        for axis, dof, limit in (("surge", 0, 0.04), ("heave", 2, 0.005),
                                 ("pitch", 4, 0.0012)):
            assert np.max(np.abs(
                columns[f"body{body}_{axis}_position"] - expected[:, 1 + dof]
            )) < limit
    expected_pto = np.loadtxt(
        Path(REFERENCE) / f"RM3_B2B_{case}_pto1.csv", delimiter=",",
    )
    assert np.max(np.abs(
        columns["pto_relative_heave_force"] - expected_pto[:, 15]
    )) < 6_000


@pytest.mark.skipif(MODEL != "OSWEC" or not (CORE and REFERENCE),
                    reason="paired MATLAB OSWEC baseline not provided")
def test_oswec_general_runner(tmp_path):
    reference = Path(REFERENCE)
    components = np.loadtxt(reference / "OSWEC_wave_components.csv", delimiter=",")
    directions = np.loadtxt(reference / "OSWEC_wave_directions.csv", delimiter=",")
    np.savetxt(tmp_path / "phase.csv", components[:, 3:], delimiter=",")
    hydro = Path(CORE) / "examples/OSWEC/hydroData/oswec.h5"
    case = {
        "simulation": {"dt": 0.1, "end_time": 400, "ramp_time": 100,
                       "radiation_memory": 30},
        "wave": {"type": "pm", "height": 2.5, "period": 8,
                 "directions": directions[:, 0].tolist(),
                 "spreading": directions[:, 1].tolist(),
                 "phase_file": "phase.csv"},
        "bodies": [{"hydro_file": str(hydro), "mass": 127_000,
                    "pitch_inertia": 1.85e6}],
        "constraint": {"kind": "fixed_hinge", "location": [0, 0, -8.9]},
        "pto": {"kind": "pitch", "damping": 12_000},
    }
    _, columns, metadata = _run(tmp_path, case)
    assert len(metadata["auxiliary_files"]) == 1
    assert len(metadata["auxiliary_files"][0]["sha256"]) == 64
    expected = np.loadtxt(reference / "OSWEC_OSWEC_body1.csv", delimiter=",")
    wave = np.loadtxt(reference / "OSWEC_wave_elevation.csv", delimiter=",")
    pto = np.loadtxt(reference / "OSWEC_OSWEC_pto1.csv", delimiter=",")
    np.testing.assert_allclose(columns["time"], expected[:, 0], rtol=0, atol=1e-10)
    assert np.max(np.abs(columns["wave_elevation"] - wave[:, 1])) < 1e-9
    assert np.max(np.abs(columns["body1_pitch_position"] - expected[:, 5])) < 0.004
    assert np.max(np.abs(columns["body1_pitch_velocity"] - expected[:, 11])) < 0.005
    assert np.max(np.abs(columns["pto_pitch_torque"] - pto[:, 17])) < 60


@pytest.mark.skipif(MODEL != "Sphere" or not (SPHERE_H5 and REFERENCE),
                    reason="paired MATLAB Sphere baseline not provided")
@pytest.mark.parametrize("case,displacement", [
    ("0m", 0), ("1m", 1), ("1m-ME", 1), ("3m", 3), ("5m", 5),
])
def test_sphere_general_runner(tmp_path, case, displacement):
    config = {
        "simulation": {"dt": 0.01, "end_time": 40,
                       "radiation_memory": 15},
        "wave": {"type": "none"},
        "bodies": [{"hydro_file": SPHERE_H5, "mass": "equilibrium"}],
        "constraint": {"kind": "heave", "initial_displacement": displacement},
    }
    _, columns, _ = _run(tmp_path, config)
    expected = np.loadtxt(Path(REFERENCE) / f"Sphere_{case}_body1.csv", delimiter=",")
    scale = max(1, displacement)
    assert np.max(np.abs(columns["body1_heave_position"] - expected[:, 3])) < 2e-4 * scale
    assert np.max(np.abs(columns["body1_heave_velocity"] - expected[:, 9])) < 2e-4 * scale
    assert np.max(np.abs(columns["body1_total_heave_force"] - expected[:, 15])) < 160 * scale
