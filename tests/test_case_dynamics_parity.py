"""Run JSON dynamics cases through the public CLI against paired MATLAB data."""

import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from wecsim import NoWave, PMWave, RegularWave, WEC, WorldPoint

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
    completed = subprocess.run(
        [sys.executable, "-m", "wecsim",
         str(case_file), "--output", str(output)],
        cwd=ROOT, check=False, capture_output=True, text=True,
    )
    assert completed.returncode == 0, completed.stderr
    values = np.loadtxt(output, delimiter=",", skiprows=1)
    columns = output.read_text(encoding="utf-8").splitlines()[0].split(",")
    metadata = json.loads(output.with_suffix(".json").read_text(encoding="utf-8"))
    assert np.isfinite(values).all()
    assert metadata["csv_columns"] == columns
    assert metadata["time_steps"] == len(values)
    assert all(len(item["sha256"]) == 64 for item in metadata["hydro_files"])
    return values, {name: values[:, i] for i, name in enumerate(columns)}, metadata


def _max_error(actual, expected, limit, signal):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, f"{signal}: mismatched time series"
    assert np.isfinite(expected).all(), f"{signal}: nonfinite MATLAB reference"
    error = np.max(np.abs(actual - expected))
    assert error < limit, f"{signal}: max error {error:.6g} exceeds {limit:.6g}"


def _assert_inactive_dofs(columns, expected, body, inactive):
    for axis, dof in inactive:
        prefix = f"body{body}_{axis}"
        _max_error(columns[f"{prefix}_position"], expected[:, 1 + dof],
                   1e-10, f"{prefix} stationary position")
        _max_error(columns[f"{prefix}_velocity"], expected[:, 7 + dof],
                   1e-10, f"{prefix} stationary velocity")


def _assert_rm3_body(columns, expected, body, *, pitch_limit):
    np.testing.assert_allclose(columns["time"], expected[:, 0], rtol=0, atol=1e-10)
    for axis, dof, position_limit, velocity_limit in (
        ("surge", 0, 0.04, 0.035),
        ("heave", 2, 0.005, 0.004),
        ("pitch", 4, pitch_limit, 0.001),
    ):
        prefix = f"body{body}_{axis}"
        _max_error(columns[f"{prefix}_position"], expected[:, 1 + dof],
                   position_limit, f"{prefix} position")
        _max_error(columns[f"{prefix}_velocity"], expected[:, 7 + dof],
                   velocity_limit, f"{prefix} velocity")
    _assert_inactive_dofs(columns, expected, body,
                          (("sway", 1), ("roll", 3), ("yaw", 5)))


def _assert_rm3_pto(columns, expected, *, force_limit):
    np.testing.assert_allclose(columns["time"], expected[:, 0], rtol=0, atol=1e-10)
    # Undo the two body-center rotation terms to recover joint-relative PTO
    # motion. The HDF5 centers and the floating joint lie on the z axis.
    center_gap = (columns["body1_heave_position"][0]
                  - columns["body2_heave_position"][0])
    pitch = columns["body1_pitch_position"]
    pitch_speed = columns["body1_pitch_velocity"]
    cosine = np.cos(pitch)
    stroke = ((columns["body1_heave_position"]
               - columns["body2_heave_position"]) / cosine - center_gap)
    speed = ((columns["body1_heave_velocity"]
              - columns["body2_heave_velocity"]
              + (center_gap + stroke) * np.sin(pitch) * pitch_speed) / cosine)
    force = columns["pto_relative_heave_force"]
    _max_error(stroke, expected[:, 3], 0.009, "RM3 PTO stroke")
    _max_error(speed, expected[:, 9], 0.0055, "RM3 PTO speed")
    _max_error(force, expected[:, 15], force_limit, "RM3 PTO force")
    _max_error(force * speed, expected[:, 21], 7_500, "RM3 PTO power")


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
    hydro = (Path(CORE) / "examples/RM3/hydroData/rm3.h5").resolve()
    _, columns, _ = _run(tmp_path, _rm3_case(hydro))
    for body in (1, 2):
        expected = np.loadtxt(
            Path(REFERENCE) / f"RM3_RM3_body{body}.csv", delimiter=",",
        )
        _assert_rm3_body(columns, expected, body, pitch_limit=0.001)
    expected_pto = np.loadtxt(Path(REFERENCE) / "RM3_RM3_pto1.csv", delimiter=",")
    _assert_rm3_pto(columns, expected_pto, force_limit=5_000)


@pytest.mark.skipif(MODEL != "RM3" or not (CORE and REFERENCE),
                    reason="paired MATLAB RM3 heave baseline not provided")
def test_rm3_linear_subspace_general_runner(tmp_path):
    hydro = (Path(CORE) / "examples/RM3/hydroData/rm3.h5").resolve()
    first = [[0, 0], [0, 0], [1, 0], [0, 0], [0, 0], [0, 0]]
    second = [[0, 0], [0, 0], [0, 1], [0, 0], [0, 0], [0, 0]]
    case = {
        "simulation": {"dt": 0.1, "end_time": 400, "ramp_time": 100},
        "wave": {"type": "regular", "height": 2.5, "period": 8},
        "bodies": [
            {"hydro_file": str(hydro), "coordinate_map": first},
            {"hydro_file": str(hydro), "coordinate_map": second},
        ],
        "constraint": {"kind": "linear_subspace"},
        "pto": {"kind": "linear", "damping_matrix": [
            [1_200_000, -1_200_000], [-1_200_000, 1_200_000],
        ]},
    }
    _, columns, _ = _run(tmp_path, case)
    for body in (1, 2):
        expected = np.loadtxt(
            Path(REFERENCE) / f"RM3_RM3_body{body}.csv", delimiter=",",
        )
        np.testing.assert_allclose(columns["time"], expected[:, 0], rtol=0, atol=1e-10)
        assert np.max(np.abs(
            columns[f"body{body}_heave_position"] - expected[:, 3]
        )) < 0.0065
        assert np.max(np.abs(
            columns[f"body{body}_heave_velocity"] - expected[:, 9]
        )) < 0.0065
    expected_pto = np.loadtxt(Path(REFERENCE) / "RM3_RM3_pto1.csv", delimiter=",")
    assert np.max(np.abs(columns["pto_coordinate1_force"] - expected_pto[:, 15])) < 9_000
    stroke = ((columns["body1_heave_position"] - columns["body1_heave_position"][0])
              - (columns["body2_heave_position"] - columns["body2_heave_position"][0]))
    speed = columns["body1_heave_velocity"] - columns["body2_heave_velocity"]
    _max_error(stroke, expected_pto[:, 3], 0.01, "linear RM3 PTO stroke")
    _max_error(speed, expected_pto[:, 9], 0.008, "linear RM3 PTO speed")
    _max_error(columns["pto_coordinate1_force"] * speed, expected_pto[:, 21],
               12_000, "linear RM3 PTO power")
    np.testing.assert_allclose(
        columns["pto_coordinate1_force"], -columns["pto_coordinate2_force"],
        rtol=0, atol=1e-9,
    )
    wec = WEC("RM3 heave in Python")
    float_body = wec.body("float", hydro)
    spar = wec.body("spar", hydro)
    wec.coordinate("float_heave", float_body.move("heave"))
    wec.coordinate("spar_heave", spar.move("heave"))
    wec.pto("main", float_body.at(0, 0, 0), spar.at(0, 0, 0),
            axis=(0, 0, 1), damping=1_200_000)
    python = wec.run(RegularWave(2.5, 8), dt=0.1, end_time=400, ramp_time=100)
    for body, name in ((1, "float"), (2, "spar")):
        expected = np.loadtxt(
            Path(REFERENCE) / f"RM3_RM3_body{body}.csv", delimiter=",",
        )
        assert np.max(np.abs(
            python.bodies[name].position[:, 2] - expected[:, 3]
        )) < 0.0065
    assert np.max(np.abs(
        -python.ptos["main"].force - expected_pto[:, 15]
    )) < 9_000
    _max_error(-python.ptos["main"].stroke, expected_pto[:, 3],
               0.01, "Python API RM3 PTO stroke")
    _max_error(-python.ptos["main"].velocity, expected_pto[:, 9],
               0.008, "Python API RM3 PTO speed")
    _max_error(-python.ptos["main"].absorbed_power, expected_pto[:, 21],
               12_000, "Python API RM3 PTO power")


@pytest.mark.skipif(MODEL != "RM3_B2B" or not (APPLICATIONS and REFERENCE),
                    reason="paired MATLAB RM3 body-to-body baseline not provided")
@pytest.mark.parametrize("case,b2b", [("B2B_Case1", False), ("B2B_Case2", True)])
def test_rm3_body_to_body_general_runner(tmp_path, case, b2b):
    hydro = (Path(APPLICATIONS) / "_Common_Input_Files/RM3/hydroData/rm3.h5").resolve()
    _, columns, _ = _run(tmp_path, _rm3_case(hydro, b2b=b2b))
    for body in (1, 2):
        expected = np.loadtxt(
            Path(REFERENCE) / f"RM3_B2B_{case}_body{body}.csv", delimiter=",",
        )
        _assert_rm3_body(columns, expected, body, pitch_limit=0.0012)
    expected_pto = np.loadtxt(
        Path(REFERENCE) / f"RM3_B2B_{case}_pto1.csv", delimiter=",",
    )
    _assert_rm3_pto(columns, expected_pto, force_limit=6_000)


@pytest.mark.skipif(MODEL != "OSWEC" or not (CORE and REFERENCE),
                    reason="paired MATLAB OSWEC baseline not provided")
def test_oswec_general_runner(tmp_path):
    reference = Path(REFERENCE)
    components = np.loadtxt(reference / "OSWEC_wave_components.csv", delimiter=",")
    directions = np.loadtxt(reference / "OSWEC_wave_directions.csv", delimiter=",")
    np.savetxt(tmp_path / "phase.csv", components[:, 3:], delimiter=",")
    hydro = (Path(CORE) / "examples/OSWEC/hydroData/oswec.h5").resolve()
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
    np.testing.assert_allclose(columns["time"], wave[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(columns["time"], pto[:, 0], rtol=0, atol=1e-10)
    _max_error(columns["wave_elevation"], wave[:, 1], 1e-11, "OSWEC wave elevation")
    for axis, dof, position_limit, velocity_limit in (
        ("surge", 0, 0.02, 0.03),
        ("heave", 2, 0.005, 0.005),
        ("pitch", 4, 0.004, 0.005),
    ):
        prefix = f"body1_{axis}"
        _max_error(columns[f"{prefix}_position"], expected[:, 1 + dof],
                   position_limit, f"OSWEC {axis} position")
        _max_error(columns[f"{prefix}_velocity"], expected[:, 7 + dof],
                   velocity_limit, f"OSWEC {axis} velocity")
    _assert_inactive_dofs(columns, expected, 1,
                          (("sway", 1), ("roll", 3), ("yaw", 5)))
    _max_error(columns["body1_pitch_position"], pto[:, 5], 0.004,
               "OSWEC PTO angle")
    _max_error(columns["body1_pitch_velocity"], pto[:, 11], 0.005,
               "OSWEC PTO angular speed")
    torque = columns["pto_pitch_torque"]
    _max_error(torque, pto[:, 17], 60, "OSWEC PTO torque")
    _max_error(torque * columns["body1_pitch_velocity"], pto[:, 23],
               max(25, 0.06 * np.max(np.abs(pto[:, 23]))), "OSWEC PTO power")

    configured = WEC("Published OSWEC")
    flap = configured.body("flap", hydro, mass=127_000,
                           inertia=(1.85e6,) * 3)
    base = configured.body("base", hydro, mass=999,
                           inertia=(999,) * 3)
    configured.fixed_hinge(
        flap, base, location=WorldPoint(0, 0, -10),
        pto_location=WorldPoint(0, 0, -8.9), damping=12_000,
    )
    public = configured.run(
        PMWave(2.5, 8, directions=tuple(directions[:, 0]),
               spreading=tuple(directions[:, 1]), phase_file="phase.csv"),
        dt=0.1, end_time=400, ramp_time=100,
        radiation_memory=30, base_dir=tmp_path,
    )
    for dof, label, position_limit, velocity_limit in (
        (0, "surge", 0.02, 0.03),
        (2, "heave", 0.005, 0.005),
        (4, "pitch", 0.004, 0.005),
    ):
        _max_error(public.bodies["flap"].position[:, dof],
                   expected[:, 1 + dof], position_limit,
                   f"configured OSWEC {label} position")
        _max_error(public.bodies["flap"].velocity[:, dof],
                   expected[:, 7 + dof], velocity_limit,
                   f"configured OSWEC {label} velocity")
    _max_error(public.ptos["hinge"].force, pto[:, 17],
               60, "configured OSWEC PTO torque")
    _max_error(-public.ptos["hinge"].absorbed_power, pto[:, 23],
               max(25, 0.06 * np.max(np.abs(pto[:, 23]))),
               "configured OSWEC PTO power")
    expected_base = np.loadtxt(
        reference / "OSWEC_OSWEC_body2.csv", delimiter=",",
    )
    np.testing.assert_allclose(public.bodies["base"].position,
                               expected_base[:, 1:7], rtol=0, atol=1e-10)
    np.testing.assert_allclose(public.bodies["base"].velocity,
                               expected_base[:, 7:13], rtol=0, atol=1e-10)
    np.testing.assert_allclose(public.coordinates["pitch"].position,
                               public.ptos["hinge"].stroke, rtol=0, atol=0)


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
        "bodies": [{"hydro_file": str(Path(SPHERE_H5).resolve()),
                    "mass": "equilibrium"}],
        "constraint": {"kind": "heave", "initial_displacement": displacement},
    }
    _, columns, _ = _run(tmp_path, config)
    expected = np.loadtxt(Path(REFERENCE) / f"Sphere_{case}_body1.csv", delimiter=",")
    np.testing.assert_allclose(columns["time"], expected[:, 0], rtol=0, atol=1e-10)
    _assert_inactive_dofs(columns, expected, 1,
                          (("surge", 0), ("sway", 1), ("roll", 3),
                           ("pitch", 4), ("yaw", 5)))
    if displacement == 0:
        _max_error(columns["body1_heave_position"], expected[:, 3],
                   1e-10, "stationary Sphere center position")
        _max_error(expected[:, 3], np.full_like(expected[:, 3], expected[0, 3]),
                   1e-10, "MATLAB stationary Sphere center")
        for name, dof in (("body1_heave_velocity", 9),
                          ("body1_total_heave_force", 15)):
            _max_error(expected[:, dof], np.zeros_like(expected[:, dof]),
                       1e-10, f"MATLAB stationary {name}")
            _max_error(columns[name], expected[:, dof], 1e-10, name)
        return
    scale = max(1, displacement)
    _max_error(columns["body1_heave_position"], expected[:, 3],
               1.2e-4 * scale, "Sphere heave position")
    _max_error(columns["body1_heave_velocity"], expected[:, 9],
               1.6e-4 * scale, "Sphere heave velocity")
    _max_error(columns["body1_total_heave_force"], expected[:, 15],
               100 * scale, "Sphere total heave force")


@pytest.mark.skipif(MODEL != "Sphere" or not (SPHERE_H5 and REFERENCE),
                    reason="paired MATLAB Sphere baseline not provided")
def test_sphere_linear_subspace_general_runner(tmp_path):
    case = {
        "simulation": {"dt": 0.01, "end_time": 40,
                       "radiation_memory": 15},
        "wave": {"type": "none"},
        "bodies": [{"hydro_file": str(Path(SPHERE_H5).resolve()),
                    "coordinate_map": [[0], [0], [1], [0], [0], [0]]}],
        "constraint": {"kind": "linear_subspace", "initial_coordinate": [5]},
    }
    _, columns, _ = _run(tmp_path, case)
    expected = np.loadtxt(Path(REFERENCE) / "Sphere_5m_body1.csv", delimiter=",")
    assert not any(name.startswith("pto_") for name in columns)
    assert np.max(np.abs(columns["body1_heave_position"] - expected[:, 3])) < 0.001
    assert np.max(np.abs(columns["body1_heave_velocity"] - expected[:, 9])) < 0.001
    wec = WEC("Sphere free decay in Python")
    sphere = wec.body("sphere", Path(SPHERE_H5).resolve())
    wec.coordinate("heave", sphere.move("heave"))
    python = wec.run(
        NoWave(), dt=0.01, end_time=40, radiation_memory=15,
        initial_coordinate={"heave": 5},
    )
    assert np.max(np.abs(python.bodies["sphere"].position[:, 2]
                         - expected[:, 3])) < 0.001
    assert np.max(np.abs(python.bodies["sphere"].velocity[:, 2]
                         - expected[:, 9])) < 0.001
