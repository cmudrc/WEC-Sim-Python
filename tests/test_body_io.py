"""Read the RM3 hydrodynamic file included with the original port."""

from pathlib import Path
import subprocess
import sys

import numpy as np
import scipy.io

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "source" / "objects"))

from bodyClass import BodyClass  # noqa: E402


def test_rm3_hydrodynamic_bodies_load():
    h5_file = ROOT / "source/objects/rm3.h5"
    for number, expected_name in ((1, "float"), (2, "spar")):
        body = BodyClass(str(h5_file))
        body.bodyNumber = number
        body.readH5file()
        assert body.name == expected_name
        assert body.hydroData["properties"]["name"] == expected_name
        assert body.hydroData["simulation_parameters"]["water_depth"] == "infinite"
        assert body.hydroData["hydro_coeffs"]["excitation"]["re"].shape[0] == 6


def test_rm3_regular_cic_force_preprocessing_matches_matlab_fixture():
    fixture = ROOT / "tests/test_objects/test_bodyclass/testData"
    case = fixture / "body_1_test"
    body = BodyClass(str(fixture / "hydroData/rm3.h5"))
    body.bodyNumber = 1
    body.readH5file()
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {
        "Drag": np.zeros((6, 6)),
        "cd": np.zeros(6),
        "characteristicArea": np.zeros(6),
    }
    body.linearDamping = np.zeros((6, 6))
    body.mass = "equilibrium"
    body.hydroForcePre(
        0.785398163397448, [0], 601,
        np.loadtxt(case / "CTTime.txt").T, [], 0.1, 1000, 9.81,
        "regularCIC", np.loadtxt(case / "waveAmpTime.txt").T,
        1, [], 0, 0, 0,
    )

    for field, filename in (
        ("linearHydroRestCoef", "linearHydroRestCoef.txt"),
        ("fAddedMass", "fAddedMass.txt"),
    ):
        np.testing.assert_allclose(
            body.hydroForce[field], np.loadtxt(case / filename), rtol=1e-12,
        )
    np.testing.assert_allclose(
        body.hydroForce["fExt"]["re"],
        [18508.7581918053, 0.223532780361825, 1444829.62495620,
         6.41939080869329, 140467.443724649, 0.0412557412840770],
        rtol=1e-12,
    )
    np.testing.assert_allclose(
        body.hydroForce["fExt"]["im"],
        [546623.415011964, -0.0428753721889769, 447276.984830190,
         0.273118282507816, 4146814.04776144, 0.0138224113420678],
        rtol=1e-12,
    )
    np.testing.assert_allclose(
        body.hydroForce["irkb"], scipy.io.loadmat(case / "irkb.mat")["a"],
        rtol=1e-12, atol=1e-9,
    )


def test_oswec_directional_irregular_preprocessing_matches_matlab_fixture():
    fixture = ROOT / "tests/test_objects/test_bodyclass/testData"
    case = fixture / "body_3_test"

    def mat(name):
        values = scipy.io.loadmat(case / f"{name}.mat")
        return next(value for key, value in values.items() if not key.startswith("__"))

    body = BodyClass(str(fixture / "hydroData/oswec.h5"))
    body.bodyNumber = 1
    body.readH5file()
    body.hydroStiffness = np.zeros((6, 6))
    body.viscDrag = {
        "Drag": np.zeros((6, 6)),
        "cd": np.zeros(6),
        "characteristicArea": np.zeros(6),
    }
    body.linearDamping = np.zeros((6, 6))
    body.mass = 127000
    body.hydroForcePre(
        mat("w").T, [0, 30, 90], 301, mat("CTTime")[0], 500,
        0.1, 1000, 9.81, "irregular", mat("waveAmpTime").T,
        1, [], 0, 0, 0,
    )

    for field, expected in (
        ("linearHydroRestCoef", mat("linearHydroRestCoef")),
        ("fAddedMass", mat("fAddedMass")),
        ("irkb", mat("irkb")),
    ):
        np.testing.assert_allclose(body.hydroForce[field], expected, rtol=1e-10, atol=1e-8)
    for component in ("re", "im", "md"):
        np.testing.assert_allclose(
            body.hydroForce["fExt"][component], mat(component), rtol=1e-10, atol=1e-8,
        )


def test_rm3_entrypoint_completes_force_preprocessing():
    code = """
import runpy
import numpy as np

state = runpy.run_path('wecSimPython.py')
assert state['simu'].numWecBodies == 2
assert len(state['body']) == 2
for body in state['body']:
    for field in ('linearHydroRestCoef', 'fAddedMass', 'fDamping'):
        matrix = np.asarray(body.hydroForce[field])
        assert matrix.shape == (6, 6)
        assert np.isfinite(matrix).all()
"""
    subprocess.run(
        [sys.executable, "-c", code],
        cwd=ROOT / "source/objects", check=True, capture_output=True, text=True,
    )
