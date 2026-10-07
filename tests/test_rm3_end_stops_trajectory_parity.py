"""Compare the opt-in implicit stop solver with refined MATLAB source runs."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import LinearHardStops, solve_rm3_regular


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
MODEL = os.environ.get("WEC_SIM_REFERENCE_MODEL")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="refined MATLAB RM3 End_Stops output not provided",
)


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = np.max(np.abs(actual - expected))
    assert error < limit, f"{label}: max error {error:.6g} exceeds {limit}"


def test_refined_end_stop_motion_against_matlab():
    # These variants change only the MATLAB integration step and end time;
    # the published 0.1 s source run is tracked separately for step sensitivity.
    variants = {
        "RM3_END_STOPS_STEP_FINE": (0.025, "End_Stops_dt0025"),
        "RM3_END_STOPS_STEP_FINER": (0.0125, "End_Stops_dt00125"),
    }
    assert MODEL in variants
    dt, case = variants[MODEL]
    reference = Path(REFERENCE)
    hydro = (Path(APPLICATIONS) /
             "_Common_Input_Files/RM3/hydroData/rm3.h5").resolve()
    pto = np.loadtxt(reference / f"{MODEL}_{case}_pto1.csv", delimiter=",")
    assert pto.shape == (round(120 / dt) + 1, 49)
    response = solve_rm3_regular(
        hydro, pto_hard_stops=LinearHardStops(-0.6, 0.6, 1e8, 1e8),
        dt=dt, end_time=120,
    )
    np.testing.assert_allclose(response.time, pto[:, 0], rtol=0, atol=1e-10)
    for body in (1, 2):
        expected = np.loadtxt(
            reference / f"{MODEL}_{case}_body{body}.csv", delimiter=",",
        )
        assert expected.shape == (len(response.time), 25)
        np.testing.assert_allclose(response.time, expected[:, 0], rtol=0, atol=1e-10)
        for dof, name, position_limit, speed_limit in (
            (0, "surge", 0.005, 0.005),
            (2, "heave", 0.006, 0.03),
            (4, "pitch", 1e-4, 1e-4),
        ):
            _max_error(response.body_position[:, body - 1, dof],
                       expected[:, 1 + dof], position_limit,
                       f"body {body} {name} position")
            _max_error(response.body_velocity[:, body - 1, dof],
                       expected[:, 7 + dof], speed_limit,
                       f"body {body} {name} speed")
    _max_error(response.pto_stroke, pto[:, 3], 0.006, "PTO stroke")
    _max_error(response.pto_velocity, pto[:, 9], 0.035, "PTO speed")
    _max_error(response.pto_force, pto[:, 15], 450_000, "total PTO force")
    python_energy = np.trapezoid(response.pto_dissipated_power, response.time)
    matlab_energy = np.trapezoid(1_200_000 * pto[:, 9]**2, pto[:, 0])
    assert abs(python_energy - matlab_energy) / matlab_energy < 0.007
