"""Pair the three published imported-spectrum RM3 MCR cases with MATLAB."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import mcr_spectrum_files, run_rm3_spectrum_mcr
from wecsim.irregularWave import imported_spectrum_components, synthesize_irregular_response


APPLICATIONS = os.environ.get("WEC_SIM_APPLICATIONS_DIR")
REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not (APPLICATIONS and REFERENCE),
    reason="paired MATLAB RM3 imported-spectrum output not provided",
)


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: max error {error:.6g} exceeds {limit}"


def test_published_three_sea_state_mcr_against_matlab():
    root = Path(APPLICATIONS)
    reference = Path(REFERENCE)
    mat_file = root / "Multiple_Condition_Runs/RM3_MCROPT3_SeaState/mcrExample.mat"
    files = mcr_spectrum_files(mat_file)
    assert tuple(path.name for path in files) == tuple(
        f"spectrumData{index}.mat" for index in (1, 2, 3)
    )
    hydro = root / "_Common_Input_Files/RM3/hydroData/rm3.h5"
    result = run_rm3_spectrum_mcr(hydro, mat_file)
    assert result.spectrum_files == files
    assert len(result.traces) == 3
    summary = np.loadtxt(reference / "RM3_MCR_SEASTATE_summary.csv", delimiter=",")
    assert summary.shape == (3,)

    for index, spectrum_file in enumerate(files, start=1):
        components = imported_spectrum_components(hydro, spectrum_file)
        expected_components = np.loadtxt(
            reference / f"RM3_MCR_SEASTATE_case{index}_components.csv",
            delimiter=",",
        )
        assert expected_components.shape == (len(components.omega), 4)
        actual_components = np.column_stack((
            components.omega, components.spectral_amplitude,
            components.d_omega, components.phase[:, 0],
        ))
        np.testing.assert_allclose(actual_components, expected_components,
                                   rtol=0, atol=1e-11)
        assert len(components.directions) == 1
        assert components.directions[0] == 0

        trace = result.traces[index - 1]
        response = trace.response
        assert response is not None
        assert response.body_position.shape == (4001, 2, 6)
        expected_wave = np.loadtxt(
            reference / f"RM3_MCR_SEASTATE_case{index}_wave.csv",
            delimiter=",",
        )
        assert expected_wave.shape == (4001, 2)
        np.testing.assert_allclose(trace.time, expected_wave[:, 0],
                                   rtol=0, atol=1e-10)
        _max_error(result.wave_elevation[index - 1], expected_wave[:, 1],
                   1e-10, f"case {index} wave elevation")

        for body in (1, 2):
            expected = np.loadtxt(
                reference / f"RM3_MCR_SEASTATE_case{index}_body{body}.csv",
                delimiter=",",
            )
            assert expected.shape == (4001, 25)
            np.testing.assert_allclose(response.time, expected[:, 0],
                                       rtol=0, atol=1e-10)
            incident = synthesize_irregular_response(
                hydro, components, dt=0.1, end_time=400,
                ramp_time=100, body_number=body,
            )
            for dof, name, force_limit in (
                (0, "surge", 1), (2, "heave", 1), (4, "pitch", 1),
            ):
                _max_error(incident.excitation_force[:, dof],
                           expected[:, 19 + dof], force_limit,
                           f"case {index} body {body} {name} excitation")
            for dof, name, position_limit, velocity_limit in (
                (0, "surge", 0.12, 0.06),
                (2, "heave", 0.012, 0.008),
                (4, "pitch", 0.0045, 0.002),
            ):
                _max_error(response.body_position[:, body - 1, dof],
                           expected[:, 1 + dof], position_limit,
                           f"case {index} body {body} {name} position")
                _max_error(response.body_velocity[:, body - 1, dof],
                           expected[:, 7 + dof], velocity_limit,
                           f"case {index} body {body} {name} velocity")

        pto = np.loadtxt(
            reference / f"RM3_MCR_SEASTATE_case{index}_pto1.csv",
            delimiter=",",
        )
        assert pto.shape == (4001, 25)
        np.testing.assert_allclose(response.time, pto[:, 0],
                                   rtol=0, atol=1e-10)
        center_gap = (response.body_position[0, 0, 2]
                      - response.body_position[0, 1, 2])
        stroke = (response.body_position[:, 0, 2]
                  - response.body_position[:, 1, 2]
                  - center_gap * np.cos(response.body_position[:, 0, 4]))
        _max_error(stroke, pto[:, 3], 0.015, f"case {index} PTO stroke")
        _max_error(response.pto_velocity, pto[:, 9], 0.009,
                   f"case {index} PTO speed")
        _max_error(response.pto_force, pto[:, 15], 11_000,
                   f"case {index} PTO force")
        _max_error(trace.absorbed_power, -pto[:, 21], 10_000,
                   f"case {index} PTO absorbed power")
        np.testing.assert_allclose(-summary[index - 1],
                                   np.mean(-pto[1999:, 21]), rtol=0, atol=1e-6)

    _max_error(result.mean_absorbed_power, -summary, 3_000,
               "three sea-state MCR mean powers")
