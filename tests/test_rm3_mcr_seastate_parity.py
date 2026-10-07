"""Pair the three published imported-spectrum RM3 MCR cases with MATLAB."""

import os
from pathlib import Path

import h5py
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


def test_matlab_applied_added_mass_matches_pinned_coefficients():
    """Pair the matrix used by Simscape with the BEM matrix and mass split."""
    root = Path(APPLICATIONS)
    reference = Path(REFERENCE)
    hydro = root / "_Common_Input_Files/RM3/hydroData/rm3.h5"
    input_inertias = ((20_907_301, 21_306_090.66, 37_085_481.11),
                     (94_419_614.57, 94_407_091.24, 28_542_224.82))
    with h5py.File(hydro) as h5:
        rho = float(np.asarray(h5["simulation_parameters/rho"]).item())
        for body in (1, 2):
            prefix = f"body{body}"
            original = rho * np.asarray(
                h5[f"{prefix}/hydro_coeffs/added_mass/inf_freq"]
            )[:, 6 * (body - 1):6 * body]
            applied = original.copy()
            mass_shift = 2 * np.trace(original[:3, :3])
            applied[:3, :3] -= np.eye(3) * mass_shift
            applied[3:6, 3:6] -= np.diag(np.diag(original[3:6, 3:6]))
            for row, column in ((3, 4), (3, 5), (4, 5)):
                applied[row, column] -= original[row, column]
                applied[column, row] -= original[row, column]
            saved_original = np.loadtxt(reference / (
                f"RM3_MCR_SEASTATE_body{body}_added_mass_original.csv"),
                delimiter=",")
            saved_applied = np.loadtxt(reference / (
                f"RM3_MCR_SEASTATE_body{body}_added_mass_applied.csv"),
                delimiter=",")
            np.testing.assert_allclose(saved_original, original,
                                       rtol=1e-12, atol=1e-6)
            np.testing.assert_allclose(saved_applied, applied,
                                       rtol=1e-12, atol=1e-6)
            volume = float(np.asarray(h5[f"{prefix}/properties/disp_vol"]).item())
            nominal_mass = rho * volume
            properties = np.loadtxt(reference / (
                f"RM3_MCR_SEASTATE_body{body}_mass_properties.csv"),
                delimiter=",")
            expected_properties = np.r_[
                nominal_mass, nominal_mass, nominal_mass + mass_shift,
                input_inertias[body - 1],
                np.asarray(input_inertias[body - 1]) + np.diag(original)[3:6],
            ]
            np.testing.assert_allclose(properties, expected_properties,
                                       rtol=1e-12, atol=1e-6)


def test_matlab_sea_state_joint_and_force_balance():
    """Identify the joint motion and applied mass convention in the baseline.

    MATLAB moves twice the translational added-mass trace into each Simscape
    body. Its exported forceTotal retains the matching modified force, so the
    adjusted mass is needed to close the logged translational balance.
    """
    root = Path(APPLICATIONS)
    reference = Path(REFERENCE)
    hydro = root / "_Common_Input_Files/RM3/hydroData/rm3.h5"
    adjusted_mass = []
    centers = []
    inertias = (21_306_090.66, 94_407_091.24)
    with h5py.File(hydro) as h5:
        rho = float(np.asarray(h5["simulation_parameters/rho"]).item())
        for body in (1, 2):
            prefix = f"body{body}"
            volume = float(np.asarray(h5[f"{prefix}/properties/disp_vol"]).item())
            centers.append(float(np.asarray(h5[f"{prefix}/properties/cg"]).ravel()[2]))
            added = rho * np.asarray(
                h5[f"{prefix}/hydro_coeffs/added_mass/inf_freq"]
            )[:, 6 * (body - 1):6 * body]
            adjusted_mass.append(rho * volume + 2 * np.trace(added[:3, :3]))

    for case in (1, 2, 3):
        bodies = [np.loadtxt(reference / f"RM3_MCR_SEASTATE_case{case}_body{body}.csv",
                             delimiter=",") for body in (1, 2)]
        forces = [np.loadtxt(reference / f"RM3_MCR_SEASTATE_case{case}_body{body}_forces.csv",
                             delimiter=",") for body in (1, 2)]
        pto = np.loadtxt(reference / f"RM3_MCR_SEASTATE_case{case}_pto1.csv",
                         delimiter=",")
        angle = bodies[0][:, 5]
        sine, cosine = np.sin(angle), np.cos(angle)
        separation = centers[0] - centers[1] + pto[:, 3]
        _max_error(bodies[0][:, 1] - bodies[1][:, 1],
                   separation * sine, 1e-9, f"case {case} joint x")
        _max_error(bodies[0][:, 3] - bodies[1][:, 3],
                   separation * cosine, 1e-9, f"case {case} joint z")
        _max_error(pto[:, 15], -1_200_000 * pto[:, 9],
                   1e-6, f"case {case} PTO damper")

        residual = np.zeros((len(angle), 4))
        for index in (0, 1):
            radius = bodies[index][:, 3] / cosine
            imbalance_x = (adjusted_mass[index] * forces[index][:, 19]
                           - bodies[index][:, 13])
            imbalance_z = (adjusted_mass[index] * forces[index][:, 21]
                           - bodies[index][:, 15])
            imbalance_pitch = (inertias[index] * forces[index][:, 23]
                               - bodies[index][:, 17])
            residual[:, 0] += imbalance_x
            residual[:, index + 1] += sine * imbalance_x + cosine * imbalance_z
            residual[:, 3] += (radius * cosine * imbalance_x
                               - radius * sine * imbalance_z
                               + imbalance_pitch)
        residual[:, 1] -= pto[:, 15]
        residual[:, 2] += pto[:, 15]
        _max_error(residual, np.zeros_like(residual), 1e-3,
                   f"case {case} generalized force balance")


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
