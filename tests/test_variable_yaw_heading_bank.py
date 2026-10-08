"""Nearest-heading selection used by the published variable-hydro yaw case."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularWave, WEC, WorldPoint
from wecsim.passiveYaw import NearestHeadingExcitation, PassiveYawExcitation


def _model():
    headings = np.arange(0.0, 360.0, 10.0)
    real = np.zeros((6, len(headings)))
    real[5] = headings
    return PassiveYawExcitation(
        headings, real, np.zeros_like(real),
        incident_direction=10, omega=1, amplitude=1, ramp_time=0,
    )


def test_nearest_bank_uses_source_relative_angle_and_first_tie():
    bank = NearestHeadingExcitation(_model(), np.arange(-40, 41, 2))
    assert bank.heading(0) == 10
    assert bank.heading(np.deg2rad(9)) == 0
    assert bank.heading(np.deg2rad(-31)) == 40
    assert bank.heading(np.deg2rad(51)) == -40
    assert bank.heading(np.deg2rad(9)) == 0  # a one-degree tie chooses lower

    yaw = np.deg2rad(9.2)
    actual = bank.force(0, yaw)
    expected = _model().force(0, yaw, coefficient_heading=0)
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-13)
    assert actual[5] == 0


def test_nearest_bank_rejects_invalid_grid():
    model = _model()
    for headings in ([0], [0, 0], [2, 1], [-181, 0], [0, np.nan]):
        with pytest.raises(ValueError, match="heading bank"):
            NearestHeadingExcitation(model, headings)


@pytest.mark.skipif(
    not (os.environ.get("WEC_SIM_APPLICATIONS_DIR")
         and os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")),
    reason="fresh variable-yaw MATLAB output not provided",
)
def test_two_degree_bank_against_pinned_matlab():
    applications = Path(os.environ["WEC_SIM_APPLICATIONS_DIR"])
    source = Path(os.environ["WEC_SIM_MATLAB_MODEL_OUTPUT_DIR"])
    prefix = "OSWEC_VARIABLE_YAW_2DEG"
    flap = np.loadtxt(source / f"{prefix}_regular_2deg_body1.csv", delimiter=",")
    base = np.loadtxt(source / f"{prefix}_regular_2deg_body2.csv", delimiter=",")
    pto = np.loadtxt(source / f"{prefix}_regular_2deg_pto1.csv", delimiter=",")
    selected = np.loadtxt(source / f"{prefix}_selected_heading.csv", delimiter=",")
    wave = np.loadtxt(source / f"{prefix}_wave.csv", delimiter=",")
    hydro = applications / "_Common_Input_Files/OSWEC/hydroData/oswec.h5"

    wec = WEC("OSWEC two-degree variable yaw")
    moving = wec.body(
        "flap", hydro, mass=12700, inertia=(1.85e6,) * 3,
        passive_yaw=True, yaw_heading_bank=np.arange(-40, 41, 2),
    )
    wec.body("base", hydro, mass=999, inertia=(999,) * 3)
    yaw = wec.coordinate("yaw", moving.move("yaw", pivot=WorldPoint(0, 0, -8.9)))
    wec.rotational_pto("hinge", yaw, damping=120000)
    result = wec.run(RegularWave(2.5, 8, 10), dt=0.01,
                     end_time=600, ramp_time=100)

    assert flap.shape == base.shape == pto.shape == (60001, 25)
    assert selected.shape == (60001, 3)
    assert wave.shape == (60001, 2)
    np.testing.assert_allclose(result.time, flap[:, 0], rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.wave_elevation, wave[:, 1], rtol=0, atol=1e-12)
    np.testing.assert_allclose(result.bodies["base"].position,
                               base[:, 1:7], rtol=0, atol=1e-10)
    np.testing.assert_allclose(result.bodies["base"].velocity,
                               base[:, 7:13], rtol=0, atol=1e-10)

    model = PassiveYawExcitation.from_hydro_data(
        _hydro_data(hydro), omega=2 * np.pi / 8,
        incident_direction=10, amplitude=1.25, ramp_time=100,
        rho=1000, g=9.81,
    )
    bank = NearestHeadingExcitation(model, np.arange(-40, 41, 2))
    source_heading = np.array([bank.heading(angle) for angle in flap[:, 6]])
    heading_error = np.max(np.abs(source_heading - selected[:, 2]))
    print(f"source heading selection max error: {heading_error:.6g} deg")
    assert heading_error == 0

    source_path_force = np.stack([
        bank.force(t, angle) for t, angle in zip(result.time, flap[:, 6])
    ])
    force_error = np.max(np.abs(source_path_force[:, 5] - flap[:, 24]))
    yaw_error = np.max(np.abs(result.bodies["flap"].position[:, 5] - flap[:, 6]))
    speed_error = np.max(np.abs(result.bodies["flap"].velocity[:, 5] - flap[:, 12]))
    torque_error = np.max(np.abs(result.ptos["hinge"].force - pto[:, 17]))
    print(f"source-path yaw excitation max error: {force_error:.6g} N m")
    print(f"trajectory max errors: {yaw_error:.6g} rad, "
          f"{speed_error:.6g} rad/s, {torque_error:.6g} N m")
    assert np.isfinite([force_error, yaw_error, speed_error, torque_error]).all()
    assert force_error < 1000

    # Before the first heading event separates the two paths, the motion
    # remains paired. The source's discontinuous selector magnifies the
    # subsequent one-sample event offset; full-run trajectory parity is open.
    before_split = result.time < 64.57
    assert np.max(np.abs(result.bodies["flap"].position[before_split, 5]
                         - flap[before_split, 6])) < 5e-5
    assert np.max(np.abs(result.bodies["flap"].velocity[before_split, 5]
                         - flap[before_split, 12])) < 5e-5
    assert np.max(np.abs(result.ptos["hinge"].force[before_split]
                         - pto[before_split, 17])) < 5


def _hydro_data(path):
    from wecsim.bodyClass import BodyClass

    body = BodyClass(str(path))
    body.bodyNumber = 1
    body.bodyTotal = 2
    body.readH5file()
    return body.hydroData
