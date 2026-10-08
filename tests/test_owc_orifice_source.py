"""Pair the published OWC orifice block with its saved MATLAB signal path."""

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from wecsim import OrificePTO


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="fresh MATLAB OWC orifice output absent",
)


def _max_error(actual, expected, limit, label):
    assert actual.shape == expected.shape, label
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    error = float(np.max(np.abs(actual - expected)))
    assert error < limit, f"{label}: {error:.6g} exceeds {limit}"


def test_published_orifice_force_on_matlab_piston_motion():
    reference = Path(REFERENCE)
    logged = np.loadtxt(reference / "OWC_ORIFICE_orifice.csv", delimiter=",")
    parameters = np.loadtxt(reference / "OWC_ORIFICE_parameters.csv", delimiter=",")
    flex = loadmat(reference / "OWC_ORIFICE_Flex_out.mat",
                   simplify_cells=True)["Flex_out"]
    time = np.asarray(flex["time"])
    values = np.asarray(flex["signals"]["values"])
    assert logged.ndim == 2 and logged.shape[1] == 6
    assert values.ndim == 2 and values.shape[1] == 10
    assert len(time) == len(values) and np.all(np.diff(time) > 0)
    assert np.all(np.diff(logged[:, 0]) > 0)
    model = OrificePTO(*parameters)
    # The orifice block logs every ode23t internal evaluation (418,828 rows),
    # while Flex_out is sampled at 0.005 s (26,001 rows). The logged flow gives
    # the exact internal piston speed. Its interpolated value agrees with the
    # independent flexible-mode speed on the output grid.
    piston_speed = logged[:, 4] / model.piston_area
    sampled_speed = np.interp(time, logged[:, 0], piston_speed)
    _max_error(sampled_speed, values[:, 1], 1e-10,
               "orifice input versus flexible-mode speed")
    result = model.evaluate(piston_speed)

    # MATLAB logs [force N, Mach flag, pressure kPa, flow m^3/s, power kW].
    _max_error(result.force, logged[:, 1], 1e-8, "orifice reaction")
    np.testing.assert_array_equal(result.compressibility_flag, logged[:, 2])
    _max_error(result.pressure_drop / 1000, logged[:, 3], 1e-10,
               "orifice pressure")
    _max_error(result.flow_rate, logged[:, 4], 1e-13, "orifice flow")
    _max_error(result.absorbed_power / 1000, logged[:, 5], 1e-10,
               "orifice power")
    _max_error(result.absorbed_power, -result.force * piston_speed,
               1e-8, "passive orifice energy")
