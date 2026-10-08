"""Pinned RM3 hydraulic PTO source signals and replayed incident sea.

This is an input/force diagnostic. It does not claim coupled hydraulic
trajectory parity with a Python PTO model.
"""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim.hydraulic import CompressibleCylinder
from wecsim.irregularWave import (
    pm_equal_energy_components, synthesize_irregular_response,
)


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
RM3_H5 = os.environ.get("WEC_SIM_RM3_H5")
pytestmark = pytest.mark.skipif(
    not (REFERENCE and RM3_H5),
    reason="pinned RM3 hydraulic MATLAB output and HDF5 absent",
)


def _load(name):
    values = np.loadtxt(Path(REFERENCE) / f"RM3_HYDRAULIC_SOURCE_{name}.csv",
                        delimiter=",")
    assert values.ndim == 2 and np.isfinite(values).all(), name
    return values


def test_hydraulic_source_wave_and_component_traces():
    components = _load("components")
    matlab_wave = _load("wave")
    assert components.shape == (500, 4)
    assert matlab_wave.shape == (40_001, 2)
    generated = pm_equal_energy_components(
        RM3_H5, significant_height=2.5, peak_period=8,
        directions=(0,), spreading=(1,), phase=components[:, 3:4],
    )
    for actual, saved in ((generated.omega, components[:, 0]),
                          (generated.spectral_amplitude, components[:, 1]),
                          (generated.d_omega, components[:, 2])):
        assert np.max(np.abs(actual - saved)) < 1e-12
    replay = synthesize_irregular_response(
        RM3_H5, generated, dt=.01, end_time=400, ramp_time=100,
    )
    np.testing.assert_allclose(replay.time, matlab_wave[:, 0], atol=1e-10,
                               rtol=0)
    assert np.max(np.abs(replay.elevation - matlab_wave[:, 1])) < 1e-10

    for name, columns in (("cylinder", 4), ("valve", 5),
                          ("high_accumulator", 3),
                          ("low_accumulator", 3), ("motor", 5),
                          ("generator", 5)):
        values = _load(name)
        assert values.shape == (40_001, columns), name
        np.testing.assert_allclose(values[:, 0], matlab_wave[:, 0],
                                   atol=1e-10, rtol=0)


def test_hydraulic_cylinder_pressure_and_force_balance():
    cylinder = _load("cylinder")
    valve = _load("valve")
    pto = _load("RM3_cHydraulic_PTO_pto1")
    assert pto.shape == (40_001, 25)
    np.testing.assert_allclose(pto[:, 0], cylinder[:, 0], atol=1e-10,
                               rtol=0)
    assert np.max(np.abs(cylinder[:, 2])) > 1e3
    model = CompressibleCylinder(.0378, .0378, 1.86e9, 70, 35)
    reconstructed = model.force(cylinder[:, 1], cylinder[:, 3])
    assert np.max(np.abs(reconstructed - cylinder[:, 2])) < 1e-4

    # The source piston selects heave from the PTO response bus. Its first
    # two valve outputs are chamber port flows, and pressure uses a discrete
    # integrator at the model time step.
    rate_a, rate_b = model.pressure_rates(
        pto[:, 3], pto[:, 9], valve[:, 1], valve[:, 2],
    )
    dt = np.diff(cylinder[:, 0])
    for pressure, rate in ((cylinder[:, 1], rate_a),
                           (cylinder[:, 3], rate_b)):
        predicted = pressure[:-1] + dt * rate[:-1]
        assert np.max(np.abs(predicted - pressure[1:])) < 3e-6
