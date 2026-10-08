"""Published Desalination cylinder and hydraulic junction source checks."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import IdealDoubleActingCylinder


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="pinned OSWEC Desalination output absent",
)


def _source(name):
    return np.loadtxt(
        Path(REFERENCE) / f"OSWEC_DESALINATION_SOURCE_{name}.csv",
        delimiter=",",
    )


def test_published_cylinder_force_and_port_flows():
    state = _source("cylinder")
    mechanical = _source("simout")
    pto = _source("Desalination_pto1")
    assert state.shape == (30_001, 8)
    np.testing.assert_array_equal(state[:, 0], mechanical[:, 0])
    cylinder = IdealDoubleActingCylinder(area_a=.26, area_b=.26)
    rod_force = cylinder.force(state[:, 3], state[:, 4])
    flow_a, flow_b = cylinder.port_flows(state[:, 2])
    np.testing.assert_allclose(rod_force, state[:, 7], rtol=0, atol=1e-7)
    np.testing.assert_allclose(flow_a, state[:, 5], rtol=0, atol=1e-12)
    np.testing.assert_allclose(flow_b, state[:, 6], rtol=0, atol=1e-12)
    np.testing.assert_allclose(-rod_force, mechanical[:, 2],
                               rtol=0, atol=1e-7)
    np.testing.assert_allclose(state[:, 2], -pto[:, 9],
                               rtol=0, atol=1e-12)
    assert np.max(np.abs(state[:, 1] - 6)) < 6
    # The legacy source permits chamber gauge pressures below vacuum.
    assert min(np.min(state[:, 3]), np.min(state[:, 4])) < -1e7


def test_published_desalination_junction_flow_balance():
    hydraulic = _source("simout1")
    relief = _source("relief")
    accumulator = _source("accumulator")
    assert hydraulic.shape == (30_001, 8)
    assert relief.shape == (30_001, 4)
    assert accumulator.shape == (30_001, 4)
    np.testing.assert_array_equal(hydraulic[:, 0], relief[:, 0])
    np.testing.assert_allclose(accumulator[:, 1], hydraulic[:, 7],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(accumulator[:, 3], hydraulic[:, 6],
                               rtol=0, atol=1e-8)
    np.testing.assert_allclose(relief[:, 1], hydraulic[:, 6],
                               rtol=0, atol=1e-8)
    np.testing.assert_allclose(hydraulic[:, 5],
                               hydraulic[:, 2] + hydraulic[:, 3],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(hydraulic[:, 4], .95 * hydraulic[:, 3],
                               rtol=0, atol=1e-12)
    np.testing.assert_allclose(
        relief[:, 3], hydraulic[:, 1] + hydraulic[:, 4]
        - hydraulic[:, 5] - hydraulic[:, 7], rtol=0, atol=1e-12,
    )
