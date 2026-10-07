"""Run the RM3 coupled solver with the fork's bundled HDF5 input."""

from pathlib import Path

import numpy as np
import pytest

from wecsim.rm3Regular import solve_rm3_regular


@pytest.mark.parametrize("b2b", [False, True])
def test_rm3_coupled_regular_wave_pipeline(b2b):
    h5_file = Path(__file__).resolve().parents[1] / "examples" / "data" / "rm3.h5"
    response = solve_rm3_regular(h5_file, end_time=40, ramp_time=10, b2b=b2b)
    assert response.body_position.shape == (401, 2, 6)
    assert response.body_velocity.shape == (401, 2, 6)
    assert response.pto_force.shape == (401,)
    assert np.isfinite(response.body_position).all()
    assert np.isfinite(response.body_velocity).all()
    assert np.isfinite(response.pto_force).all()
    assert np.max(np.abs(response.body_position[:, 0, 0])) > 0
    np.testing.assert_allclose(response.body_position[:, 0, 4], response.body_position[:, 1, 4])


def test_fir_requires_radiation_memory():
    h5_file = Path(__file__).resolve().parents[1] / "examples" / "data" / "rm3.h5"
    with pytest.raises(ValueError, match="convolution and FIR need it"):
        solve_rm3_regular(h5_file, radiation_method="fir", end_time=0)
