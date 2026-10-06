"""Run the RM3 coupled solver with the fork's bundled HDF5 input."""

from pathlib import Path

import numpy as np

from source.objects.rm3Regular import solve_rm3_regular


def test_rm3_coupled_regular_wave_pipeline():
    h5_file = Path(__file__).resolve().parents[1] / "source" / "objects" / "rm3.h5"
    response = solve_rm3_regular(h5_file, end_time=40, ramp_time=10)
    assert response.body_position.shape == (401, 2, 6)
    assert response.body_velocity.shape == (401, 2, 6)
    assert response.pto_force.shape == (401,)
    assert np.isfinite(response.body_position).all()
    assert np.isfinite(response.body_velocity).all()
    assert np.isfinite(response.pto_force).all()
    assert np.max(np.abs(response.body_position[:, 0, 0])) > 0
    np.testing.assert_allclose(response.body_position[:, 0, 4], response.body_position[:, 1, 4])
