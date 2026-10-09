"""Pair the published RM3 and OSWEC wave-marker elevation histories."""

import os
from pathlib import Path

import numpy as np
import pytest

from wecsim import RegularWave


REFERENCE = os.environ.get("WEC_SIM_MATLAB_MODEL_OUTPUT_DIR")
pytestmark = pytest.mark.skipif(
    not REFERENCE, reason="paired MATLAB wave-marker output not provided",
)


@pytest.mark.parametrize("family", ["RM3", "OSWEC"])
def test_published_regular_wave_markers(family):
    reference = Path(REFERENCE)
    prefix = f"WAVE_MARKERS_{family}"
    locations = np.loadtxt(reference / f"{prefix}_locations.csv", delimiter=",")
    markers = np.loadtxt(reference / f"{prefix}_markers.csv", delimiter=",")
    origin = np.loadtxt(reference / f"{prefix}_origin.csv", delimiter=",")
    height, period, direction, ramp_time, gravity, depth, source_k = np.loadtxt(
        reference / f"{prefix}_wave_parameters.csv", delimiter=",",
    )

    grid = np.arange(-20, 21, 10)
    x, y = np.meshgrid(grid, grid)
    expected_locations = np.column_stack((x.ravel(order="F"), y.ravel(order="F")))
    np.testing.assert_array_equal(locations, expected_locations)
    assert markers.shape == (4001, 26)
    assert origin.shape == (4001, 2)
    assert np.isfinite(markers).all()
    np.testing.assert_allclose(markers[:, 0], np.arange(4001) * 0.1,
                               rtol=0, atol=1e-10)
    np.testing.assert_allclose(markers[:, 0], origin[:, 0], rtol=0, atol=1e-10)

    wave = RegularWave(height, period, direction)
    actual = wave.elevation_at(
        markers[:, 0], locations, water_depth=depth,
        ramp_time=ramp_time, g=gravity,
    )
    np.testing.assert_allclose(actual, markers[:, 1:], rtol=0, atol=1e-10)
    center = np.flatnonzero(np.all(locations == [0, 0], axis=1))
    assert center.size == 1
    np.testing.assert_allclose(actual[:, center[0]], origin[:, 1],
                               rtol=0, atol=1e-10)
    # Source k is retained to distinguish a dispersion mismatch from a
    # marker-phase mismatch if the pair ever fails.
    omega = 2 * np.pi / period
    residual = omega**2 - gravity * source_k * np.tanh(source_k * depth)
    assert abs(residual) < 1e-11
