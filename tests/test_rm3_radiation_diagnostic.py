"""Check the suspect fitted radiation against the supplied source data."""

from pathlib import Path

import numpy as np
import pytest

from tools.rm3_radiation_diagnostic import common_surge_curves
from wecsim.rm3Regular import solve_rm3_regular


RM3 = (Path(__file__).parent / "test_objects" / "test_bodyclass"
       / "testData" / "hydroData" / "rm3.h5")


def test_fitted_common_surge_has_opposite_low_frequency_sign():
    source_frequency, source, fitted_frequency, fitted = common_surge_curves(RM3)
    assert fitted_frequency[0] == 0
    np.testing.assert_array_equal(fitted_frequency[1:], source_frequency)
    assert np.all(source[:, source_frequency < 0.3] >= 0)
    assert np.all(fitted[:, 0] < -10_000)
    assert np.all(fitted[:, np.argmin(abs(fitted_frequency - 2 * np.pi / 8))] > 0)


@pytest.mark.parametrize("body_to_body", [False, True])
def test_default_convolution_does_not_amplify_no_wave_common_surge(body_to_body):
    """The physical default should dissipate a free surge perturbation.

    This checks an unforced response independently of the published MATLAB
    Cases 5–6, whose fitted radiation supplies negative low-frequency damping.
    """
    initial_speed = 0.01
    response = solve_rm3_regular(
        RM3, wave_height=0, no_wave=True, b2b=body_to_body,
        radiation_memory=60, end_time=400,
        initial_speed=np.array([initial_speed, 0, 0, 0]),
    )
    surge_speed = response.body_velocity[:, 0, 0]
    np.testing.assert_allclose(surge_speed[0], initial_speed, rtol=0, atol=1e-12)
    assert np.max(np.abs(surge_speed)) <= initial_speed * (1 + 1e-6)
    assert 0 < surge_speed[-1] < initial_speed
