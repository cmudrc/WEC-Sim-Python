"""Check the suspect fitted radiation against the supplied source data."""

from pathlib import Path

import numpy as np

from tools.rm3_radiation_diagnostic import common_surge_curves


RM3 = (Path(__file__).parent / "test_objects" / "test_bodyclass"
       / "testData" / "hydroData" / "rm3.h5")


def test_fitted_common_surge_has_opposite_low_frequency_sign():
    source_frequency, source, fitted_frequency, fitted = common_surge_curves(RM3)
    assert fitted_frequency[0] == 0
    np.testing.assert_array_equal(fitted_frequency[1:], source_frequency)
    assert np.all(source[:, source_frequency < 0.3] >= 0)
    assert np.all(fitted[:, 0] < -10_000)
    assert np.all(fitted[:, np.argmin(abs(fitted_frequency - 2 * np.pi / 8))] > 0)
