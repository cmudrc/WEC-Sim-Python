"""Geometry checks for the two published OSWEC motion converters."""

import numpy as np
import pytest

from wecsim.crank import AdjustableRodCrank, FixedRodCrank


@pytest.mark.parametrize("linkage", [
    FixedRodCrank(3, 1.3, 5),
    AdjustableRodCrank(3, 1.3, 5),
])
def test_stroke_starts_at_zero_and_jacobian_matches_finite_difference(linkage):
    assert abs(linkage.stroke_and_jacobian(0)[0]) < 1e-14
    angles = np.linspace(-.8, .8, 21)
    stroke, jacobian = linkage.stroke_and_jacobian(angles)
    delta = 1e-6
    plus = linkage.stroke_and_jacobian(angles + delta)[0]
    minus = linkage.stroke_and_jacobian(angles - delta)[0]
    np.testing.assert_allclose((plus - minus) / (2 * delta), jacobian,
                               rtol=0, atol=1e-9)
    assert np.isfinite(stroke).all()


def test_fixed_rod_rejects_unreachable_angle():
    linkage = FixedRodCrank(3, 1.3, 2)
    with pytest.raises(ValueError, match="outside"):
        linkage.stroke_and_jacobian(np.pi)
