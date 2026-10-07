"""Exercise the published OSWEC settings without any MATLAB-generated output."""

from pathlib import Path

import numpy as np

from wecsim.hingePitch import solve_hinged_pitch_from_excitation
from wecsim.irregularWave import (
    pm_equal_energy_components, synthesize_irregular_response,
)


def test_seeded_oswec_wave_and_pitch_pipeline():
    h5_file = Path(__file__).parent / "test_objects" / "test_bodyclass" / "testData" / "hydroData" / "oswec.h5"
    settings = dict(
        significant_height=2.5, peak_period=8,
        directions=[0, 30, 90], spreading=[0.1, 0.2, 0.7], seed=7,
    )
    components = pm_equal_energy_components(h5_file, **settings)
    replay = pm_equal_energy_components(h5_file, **settings)
    np.testing.assert_array_equal(components.phase, replay.phase)
    wave = synthesize_irregular_response(
        h5_file, components, dt=0.1, end_time=40, ramp_time=10,
    )
    motion = solve_hinged_pitch_from_excitation(
        h5_file, wave.excitation_force, hinge_z=-8.9,
        body_mass=127_000, pitch_inertia=1.85e6, pto_damping=12_000,
    )
    assert len(motion.time) == 401
    assert np.isfinite(wave.elevation).all()
    assert np.isfinite(wave.excitation_force).all()
    assert np.isfinite(motion.angle).all()
    assert 0 < np.max(np.abs(motion.angle)) < 1
    np.testing.assert_allclose(motion.time, wave.time, rtol=0, atol=1e-12)
