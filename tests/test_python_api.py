"""Exercise the public Python WEC builder without JSON or CSV files."""

import json
from pathlib import Path

import numpy as np
import pytest

from examples.configurable_rm3_pto import HYDRO, build_wec
from wecsim.caseDynamics import run_case
from wecsim import LatchingControl, NoWave, RegularWave, WEC, WorldPoint


ROOT = Path(__file__).resolve().parents[1]
JSON_EXAMPLE = ROOT / "examples/configurable_rm3_pto.json"


def test_python_builder_matches_equivalent_case():
    wec = build_wec()
    wave = RegularWave(height=2.5, period=8)
    result = wec.run(wave, dt=0.1, end_time=0.2, ramp_time=5)
    case = json.loads(JSON_EXAMPLE.read_text(encoding="utf-8"))
    case["simulation"]["end_time"] = 0.2
    reference = run_case(case, base_dir=JSON_EXAMPLE.parent)
    np.testing.assert_allclose(result.time, reference.time)
    np.testing.assert_allclose(result.raw.body_position, reference.body_position)
    np.testing.assert_allclose(
        result.ptos["main"].force,
        dict(reference.extra_outputs)["pto_main_force"],
    )
    np.testing.assert_allclose(
        result.coordinates["shared_pitch"].position,
        dict(reference.extra_outputs)["coordinate_shared_pitch_position"],
    )
    assert result.bodies["float"].position.shape == (3, 6)
    assert result.ptos["main"].absorbed_power.shape == (3,)
    case_mapping = wec.to_case(wave, dt=0.1, end_time=0.2, ramp_time=5)
    assert "coordinates" in case_mapping["constraint"]
    assert "coordinate_map" not in case_mapping["bodies"][0]
    json.dumps(case_mapping)


def test_python_builder_supports_ground_anchor_and_named_initial_state():
    wec = WEC("Anchored float")
    body = wec.body("float", HYDRO)
    wec.coordinate("heave", body.move("heave"))
    wec.pto("spring", body.at(0, 0, 0), WorldPoint(0, 0, 10),
            stiffness=1000)
    result = wec.run(
        NoWave(), dt=0.1, end_time=0.2, radiation_memory=0.2,
        initial_coordinate={"heave": 1},
    )
    np.testing.assert_allclose(result.ptos["spring"].stroke[0], -1)
    np.testing.assert_allclose(result.ptos["spring"].force[0], 1000)
    assert np.isfinite(result.bodies["float"].position).all()


def test_python_builder_rejects_foreign_body_attachment():
    first = WEC("first")
    second = WEC("second")
    first.body("body", HYDRO)
    other = second.body("body", HYDRO)
    with pytest.raises(ValueError, match="must belong to this WEC"):
        first.pto("bad", other.at(0, 0, 0), WorldPoint(0, 0, 2),
                  damping=1)


def test_python_builder_exposes_latching_pto_settings_and_force():
    wec = WEC("Latching float")
    body = wec.body("float", HYDRO)
    wec.coordinate("heave", body.move("heave"))
    control = LatchingControl(gain=100, latch_damping=10_000,
                             latch_time=0.2)
    wec.pto("latch", body.at(0, 0, 0), WorldPoint(0, 0, 10),
            axis=(0, 0, 1), control=control)
    wave = RegularWave(2.5, 8)
    case = wec.to_case(wave, dt=0.01, end_time=0.1,
                       initial_speed={"heave": 0.1})
    assert case["ptos"][0]["control"] == {
        "kind": "latching", "latch_time": 0.2, "latch_damping": 10_000,
        "minimum_normal_time": 0.2,
    }
    result = wec.run(wave, dt=0.01, end_time=0.1, ramp_time=1,
                     initial_speed={"heave": 0.1})
    np.testing.assert_allclose(result.ptos["latch"].force[0],
                               -100 * result.ptos["latch"].velocity[0])
    assert np.isfinite(result.bodies["float"].position).all()
