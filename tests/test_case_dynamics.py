"""Exercise the case input boundary independently of MATLAB availability."""

import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from source.objects.caseDynamics import run_case


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples/python/rm3.json"


def test_example_case_runs_and_records_provenance(tmp_path):
    case = json.loads(EXAMPLE.read_text(encoding="utf-8"))
    case["simulation"]["end_time"] = 2
    case_file = tmp_path / "case.json"
    for body in case["bodies"]:
        body["hydro_file"] = str((EXAMPLE.parent / body["hydro_file"]).resolve())
    case_file.write_text(json.dumps(case), encoding="utf-8")
    output = tmp_path / "motion.csv"
    subprocess.run(
        [sys.executable, "-m", "source.objects.wecSimPython",
         str(case_file), "--output", str(output)],
        cwd=ROOT, check=True, capture_output=True, text=True,
    )
    values = np.loadtxt(output, delimiter=",", skiprows=1)
    metadata = json.loads(output.with_suffix(".json").read_text(encoding="utf-8"))
    assert values.shape == (21, 27)
    assert np.isfinite(values).all()
    assert metadata["time_steps"] == 21
    assert len(metadata["hydro_files"]) == 1
    assert len(metadata["hydro_files"][0]["sha256"]) == 64


@pytest.mark.parametrize("change,explanation", [
    ({"constraint": {"kind": "free_six_dof"}}, "unsupported constraint layout"),
    ({"pto": {"kind": "pitch", "damping": 1}}, "relative_heave PTO"),
    ({"wave": {"type": "irregular"}}, "regular waves"),
])
def test_unsupported_physics_fails_explicitly(change, explanation):
    case = json.loads(EXAMPLE.read_text(encoding="utf-8"))
    case["simulation"]["end_time"] = 0
    case.update(change)
    with pytest.raises(ValueError, match=explanation):
        run_case(case, base_dir=EXAMPLE.parent)
