"""Run a supported WEC-Sim-Python dynamics case from a JSON input file.

This replaces the original unfinished script's commented-out simulation
stage. The case schema is documented in README.md and deliberately rejects
layouts or physics that have not been implemented and compared with MATLAB.

The project originated with Sungjun Won's WEC-Sim-Python port (2020).
"""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np

from .caseDynamics import AXES, run_case


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_state():
    root = Path(__file__).resolve().parents[2]
    try:
        revision = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=root,
            capture_output=True, text=True, check=False,
        )
        changes = subprocess.run(
            ["git", "status", "--porcelain"], cwd=root,
            capture_output=True, text=True, check=False,
        )
    except OSError:
        return None, None
    return (
        revision.stdout.strip() if revision.returncode == 0 else None,
        bool(changes.stdout) if changes.returncode == 0 else None,
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path, help="JSON dynamics input file")
    parser.add_argument("--output", type=Path, required=True, help="numeric CSV output")
    args = parser.parse_args(argv)
    case_file = args.case.expanduser().resolve(strict=True)
    case = json.loads(case_file.read_text(encoding="utf-8"))
    response = run_case(case, base_dir=case_file.parent)

    columns = ["time"]
    arrays = [response.time]
    if response.wave_elevation is not None:
        columns.append("wave_elevation")
        arrays.append(response.wave_elevation)
    for body in range(response.body_position.shape[1]):
        for quantity, data in (
            ("position", response.body_position),
            ("velocity", response.body_velocity),
        ):
            for dof, axis in enumerate(AXES):
                columns.append(f"body{body + 1}_{axis}_{quantity}")
                arrays.append(data[:, body, dof])
    if response.pto_force is not None:
        columns.append(response.pto_label)
        arrays.append(response.pto_force)
    if response.pto_generalized_force is not None:
        for coordinate in range(response.pto_generalized_force.shape[1]):
            columns.append(f"pto_coordinate{coordinate + 1}_force")
            arrays.append(response.pto_generalized_force[:, coordinate])
    if response.total_heave_force is not None:
        columns.append("body1_total_heave_force")
        arrays.append(response.total_heave_force)
    for name, values in response.extra_outputs:
        columns.append(name)
        arrays.append(values)
    values = np.column_stack(arrays)
    if not np.isfinite(values).all():
        raise RuntimeError("the dynamics produced nonfinite output")

    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(output, values, delimiter=",", header=",".join(columns),
               comments="", fmt="%.17g")
    revision, dirty = _git_state()
    metadata = {
        "case_file": str(case_file),
        "case_sha256": _sha256(case_file),
        "case": case,
        "hydro_files": [
            {"path": str(path), "sha256": _sha256(path)}
            for path in dict.fromkeys(response.hydro_files)
        ],
        "auxiliary_files": [
            {"path": str(path), "sha256": _sha256(path)}
            for path in response.auxiliary_files
        ],
        "git_commit": revision,
        "git_dirty": dirty,
        "numpy_version": np.__version__,
        "csv_columns": columns,
        "time_steps": len(values),
    }
    output.with_suffix(".json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8",
    )
    print(output)


if __name__ == "__main__":
    main()
