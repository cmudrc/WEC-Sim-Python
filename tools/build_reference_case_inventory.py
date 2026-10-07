"""Inventory the published WEC-Sim input files at pinned source revisions.

Usage: python tools/build_reference_case_inventory.py CORE_CHECKOUT APPLICATIONS_CHECKOUT
"""

import csv
from pathlib import Path
import re
import subprocess
import sys


MODELS = ("RM3", "OSWEC", "Sphere")


def revision(root):
    return subprocess.check_output(
        ["git", "-C", str(root), "rev-parse", "HEAD"], text=True,
    ).strip()


def describe(source, root, input_file, source_revision):
    active = "\n".join(line.split("%", 1)[0] for line in input_file.read_text().splitlines())
    hydro_files = sorted(set(re.findall(r"bodyClass\(['\"]([^'\"]*)['\"]\)", active)))
    wave_types = sorted(set(re.findall(r"waveClass\(['\"]([^'\"]+)['\"]\)", active)))
    families = [
        model for model in MODELS
        if any(model.lower() in path.lower() for path in hydro_files)
    ]
    return {
        "source": source,
        "case": str(input_file.relative_to(root).parent),
        "family": "+".join(families) if families else "Other or dynamic",
        "wave_types": ";".join(wave_types),
        "hydro_files": ";".join(hydro_files),
        "revision": source_revision,
    }


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    core, applications = map(Path, sys.argv[1:])
    core_revision = revision(core)
    applications_revision = revision(applications)
    files = sorted((core / "examples").glob("*/wecSimInputFile.m"))
    files = [path for path in files if path.parent.name in ("RM3", "OSWEC")]
    rows = [describe("WEC-Sim", core, path, core_revision) for path in files]
    rows += [
        describe("WEC-Sim_Applications", applications, path, applications_revision)
        for path in sorted(applications.rglob("wecSimInputFile.m"))
    ]
    writer = csv.DictWriter(sys.stdout, fieldnames=rows[0].keys(), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)


if __name__ == "__main__":
    main()
