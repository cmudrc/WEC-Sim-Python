"""Extract the pinned MOST MATLAB Function source from its Simulink library."""

from pathlib import Path
import sys
from xml.etree import ElementTree
from zipfile import ZipFile


def main(library: str, output: str) -> None:
    with ZipFile(library) as archive:
        chart = ElementTree.fromstring(
            archive.read("simulink/stateflow/chart_27.xml"),
        )
    script = chart.findtext(".//P[@Name='script']")
    if not script or not script.startswith("function F_aero_r_bl = BEM("):
        raise ValueError("pinned MOST library BEM function was not found")
    destination = Path(output)
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "BEM.m").write_text(script)


if __name__ == "__main__":
    main(*sys.argv[1:])
