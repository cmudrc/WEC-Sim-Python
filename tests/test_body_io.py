"""Read the RM3 hydrodynamic file included with the original port."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "source" / "objects"))

from bodyClass import BodyClass  # noqa: E402


def test_rm3_hydrodynamic_bodies_load():
    h5_file = ROOT / "source/objects/rm3.h5"
    for number, expected_name in ((1, "float"), (2, "spar")):
        body = BodyClass(str(h5_file))
        body.bodyNumber = number
        body.readH5file()
        assert body.name == expected_name
        assert body.hydroData["properties"]["name"] == expected_name
        assert body.hydroData["simulation_parameters"]["water_depth"] == "infinite"
        assert body.hydroData["hydro_coeffs"]["excitation"]["re"].shape[0] == 6
