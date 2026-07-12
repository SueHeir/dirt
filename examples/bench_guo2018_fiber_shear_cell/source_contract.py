#!/usr/bin/env python3
"""Check that a DIRT case matches Guo et al.'s *numerical* control cell.

Guo et al. used a Schulze ring shear tester experimentally, but deliberately
did **not** simulate its annulus.  Their numerical comparison (Computational
Set-up, pp. 5--7, doi:10.1002/aic.16397) uses a small plane-driven cell with
periodic x/z boundaries, a vertically mobile loaded lid, and 4/2-mm lid/base
blades.  Conflating the apparatus with this published DEM reduction was the
previous branch's central source error.

This checker is intentionally a protocol gate, not a result gate. It cannot
turn a topology audit or a solver history into a successful Fig. 6/7
replication. The primary source specifies walls made from rigidly connected
spheres; smooth finite planes are not an equivalent rough-wall model.
"""
import argparse
import json
import tempfile
import tomllib
import unittest
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
DOI = "10.1002/aic.16397"
CROSSREF = f"https://api.crossref.org/works/{DOI}"


def published_record() -> dict:
    """Fetch bibliographic identity only; geometry comes from the primary paper."""
    with urllib.request.urlopen(CROSSREF, timeout=30) as response:
        record = json.load(response)["message"]
    title = " ".join(record.get("title", []))
    if record.get("DOI", "").lower() != DOI or "Flexible Fiber Flows" not in title:
        raise ValueError("Crossref record does not identify Guo et al. doi:10.1002/aic.16397")
    return {"title": title, "doi": record["DOI"]}


def require_published_control_cell(config: Path) -> None:
    """Reject only deviations from the paper's stated DEM reduction.

    The current cylinder/fixed-boundary input is neither the physical SRST nor
    the periodic planar control cell used for Fig. 6/7.  Do not silently treat
    its output as an experimental comparison.
    """
    with config.open("rb") as stream:
        input_config = tomllib.load(stream)
    domain = input_config.get("domain", {})
    walls = input_config.get("wall", [])
    failures = []
    if domain.get("boundary_x") != "periodic" or domain.get("boundary_z") != "periodic":
        failures.append("the published numerical cell requires periodic x and z boundaries")
    if any(wall.get("type") == "cylinder" for wall in walls):
        failures.append("the published numerical cell has no cylindrical sidewall")
    if domain.get("x_high") != 0.064 or domain.get("z_high") != 0.036:
        failures.append("the baseline control cell must be Lx=64 mm and Lz=36 mm")
    names = {wall.get("name") for wall in walls}
    blade_names = {f"{level}_blade_{side}" for level in ("upper", "lower") for side in ("pos", "neg")}
    if not blade_names <= names:
        failures.append("the published numerical cell requires paired 4-mm upper and 2-mm lower blade faces")
    by_name = {wall.get("name"): wall for wall in walls}
    for prefix, height in (("upper_blade", 0.004), ("lower_blade", 0.002)):
        for side in ("pos", "neg"):
            wall = by_name.get(f"{prefix}_{side}", {})
            if abs((wall.get("bound_y_high", 0.0) - wall.get("bound_y_low", 0.0)) - height) > 1e-12:
                failures.append(f"{prefix} must have {height * 1000:g}-mm vertical extent")
    # The source topology must be present in solver input, not created later
    # by benchmark code. The generic input loader expands these four entries
    # into 32 finite faces before contact integration.
    for name in blade_names:
        wall = by_name.get(name, {})
        if wall.get("repeat_x_count") != 8 or abs(wall.get("repeat_x_pitch", 0.0) - 0.008) > 1e-12:
            failures.append("the source requires eight 8-mm-pitch blade positions per wall")
            break
    expected_assembly = {"lower_blade_pos": "lower_plate", "lower_blade_neg": "lower_plate",
                         "upper_blade_pos": "upper_plate", "upper_blade_neg": "upper_plate"}
    if any(by_name.get(name, {}).get("assembly") != assembly for name, assembly in expected_assembly.items()):
        failures.append("blade arrays must be rigidly attached to their translating or servo-controlled plate")
    if any(wall.get("type") == "plane" for wall in walls):
        failures.append(
            "the source uses rigidly connected sphere walls and blades; DIRT's smooth plane-wall assembly is not source-equivalent"
        )
    if failures:
        raise RuntimeError("BLOCKED: current DIRT input is not the published periodic control cell: " + "; ".join(failures))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=HERE / "config.toml")
    parser.add_argument("--verify-doi", action="store_true", help="query Crossref for bibliographic identity")
    parser.add_argument("--require-runnable", action="store_true", help="raise unless the configured numerical cell matches the source protocol")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(SourceContractTests)
        raise SystemExit(not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful())
    if not args.verify_doi and not args.require_runnable:
        parser.error("choose --verify-doi and/or --require-runnable")
    if args.verify_doi:
        record = published_record()
        print(f"EXTERNAL BIBLIOGRAPHY: {record['doi']} — {record['title']}")
        print("PRIMARY-SOURCE METHOD: pp. 5--7 use periodic x/z planar control cell, loaded mobile lid, and 4/2-mm blade arrays.")
    if args.require_runnable:
        require_published_control_cell(args.config)


class SourceContractTests(unittest.TestCase):
    def test_rejects_smooth_plane_surrogate_for_sphere_built_walls(self):
        with self.assertRaisesRegex(RuntimeError, "rigidly connected sphere walls"):
            require_published_control_cell(HERE / "config.toml")

    def test_rejects_a_blade_array_that_is_not_materialized_in_solver_input(self):
        """A four-template input must not be accepted as the 32-face source cell."""
        text = (HERE / "config.toml").read_text().replace("repeat_x_count = 8", "repeat_x_count = 1", 1)
        with tempfile.NamedTemporaryFile("w", suffix=".toml") as config:
            config.write(text)
            config.flush()
            with self.assertRaisesRegex(RuntimeError, "eight 8-mm-pitch"):
                require_published_control_cell(Path(config.name))


if __name__ == "__main__":
    main()
