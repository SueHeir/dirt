#!/usr/bin/env python3
"""Guard against deriving a wall mesh from the non-metric Figure 2 rendering.

Figure 2(b) identifies the wall material, but it provides no scale bar, camera
calibration, or sphere count/layout annotation.  It is qualitative evidence
for the boundary class, not a measurement of the constituent wall spheres.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import unittest


ROOT = pathlib.Path(__file__).resolve().parent
CONTRACT = ROOT / "data" / "source_geometry_contract.json"
UNREPORTED = "NOT_REPORTED"


def load_contract(path: pathlib.Path = CONTRACT) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def audit(contract: dict) -> None:
    """Accept only the source-supported conclusion for Figure 2 wall geometry."""
    reported = contract["reported"]
    missing = contract["not_reported"]
    if reported.get("sphere_built_walls") != "Walls are rigidly connected spheres.":
        raise ValueError("Figure 2 boundary class must remain rigidly connected spheres")
    if missing.get("figure_2_wall_mesh_scale") != UNREPORTED:
        raise ValueError(
            "Figure 2 has no metric wall-mesh scale; do not infer a wall sphere diameter "
            "or lattice from rendered pixels"
        )
    for key in ("wall_sphere_diameter_mm", "wall_sphere_layout"):
        if missing.get(key) != UNREPORTED:
            raise ValueError(
                f"{key} is not reported by the primary source; a Figure 2 visual estimate "
                "cannot certify a source-equivalent wall"
            )


class FigureScaleAuditTests(unittest.TestCase):
    def test_receipted_contract_keeps_the_figure_nonmetric(self) -> None:
        audit(load_contract())

    def test_rejects_pixel_derived_wall_diameter(self) -> None:
        contract = load_contract()
        contract["not_reported"]["figure_2_wall_mesh_scale"] = "1.2 mm inferred from image pixels"
        with self.assertRaisesRegex(ValueError, "no metric wall-mesh scale"):
            audit(contract)

    def test_rejects_layout_inferred_from_the_rendering(self) -> None:
        contract = load_contract()
        contract["not_reported"]["wall_sphere_layout"] = "square lattice inferred from Figure 2"
        with self.assertRaisesRegex(ValueError, "not reported"):
            audit(contract)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true", help="verify the receipted source contract")
    parser.add_argument("--self-test", action="store_true", help="run adversarial contract tests")
    args = parser.parse_args()
    if not (args.verify or args.self_test):
        parser.error("choose --verify and/or --self-test")
    if args.verify:
        audit(load_contract())
        print("PASS: Figure 2 establishes sphere-built walls but no metric wall-mesh geometry.")
    if args.self_test:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(FigureScaleAuditTests)
        result = unittest.TextTestRunner(verbosity=2).run(suite)
        if not result.wasSuccessful():
            raise SystemExit(1)


if __name__ == "__main__":
    main()
