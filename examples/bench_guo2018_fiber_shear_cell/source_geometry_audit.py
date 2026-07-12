#!/usr/bin/env python3
"""Keep unreported wall discretisation separate from source-reported physics.

Guo et al. specify a *class* of boundary (rigidly connected spheres, with
specified blades and motion), but do not give a sphere radius or tessellation.
That makes an exact byte-for-byte reconstruction unidentifiable.  It does not
license an arbitrary value or a calibration to Fig. 6/7.  This module makes
the useful middle ground executable: a future DIRT campaign may use only an
independently justified realisation or a predeclared wall-resolution
sensitivity ensemble, and neither is labelled source-equivalent.
"""
import argparse
import json
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONTRACT = HERE / "data" / "source_geometry_contract.json"


def load_contract(path=CONTRACT):
    record = json.loads(Path(path).read_text())
    required = {"reported", "not_reported", "source_pages"}
    if set(record) != required:
        raise ValueError("geometry contract schema changed")
    reported = {"sphere_built_walls", "upper_blade_length_mm", "lower_blade_length_mm",
                "blade_pitch_mm", "upper_wall_constraints", "periodic_cell_mm"}
    if not reported <= set(record["reported"]):
        raise ValueError("geometry contract omits a reported source constraint")
    missing = {"wall_sphere_diameter_mm", "wall_sphere_layout"}
    if set(record["not_reported"]) != missing:
        raise ValueError("geometry contract must name every unresolved boundary quantity")
    if any(record["not_reported"][key] != "NOT_REPORTED" for key in missing):
        raise ValueError("unpublished wall data must not be guessed as source facts")
    return record


def validate_wall_realisation(manifest):
    """Validate provenance, not physical agreement, for one proposed wall mesh.

    The only allowed ways to choose an unreported discretisation are a
    separately citable measurement/artifact or a resolution sweep declared
    before looking at any Guo observable.  A solver result is intentionally
    not accepted as a selection basis.
    """
    item = json.loads(Path(manifest).read_text())
    required = {"label", "source_equivalent", "selection_basis", "diameter_mm",
                "layout", "reference_observables_consulted"}
    if set(item) != required:
        raise ValueError("wall realisation schema changed")
    if item["source_equivalent"]:
        raise ValueError("unreported wall discretisation cannot be labelled source-equivalent")
    if item["selection_basis"] not in {"independent_measurement", "archival_artifact", "predeclared_sensitivity"}:
        raise ValueError("wall choice needs independent evidence or a predeclared sensitivity basis")
    if not isinstance(item["diameter_mm"], (int, float)) or item["diameter_mm"] <= 0:
        raise ValueError("wall sphere diameter must be a positive physical value")
    if not isinstance(item["layout"], str) or not item["layout"].strip():
        raise ValueError("wall layout must be declared")
    if item["reference_observables_consulted"]:
        raise ValueError("wall choice must not consult Fig. 6/7 observables")
    return item


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--wall-realisation", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        raise SystemExit(not unittest.main(argv=["audit"], exit=False).result.wasSuccessful())
    if not args.verify:
        parser.error("--verify is required")
    record = load_contract()
    print("PASS: source boundary class is recorded; wall discretisation is explicitly unreported.")
    print("reported:", ", ".join(sorted(record["reported"])))
    if args.wall_realisation:
        item = validate_wall_realisation(args.wall_realisation)
        print(f"PASS: {item['label']} is a non-equivalent, non-calibrated wall realisation.")


class AuditTests(unittest.TestCase):
    def test_contract_is_explicit_about_missing_discretisation(self):
        record = load_contract()
        self.assertEqual(set(record["not_reported"]), {"wall_sphere_diameter_mm", "wall_sphere_layout"})

    def test_observable_fitting_cannot_select_a_wall(self):
        import tempfile
        candidate = {"label": "bad", "source_equivalent": False,
                     "selection_basis": "predeclared_sensitivity", "diameter_mm": 2.4,
                     "layout": "square", "reference_observables_consulted": ["Fig. 6"]}
        with tempfile.NamedTemporaryFile("w") as f:
            json.dump(candidate, f); f.flush()
            with self.assertRaisesRegex(ValueError, "must not consult"):
                validate_wall_realisation(f.name)

    def test_source_equivalence_is_not_claimed_for_unknown_mesh(self):
        import tempfile
        candidate = {"label": "bad", "source_equivalent": True,
                     "selection_basis": "archival_artifact", "diameter_mm": 2.4,
                     "layout": "square", "reference_observables_consulted": []}
        with tempfile.NamedTemporaryFile("w") as f:
            json.dump(candidate, f); f.flush()
            with self.assertRaisesRegex(ValueError, "cannot be labelled"):
                validate_wall_realisation(f.name)


if __name__ == "__main__":
    main()
