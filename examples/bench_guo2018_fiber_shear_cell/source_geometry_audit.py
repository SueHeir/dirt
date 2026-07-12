#!/usr/bin/env python3
"""Audit what the Guo paper does--and does not--specify for its wall bodies.

The article identifies sphere-built walls, blade lengths/pitch, their motion
constraints, and the periodic control-cell dimensions.  It does not publish a
wall-sphere radius, tessellation, or a construction drawing
from which those quantities can be recovered.  Those missing quantities alter
the contact boundary and the gravity-loaded lid dynamics.  This script makes
that limit executable: a configuration cannot be certified as a reproduction
by filling those values from a convenient DIRT or LAMMPS default.
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
    if not {"sphere_built_walls", "upper_blade_length_mm", "lower_blade_length_mm",
            "blade_pitch_mm", "upper_wall_constraints", "periodic_cell_mm"} <= set(record["reported"]):
        raise ValueError("geometry contract omits a reported source constraint")
    missing = {"wall_sphere_diameter_mm", "wall_sphere_layout"}
    if set(record["not_reported"]) != missing:
        raise ValueError("geometry contract must name every unresolved boundary quantity")
    if any(record["not_reported"][key] != "NOT_REPORTED" for key in missing):
        raise ValueError("unpublished wall data must not be guessed")
    return record


def require_reproducible_wall_body(path=CONTRACT):
    record = load_contract(path)
    absent = ", ".join(sorted(record["not_reported"]))
    raise RuntimeError(
        "BLOCKED: the primary article does not report " + absent +
        "; a sphere-built, gravity-loaded DIRT wall body cannot be uniquely "
        "transcribed. Obtain author input or an archival supplement before "
        "calling any generated geometry a Guo replication."
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--require-reproducible-wall-body", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        raise SystemExit(not unittest.main(argv=["audit"], exit=False).result.wasSuccessful())
    if not args.verify and not args.require_reproducible_wall_body:
        parser.error("choose --verify and/or --require-reproducible-wall-body")
    if args.verify:
        record = load_contract()
        print("PASS: source geometry contract records reported constraints and unresolved wall-body data")
        print("reported:", ", ".join(sorted(record["reported"])))
    if args.require_reproducible_wall_body:
        require_reproducible_wall_body()


class AuditTests(unittest.TestCase):
    def test_contract_is_complete_and_explicit_about_missing_geometry(self):
        record = load_contract()
        self.assertEqual(set(record["not_reported"]), {"wall_sphere_diameter_mm", "wall_sphere_layout"})

    def test_no_default_can_turn_missing_wall_data_into_a_reproduction(self):
        with self.assertRaisesRegex(RuntimeError, "cannot be uniquely transcribed"):
            require_reproducible_wall_body()


if __name__ == "__main__":
    main()
