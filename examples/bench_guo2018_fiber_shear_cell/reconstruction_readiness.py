#!/usr/bin/env python3
"""Audit whether the Guo shear-cell source is sufficient to reconstruct a DIRT case.

This is an independent, source-ledger check.  It verifies that the committed
candidate carries the directly reported rubber-cord parameters, distinguishes
the 96-mm population inference from a source fact, and reports the unresolved
wall-body inputs that prevent a source-equivalent solver run.  It never emits a
replication pass and it does not synthesize a wall geometry or a lid mass.
"""
import argparse
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
LEDGER = HERE / "data" / "reconstruction_ledger.json"


def load_ledger(path=LEDGER):
    ledger = json.loads(Path(path).read_text())
    required = {"source", "directly_reported", "derived_inputs",
                "unresolved_required_boundary_inputs", "page_receipt"}
    if set(ledger) != required:
        raise ValueError("reconstruction ledger schema changed")
    rubber = ledger["directly_reported"].get("rubber_cord", {})
    required_rubber = {"fiber_diameter_mm", "fiber_length_mm", "fiber_density_kg_m3",
                       "spheres_per_fiber", "bonds_per_fiber", "bond_length_mm",
                       "contact_and_bond_modulus_pa", "fiber_fiber_and_wall_friction",
                       "fiber_count", "time_step_s"}
    if set(rubber) != required_rubber:
        raise ValueError("ledger does not contain the complete Table-2 rubber-cord transcription")
    unresolved = {"wall_sphere_diameter_mm", "wall_sphere_layout", "upper_wall_mass_kg"}
    if set(ledger["unresolved_required_boundary_inputs"]) != unresolved:
        raise ValueError("ledger must name every unresolved wall-body input")
    if any(value != "not reported" for value in ledger["unresolved_required_boundary_inputs"].values()):
        raise ValueError("unresolved inputs must not be replaced by defaults")
    if ledger["derived_inputs"]["96_mm_fiber_count"]["value"] != 750:
        raise ValueError("derived size-sensitivity population changed")
    return ledger


def verify_pdf_receipt(source_pdf, ledger):
    path = Path(source_pdf)
    if not path.is_file() or path.read_bytes()[:5] != b"%PDF-":
        raise ValueError("source artifact is not a readable PDF")
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != ledger["source"]["primary_pdf_sha256"]:
        raise ValueError("source PDF does not match the independently transcribed ledger")


def candidate_parameter_mismatches(ledger):
    """Compare only direct Table-2 facts to candidate topology constants."""
    import prepare
    rubber = ledger["directly_reported"]["rubber_cord"]
    expected = {
        "DIAMETER": rubber["fiber_diameter_mm"] / 1000,
        "LENGTH": rubber["fiber_length_mm"] / 1000,
        "BEADS": rubber["spheres_per_fiber"],
        "SPACING": rubber["bond_length_mm"] / 1000,
        "DENSITY": rubber["fiber_density_kg_m3"],
        "YOUNGS_MODULUS": rubber["contact_and_bond_modulus_pa"],
        "TIMESTEP": rubber["time_step_s"],
        "FIBRES_AT_64_MM": rubber["fiber_count"],
        "CELL_Z": ledger["directly_reported"]["control_cell_xz_mm"][1] / 1000,
    }
    return [f"{name}={getattr(prepare, name)!r}, expected {value!r}"
            for name, value in expected.items() if getattr(prepare, name) != value]


def audit(source_pdf=None):
    ledger = load_ledger()
    if source_pdf is not None:
        verify_pdf_receipt(source_pdf, ledger)
    mismatches = candidate_parameter_mismatches(ledger)
    if mismatches:
        raise ValueError("candidate topology differs from direct source facts: " + "; ".join(mismatches))
    return ledger


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-pdf", type=Path, help="hash-matched primary PDF")
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        raise SystemExit(not unittest.main(argv=["readiness"], exit=False).result.wasSuccessful())
    if not args.verify:
        parser.error("--verify is required")
    ledger = audit(args.source_pdf)
    unresolved = ", ".join(sorted(ledger["unresolved_required_boundary_inputs"]))
    print("AUDIT OK: candidate constants agree with direct Table-2 and cell facts.")
    print("SOURCE LIMIT: a source-equivalent replication remains BLOCKED; missing " + unresolved + ".")
    print("The 96-mm/750-fibre input is an explicitly derived sensitivity case, not a reported population.")


class LedgerTests(unittest.TestCase):
    def test_candidate_constants_match_the_independently_transcribed_table(self):
        self.assertEqual(candidate_parameter_mismatches(load_ledger()), [])

    def test_missing_wall_body_data_cannot_be_silently_promoted_to_a_default(self):
        with tempfile.TemporaryDirectory() as directory:
            altered = Path(directory) / "ledger.json"
            ledger = load_ledger()
            ledger["unresolved_required_boundary_inputs"]["upper_wall_mass_kg"] = 1.0
            altered.write_text(json.dumps(ledger))
            with self.assertRaisesRegex(ValueError, "must not be replaced"):
                load_ledger(altered)

    def test_unmatched_pdf_is_not_evidence_for_this_transcription(self):
        with tempfile.NamedTemporaryFile() as source:
            source.write(b"%PDF-not-the-guo-paper")
            source.flush()
            with self.assertRaisesRegex(ValueError, "does not match"):
                verify_pdf_receipt(source.name, load_ledger())


if __name__ == "__main__":
    main()
