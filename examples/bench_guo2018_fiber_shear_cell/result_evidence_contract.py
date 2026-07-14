#!/usr/bin/env python3
"""Authenticate future non-equivalent Guo sensitivity artifacts.

This is a provenance gate, not a numerical validator.  It cannot establish
that a solver is physically correct or that a digitised point is raw data.
It does make a later campaign reproducible enough to audit: every declared
case must bind its input, history, and summary bytes to one receipt, retain the
primary-paper hash, and identify the generated non-equivalent wall mesh.
"""
import argparse
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROVENANCE = HERE / "data" / "reference_provenance.json"
ARTIFACTS = ("solver_input", "cell_history", "observable_summary")


def digest(path: Path) -> str:
    if not path.is_file():
        raise ValueError(f"required artifact is missing: {path}")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def primary_pdf_hash(provenance: Path = PROVENANCE) -> str:
    record = json.loads(provenance.read_text(encoding="utf-8"))
    value = record.get("primary_pdf_sha256")
    if not isinstance(value, str) or len(value) != 64:
        raise ValueError("primary-source receipt lacks a SHA-256")
    return value


def verify_case(case: Path, wall: dict, stress: int, provenance: Path = PROVENANCE) -> dict:
    """Verify bytes and identity for one predeclared non-equivalent case."""
    receipt_path = case / "provenance_receipt"
    if not receipt_path.is_file():
        raise ValueError(f"required artifact is missing: {receipt_path}")
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    required = {"schema", "case_label", "normal_stress_pa", "source_equivalent",
                "primary_pdf_sha256", "wall_manifest_sha256", "artifacts",
                "solver_identity", "ai_authorship_and_limits"}
    if set(receipt) != required or receipt["schema"] != 1:
        raise ValueError("provenance receipt schema changed")
    if receipt["case_label"] != wall["label"] or receipt["normal_stress_pa"] != stress:
        raise ValueError("provenance receipt does not identify the preregistered case")
    if receipt["source_equivalent"]:
        raise ValueError("a wall-resolution sensitivity result cannot claim source equivalence")
    if receipt["primary_pdf_sha256"] != primary_pdf_hash(provenance):
        raise ValueError("result receipt is not bound to the receipted primary paper")
    if not isinstance(receipt["solver_identity"], dict) or set(receipt["solver_identity"]) != {"engine", "revision"} \
            or not all(isinstance(value, str) and value.strip() for value in receipt["solver_identity"].values()):
        raise ValueError("result receipt needs a nonempty solver engine and revision")
    if not isinstance(receipt["ai_authorship_and_limits"], str) or not receipt["ai_authorship_and_limits"].strip():
        raise ValueError("result receipt must disclose AI authorship and validation limits")
    expected = {name: digest(case / name) for name in ARTIFACTS}
    if receipt["artifacts"] != expected:
        raise ValueError("result receipt hashes do not match the retained case artifacts")
    wall_manifest = case / "wall_realisation.json"
    if receipt["wall_manifest_sha256"] != digest(wall_manifest):
        raise ValueError("result receipt is not bound to its retained wall manifest")
    wall_record = json.loads(wall_manifest.read_text(encoding="utf-8"))
    for key in ("label", "diameter_mm", "layout", "selection_basis", "source_equivalent"):
        if wall_record.get(key) != wall.get(key):
            raise ValueError("retained wall manifest differs from the preregistered wall")
    return receipt


class ResultEvidenceTests(unittest.TestCase):
    def make_case(self, root: Path, *, alter_hash=False, source_equivalent=False):
        wall = {"label": "wall", "diameter_mm": 1.2, "layout": "square",
                "selection_basis": "predeclared_sensitivity", "source_equivalent": False}
        case = root / "wall" / "p651"; case.mkdir(parents=True)
        for name in ARTIFACTS:
            (case / name).write_text(name, encoding="utf-8")
        (case / "wall_realisation.json").write_text(json.dumps(wall), encoding="utf-8")
        provenance = root / "provenance.json"
        provenance.write_text(json.dumps({"primary_pdf_sha256": "a" * 64}), encoding="utf-8")
        receipt = {"schema": 1, "case_label": "wall", "normal_stress_pa": 651,
                   "source_equivalent": source_equivalent, "primary_pdf_sha256": "a" * 64,
                   "wall_manifest_sha256": digest(case / "wall_realisation.json"),
                   "artifacts": {name: digest(case / name) for name in ARTIFACTS},
                   "solver_identity": {"engine": "DIRT", "revision": "abc"},
                   "ai_authorship_and_limits": "AI assisted setup; no independent physics validation."}
        if alter_hash:
            receipt["artifacts"]["cell_history"] = "0" * 64
        (case / "provenance_receipt").write_text(json.dumps(receipt), encoding="utf-8")
        return case, wall, provenance

    def test_receipt_binds_retained_artifacts(self):
        with tempfile.TemporaryDirectory() as directory:
            case, wall, provenance = self.make_case(Path(directory))
            self.assertEqual(verify_case(case, wall, 651, provenance)["case_label"], "wall")

    def test_changed_history_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            case, wall, provenance = self.make_case(Path(directory), alter_hash=True)
            with self.assertRaisesRegex(ValueError, "hashes do not match"):
                verify_case(case, wall, 651, provenance)

    def test_source_equivalence_claim_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            case, wall, provenance = self.make_case(Path(directory), source_equivalent=True)
            with self.assertRaisesRegex(ValueError, "cannot claim source equivalence"):
                verify_case(case, wall, 651, provenance)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path)
    parser.add_argument("--wall-manifest", type=Path)
    parser.add_argument("--normal-stress-pa", type=int)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        raise SystemExit(not unittest.main(argv=["result-evidence"], exit=False).result.wasSuccessful())
    if not (args.case and args.wall_manifest and args.normal_stress_pa is not None):
        parser.error("--case, --wall-manifest, and --normal-stress-pa are required")
    wall = json.loads(args.wall_manifest.read_text(encoding="utf-8"))
    verify_case(args.case, wall, args.normal_stress_pa)
    print("PASS: case artifacts are receipt-bound; this is provenance only, not physics validation.")


if __name__ == "__main__":
    main()
