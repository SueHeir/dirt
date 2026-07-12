#!/usr/bin/env python3
"""Bind Guo Fig. 6/7 numbers to a primary-paper artifact before comparison.

This is deliberately an evidence gate, not a physics or solver gate.  A CSV
that is only internally consistent with a printed regression fit is not an
independent experimental reference.  The gate requires the actual PDF used to
digitize the points and records its content hash in version control.
"""
import argparse
import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROVENANCE = HERE / "data" / "reference_provenance.json"
REFERENCE = HERE / "data" / "guo_2019_rubber_cord.csv"
REQUIRED_COLUMNS = {"observable", "normal_stress_pa", "value", "figure", "material", "source"}
REQUIRED_FIGURES = {"6", "7"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_reference_data(source_pdf: Path, provenance: Path = PROVENANCE,
                          reference: Path = REFERENCE) -> dict:
    """Return audited provenance or raise; never infer it from DIRT output."""
    if not source_pdf.is_file():
        raise ValueError(f"primary PDF is unavailable: {source_pdf}")
    record = json.loads(provenance.read_text())
    missing = [key for key in ("doi", "primary_pdf_sha256", "digitizer", "digitized_on")
               if not record.get(key) or record[key] == "UNVERIFIED"]
    if missing:
        raise ValueError("reference provenance is incomplete: " + ", ".join(missing))
    if record["doi"].lower() != "10.1002/aic.16397":
        raise ValueError("reference provenance names a different DOI")
    if sha256(source_pdf) != record["primary_pdf_sha256"]:
        raise ValueError("primary PDF hash does not match the committed digitization record")
    with reference.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows or set(rows[0]) != REQUIRED_COLUMNS:
        raise ValueError("reference CSV does not have the declared schema")
    if {row["figure"] for row in rows} != REQUIRED_FIGURES:
        raise ValueError("reference CSV is not explicitly tied to both Fig. 6 and Fig. 7")
    if any(row["source"] != "experiment" or row["material"] != "rubber_cord" for row in rows):
        raise ValueError("reference CSV includes a non-experimental or wrong-material point")
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-pdf", type=Path, required=False)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(EvidenceContractTests)
        raise SystemExit(not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful())
    if args.source_pdf is None:
        parser.error("--source-pdf is required; a DOI or fitted line is not primary figure evidence")
    record = verify_reference_data(args.source_pdf)
    print(f"PASS: Fig. 6/7 reference is bound to {record['doi']} PDF SHA-256")


class EvidenceContractTests(unittest.TestCase):
    def test_rejects_inherited_unverified_transcription(self):
        with tempfile.NamedTemporaryFile() as pdf:
            with self.assertRaisesRegex(ValueError, "provenance is incomplete"):
                verify_reference_data(Path(pdf.name))

    def test_rejects_pdf_that_does_not_match_a_completed_record(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pdf = root / "paper.pdf"; pdf.write_bytes(b"primary")
            wrong = root / "wrong.pdf"; wrong.write_bytes(b"other")
            provenance = root / "provenance.json"
            provenance.write_text(json.dumps({"doi": "10.1002/aic.16397",
                "primary_pdf_sha256": sha256(pdf), "digitizer": "initials", "digitized_on": "2026-07-12"}))
            reference = root / "reference.csv"
            reference.write_text("observable,normal_stress_pa,value,figure,material,source\n"
                "shear_stress_pa,651,480,6,rubber_cord,experiment\n"
                "solid_fraction,651,0.33,7,rubber_cord,experiment\n")
            with self.assertRaisesRegex(ValueError, "hash"):
                verify_reference_data(wrong, provenance, reference)


if __name__ == "__main__":
    main()
