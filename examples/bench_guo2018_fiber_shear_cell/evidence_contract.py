#!/usr/bin/env python3
"""Bind Guo methods and Fig. 6/7 numbers to a primary-paper artifact.

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


def is_pdf(path: Path) -> bool:
    """Reject a renamed arbitrary file before accepting its recorded hash."""
    return path.read_bytes()[:5] == b"%PDF-"


def verify_primary_reference(source_pdf: Path, provenance: Path = PROVENANCE) -> dict:
    """Return an authenticated primary-source record or raise.

    A DOI proves bibliographic identity only. It cannot authenticate either a
    digitized result or a transcription of the numerical geometry/material
    table, so every route that calls something a Guo case must cross this
    boundary first.
    """
    if not source_pdf.is_file():
        raise ValueError(f"primary PDF is unavailable: {source_pdf}")
    record = json.loads(provenance.read_text())
    missing = [key for key in ("doi", "primary_pdf_sha256", "digitizer", "digitized_on",
                               "method_pages", "figure_pages")
               if not record.get(key) or record[key] == "UNVERIFIED"]
    if missing:
        raise ValueError("reference provenance is incomplete: " + ", ".join(missing))
    if record["doi"].lower() != "10.1002/aic.16397":
        raise ValueError("reference provenance names a different DOI")
    if sha256(source_pdf) != record["primary_pdf_sha256"]:
        raise ValueError("primary PDF hash does not match the committed digitization record")
    if not is_pdf(source_pdf):
        raise ValueError("recorded primary artifact is not a PDF")
    if not all(record["method_pages"].get(key) for key in
               ("control_cell", "wall_construction", "normal_load", "shear_protocol")):
        raise ValueError("receipt lacks page evidence for every source-derived method claim")
    if not all(record["figure_pages"].get(str(figure)) for figure in (6, 7)):
        raise ValueError("receipt lacks page evidence for both digitized figures")
    return record


def verify_reference_data(source_pdf: Path, provenance: Path = PROVENANCE,
                          reference: Path = REFERENCE) -> dict:
    """Verify the result digitization after authenticating its primary PDF."""
    record = verify_primary_reference(source_pdf, provenance)
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
    def test_rejects_arbitrary_pdf_instead_of_receipted_primary(self):
        with tempfile.NamedTemporaryFile() as pdf:
            with self.assertRaisesRegex(ValueError, "hash"):
                verify_reference_data(Path(pdf.name))

    def test_rejects_arbitrary_pdf_before_any_method_transcription(self):
        with tempfile.NamedTemporaryFile() as pdf:
            with self.assertRaisesRegex(ValueError, "hash"):
                verify_primary_reference(Path(pdf.name))

    def test_rejects_pdf_that_does_not_match_a_completed_record(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pdf = root / "paper.pdf"; pdf.write_bytes(b"primary")
            wrong = root / "wrong.pdf"; wrong.write_bytes(b"other")
            provenance = root / "provenance.json"
            provenance.write_text(json.dumps({"doi": "10.1002/aic.16397",
                "primary_pdf_sha256": sha256(pdf), "digitizer": "initials", "digitized_on": "2026-07-12",
                "method_pages": {"control_cell": [5], "wall_construction": [5], "normal_load": [6], "shear_protocol": [6]},
                "figure_pages": {"6": [8], "7": [9]}}))
            reference = root / "reference.csv"
            reference.write_text("observable,normal_stress_pa,value,figure,material,source\n"
                "shear_stress_pa,651,480,6,rubber_cord,experiment\n"
                "solid_fraction,651,0.33,7,rubber_cord,experiment\n")
            with self.assertRaisesRegex(ValueError, "hash"):
                verify_reference_data(wrong, provenance, reference)

    def test_rejects_hash_matched_nonpaper_artifact(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            artifact = root / "not_the_paper.pdf"
            artifact.write_bytes(b"not a PDF")
            provenance = root / "provenance.json"
            provenance.write_text(json.dumps({"doi": "10.1002/aic.16397",
                "primary_pdf_sha256": sha256(artifact), "digitizer": "initials", "digitized_on": "2026-07-12",
                "method_pages": {"control_cell": [5], "wall_construction": [5], "normal_load": [6], "shear_protocol": [6]},
                "figure_pages": {"6": [8], "7": [9]}}))
            with self.assertRaisesRegex(ValueError, "not a PDF"):
                verify_primary_reference(artifact, provenance)


if __name__ == "__main__":
    main()
