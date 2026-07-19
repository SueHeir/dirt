#!/usr/bin/env python3
"""Focused, offline checks for reference-record admission rules."""

import json
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import reference_audit
import sweep


class _Response:
    def __init__(self, payload):
        self.payload = payload

    def __enter__(self):
        return self

    def __exit__(self, *unused):
        return False

    def read(self, *unused):
        return json.dumps(self.payload).encode("utf-8")


def _crossref_opener(request, timeout):
    del request, timeout
    return _Response({"message": {
        "DOI": "10.1016/s0378-4371(99)00183-1",
        "title": ["Rolling friction in the dynamic simulation of sandpile formation"],
        "author": [{"family": "Zhou"}],
        "published-print": {"date-parts": [[1999]]},
    }})


class ReferenceAuditTests(unittest.TestCase):
    def _write(self, record):
        handle = tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False)
        with handle:
            json.dump(record, handle)
        self.addCleanup(Path(handle.name).unlink, missing_ok=True)
        return handle.name

    @staticmethod
    def _source_locator(fragment):
        return "https://doi.org/10.1000/example#" + fragment

    def test_negative_control_has_independent_bibliographic_identity(self):
        record = reference_audit.audit_record(
            Path(__file__).parent / "external_records" / "zhou_1999_rolling_sandpile.json",
            opener=_crossref_opener)
        self.assertEqual(record["comparability"], "incompatible")

    def test_missing_external_protocol_record_is_not_admitted(self):
        """Missing provenance is rejected without fabricating a solver campaign."""
        with mock.patch.object(sweep, "PROTOCOL_REFERENCE", "/nonexistent/record.json"):
            self.assertFalse(sweep._protocol_reference_ok())

    def test_lammps_replay_export_preserves_complete_mechanical_state(self):
        """A replay must not silently discard the residual spin or velocity."""
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "prelift.csv"
            data = Path(directory) / "prelift.data"
            fields = ["x", "y", "z", "radius", "vx", "vy", "vz",
                      "omega_x", "omega_y", "omega_z"]
            with source.open("w", newline="") as handle:
                writer = __import__("csv").DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                for _ in range(sweep.HEAP_COUNT):
                    writer.writerow(dict(zip(fields, [
                        1.25e-3, -2.5e-3, 3.75e-3, 2.0e-3, 1.0e-4,
                        -2.0e-4, 3.0e-4, 4.0, -5.0, 6.0])))
            self.assertEqual(sweep._write_lammps_replay_data(source, data), sweep.HEAP_COUNT)
            self.assertLessEqual(sweep._lammps_replay_state_error(source, data), 1.0e-9)

    def test_matched_record_needs_reviewable_source_extraction(self):
        record = {
            "doi": "10.1000/example", "title": "Example", "first_author_family": "Doe",
            "published_year": 2020, "source_url": "https://doi.org/10.1000/example",
            "comparability": "matched", "sphere_radius_m": 0.002,
            "density_kg_m3": 2500.0, "restitution": 0.4, "sliding_friction": 0.16,
            "container": "cylinder", "deposition": "settled", "release": "lift",
            "angle_definition": "fit", "band_deg": [22.0, 26.0],
        }
        with self.assertRaisesRegex(ValueError, "human-auditable evidence"):
            reference_audit.load_record(self._write(record))

    def test_matched_record_cannot_attach_a_band_unrelated_to_its_measurements(self):
        record = {
            "doi": "10.1000/example", "title": "Example", "first_author_family": "Doe",
            "published_year": 2020, "source_url": "https://doi.org/10.1000/example",
            "comparability": "matched", "sphere_radius_m": 0.002,
            "density_kg_m3": 2500.0, "restitution": 0.4, "sliding_friction": 0.16,
            "container": "cylinder", "deposition": "settled", "release": "lift",
            "angle_definition": "fit", "band_deg": [22.0, 26.0],
            "source_locator": self._source_locator("table-2"), "extraction_method": "manual transcription",
            "observations": [{"value_deg": 31.0, "source_locator": self._source_locator("table-2-row-1")}],
        }
        with self.assertRaisesRegex(ValueError, "does not bound"):
            reference_audit.load_record(self._write(record))

    def test_matched_record_requires_a_locatable_numeric_observation(self):
        record = {
            "doi": "10.1000/example", "title": "Example", "first_author_family": "Doe",
            "published_year": 2020, "source_url": "https://doi.org/10.1000/example",
            "comparability": "matched", "sphere_radius_m": 0.002,
            "density_kg_m3": 2500.0, "restitution": 0.4, "sliding_friction": 0.16,
            "container": "cylinder", "deposition": "settled", "release": "lift",
            "angle_definition": "fit", "band_deg": [22.0, 26.0],
            "source_locator": self._source_locator("table-2"), "extraction_method": "manual transcription",
            "observations": [{"value_deg": 24.0}],
        }
        with self.assertRaisesRegex(ValueError, "source_locator"):
            reference_audit.load_record(self._write(record))

    def test_matched_record_accepts_locatable_measurements_inside_its_band(self):
        record = {
            "doi": "10.1000/example", "title": "Example", "first_author_family": "Doe",
            "published_year": 2020, "source_url": "https://doi.org/10.1000/example",
            "comparability": "matched", "sphere_radius_m": 0.002,
            "density_kg_m3": 2500.0, "restitution": 0.4, "sliding_friction": 0.16,
            "container": "cylinder", "deposition": "settled", "release": "lift",
            "angle_definition": "fit", "band_deg": [22.0, 26.0],
            "source_locator": self._source_locator("table-2"), "extraction_method": "manual transcription",
            "reported_band_deg": [22.0, 26.0],
            "reported_band_locator": self._source_locator("table-2-range"),
            "observations": [
                {"value_deg": 22.2, "source_locator": self._source_locator("table-2-row-1")},
                {"value_deg": 25.8, "source_locator": self._source_locator("table-2-row-2")},
            ],
        }
        self.assertEqual(reference_audit.load_record(self._write(record))["band_deg"], [22.0, 26.0])

    def test_matched_record_cannot_invent_a_target_band_around_observations(self):
        record = {
            "doi": "10.1000/example", "title": "Example", "first_author_family": "Doe",
            "published_year": 2020, "source_url": "https://doi.org/10.1000/example",
            "comparability": "matched", "sphere_radius_m": 0.002,
            "density_kg_m3": 2500.0, "restitution": 0.4, "sliding_friction": 0.16,
            "container": "cylinder", "deposition": "settled", "release": "lift",
            "angle_definition": "fit", "band_deg": [22.0, 26.0],
            "source_locator": self._source_locator("table-2"), "extraction_method": "manual transcription",
            "reported_band_deg": [23.0, 25.0],
            "reported_band_locator": self._source_locator("table-2-range"),
            "observations": [{"value_deg": 24.0, "source_locator": self._source_locator("table-2-row-1")}],
        }
        with self.assertRaisesRegex(ValueError, "must equal"):
            reference_audit.load_record(self._write(record))

    def test_matched_record_requires_a_locator_for_the_reported_band(self):
        record = {
            "doi": "10.1000/example", "title": "Example", "first_author_family": "Doe",
            "published_year": 2020, "source_url": "https://doi.org/10.1000/example",
            "comparability": "matched", "sphere_radius_m": 0.002,
            "density_kg_m3": 2500.0, "restitution": 0.4, "sliding_friction": 0.16,
            "container": "cylinder", "deposition": "settled", "release": "lift",
            "angle_definition": "fit", "band_deg": [22.0, 26.0],
            "source_locator": self._source_locator("table-2"), "extraction_method": "manual transcription",
            "reported_band_deg": [22.0, 26.0],
            "observations": [{"value_deg": 24.0, "source_locator": self._source_locator("table-2-row-1")}],
        }
        with self.assertRaisesRegex(ValueError, "source-reported band"):
            reference_audit.load_record(self._write(record))

    def test_matched_record_rejects_local_solver_artifact_as_source(self):
        record = {
            "doi": "10.1000/example", "title": "Example", "first_author_family": "Doe",
            "published_year": 2020, "source_url": "https://doi.org/10.1000/example",
            "comparability": "matched", "sphere_radius_m": 0.002,
            "density_kg_m3": 2500.0, "restitution": 0.4, "sliding_friction": 0.16,
            "container": "cylinder", "deposition": "settled", "release": "lift",
            "angle_definition": "fit", "band_deg": [22.0, 26.0],
            "source_locator": "data/repose_sweep.csv", "extraction_method": "manual transcription",
            "reported_band_deg": [22.0, 26.0], "reported_band_locator": "data/repose_sweep.csv",
            "observations": [{"value_deg": 24.0, "source_locator": "data/repose_sweep.csv"}],
        }
        with self.assertRaisesRegex(ValueError, "primary DOI source"):
            reference_audit.load_record(self._write(record))

    def test_crossref_metadata_disagreement_is_rejected(self):
        record = {
            "doi": "10.1000/example", "title": "Wrong title", "first_author_family": "Doe",
            "published_year": 2020, "source_url": "https://doi.org/10.1000/example",
            "comparability": "incompatible",
        }
        def same_doi_opener(request, timeout):
            del request, timeout
            return _Response({"message": {
                "DOI": "10.1000/example", "title": ["Correct title"],
                "author": [{"family": "Doe"}],
                "published-print": {"date-parts": [[2020]]},
            }})

        with self.assertRaisesRegex(ValueError, "Crossref metadata disagrees"):
            reference_audit.audit_record(self._write(record), opener=same_doi_opener)

    def test_record_must_name_its_canonical_doi_resolver(self):
        record = {
            "doi": "10.1000/example", "title": "Example", "first_author_family": "Doe",
            "published_year": 2020, "source_url": "https://example.invalid/paper",
            "comparability": "incompatible",
        }
        with self.assertRaisesRegex(ValueError, "canonical HTTPS DOI resolver"):
            reference_audit.load_record(self._write(record))

    def test_crossref_doi_substitution_is_rejected_even_when_text_matches(self):
        record = {
            "doi": "10.1000/example", "title": "Example", "first_author_family": "Doe",
            "published_year": 2020, "source_url": "https://doi.org/10.1000/example",
            "comparability": "incompatible",
        }

        def wrong_doi_opener(request, timeout):
            del request, timeout
            return _Response({"message": {
                "DOI": "10.1000/substituted", "title": ["Example"],
                "author": [{"family": "Doe"}],
                "published-print": {"date-parts": [[2020]]},
            }})

        with self.assertRaisesRegex(ValueError, "different DOI"):
            reference_audit.audit_record(self._write(record), opener=wrong_doi_opener)


if __name__ == "__main__":
    unittest.main()
