#!/usr/bin/env python3
"""Unit tests for source-derived parsing; live catalogues are not mocked."""

import importlib.util
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("audit_retired_sph_claim.py")
SPEC = importlib.util.spec_from_file_location("retired_claim_audit", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


class HistoricalClaimTests(unittest.TestCase):
    def test_archived_claim_and_quoted_reference_are_derived(self) -> None:
        band, title = audit.archived_claim_and_title(
            "measured glass repose band [22, 26]\n## References\n"
            '1. Author, "Rolling friction in the dynamic simulation of sandpile formation"'
        )
        self.assertEqual(band, "[22, 26]")
        self.assertEqual(title, "rolling friction in the dynamic simulation of sandpile formation")

    def test_complete_bibliography_is_derived_in_source_order(self) -> None:
        band, titles = audit.archived_claim_and_references(
            "measured glass repose band [22, 26]\n## References\n"
            '1. A, "First source"\n2. B, "Second source"\n3. C, "Third source"'
        )
        self.assertEqual(band, "[22, 26]")
        self.assertEqual(titles, ["first source", "second source", "third source"])

    def test_missing_historical_claim_is_rejected(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "claim not found"):
            audit.archived_claim_and_title("## References\nComputer simulation of sandpile formation")

    def test_claim_without_a_source_band_is_rejected(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "no numeric band"):
            audit.archived_claim_and_title(
                'measured glass repose band\n## References\n1. Author, "A simulation study"'
            )

    def test_source_band_is_not_a_local_expected_value(self) -> None:
        band, _ = audit.archived_claim_and_title(
            'measured glass repose band [19.5, 27]\n## References\n1. Author, "A simulation study"'
        )
        self.assertEqual(band, "[19.5, 27]")

    def test_band_is_taken_from_the_claim_line_not_an_earlier_number(self) -> None:
        band, _ = audit.archived_claim_and_references(
            'an unrelated range [1, 2]\nmeasured glass repose band [19.5, 27]\n## References\n1. A, "A source"'
        )
        self.assertEqual(band, "[19.5, 27]")

    def test_missing_reference_section_is_rejected(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "no reference section"):
            audit.archived_claim_and_title("measured glass repose band [22, 26]")

    def test_unquoted_first_reference_is_rejected(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "no quoted title"):
            audit.archived_claim_and_title(
                "measured glass repose band [22, 26]\n## References\n1. An unquoted reference"
            )

    def test_inline_numeric_citation_requires_primary_source_review(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "inline citation"):
            audit.archived_claim_and_references(
                'measured glass repose band [22, 26] [1]\n## References\n1. A, "A source"'
            )

    def test_missing_bibliography_number_is_rejected(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "numbering"):
            audit.archived_claim_and_references(
                'measured glass repose band [22, 26]\n## References\n1. A, "First"\n3. B, "Third"'
            )


if __name__ == "__main__":
    unittest.main()
