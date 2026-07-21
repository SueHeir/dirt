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

    def test_missing_historical_claim_is_rejected(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "claim not found"):
            audit.archived_claim_and_title("## References\nComputer simulation of sandpile formation")

    def test_missing_reference_section_is_rejected(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "no reference section"):
            audit.archived_claim_and_title("measured glass repose band [22, 26]")

    def test_unquoted_first_reference_is_rejected(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "no quoted title"):
            audit.archived_claim_and_title(
                "measured glass repose band [22, 26]\n## References\n1. An unquoted reference"
            )


if __name__ == "__main__":
    unittest.main()
