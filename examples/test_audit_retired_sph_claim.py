#!/usr/bin/env python3
"""Negative tests for the archived-claim audit; live catalogues are not mocked."""

import importlib.util
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("audit_retired_sph_claim.py")
SPEC = importlib.util.spec_from_file_location("retired_claim_audit", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


class HistoricalClaimTests(unittest.TestCase):
    def test_historical_claim_with_simulation_only_reference_is_auditable(self) -> None:
        audit.require_unsupported_claim(
            "measured glass repose band [22, 26]\n## References\n"
            "Rolling friction in the dynamic simulation of sandpile formation"
        )

    def test_missing_historical_claim_is_rejected(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "claim not found"):
            audit.require_unsupported_claim("## References\nComputer simulation of sandpile formation")

    def test_missing_reference_section_is_rejected(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "no reference section"):
            audit.require_unsupported_claim("measured glass repose band [22, 26]")

    def test_wrong_historical_paper_is_rejected(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "historical cited paper not found"):
            audit.require_unsupported_claim(
                "measured glass repose band [22, 26]\n## References\n"
                "A different simulation of sandpile formation"
            )

    def test_experimental_label_in_reference_section_is_rejected(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "unexpectedly claims an experiment"):
            audit.require_unsupported_claim(
                "measured glass repose band [22, 26]\n## References\n"
                "Rolling friction in the dynamic simulation of sandpile formation experiment"
            )

    def test_archived_author_attribution_mismatch_is_detected(self) -> None:
        # Crossref/OpenAlex identify the simulation paper as Zhou, Wright,
        # Yang, Xu, Yu; the archived README instead listed Zulli.
        self.assertFalse(audit.archived_attribution_matches({"zhou", "wright", "yang", "xu", "yu"}))


if __name__ == "__main__":
    unittest.main()
