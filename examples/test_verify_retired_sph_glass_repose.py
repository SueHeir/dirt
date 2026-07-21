#!/usr/bin/env python3
"""Unit tests for the retired-SPH evidence boundary.

These tests exercise negative cases only.  They do not substitute local mock
metadata for the live Crossref/OpenAlex audit required by the actual command.
"""

import importlib.util
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("verify_retired_sph_glass_repose.py")
SPEC = importlib.util.spec_from_file_location("retired_repose_audit", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


class RetiredSurfaceTests(unittest.TestCase):
    def test_complete_removed_surface_is_accepted(self) -> None:
        historical = {"case/README.md", "case/main.rs", "case/sweep.py"}
        audit.require_retired_surface_absent(historical, historical, set())

    def test_partial_removal_is_rejected(self) -> None:
        historical = {"case/README.md", "case/main.rs", "case/sweep.py"}
        with self.assertRaisesRegex(RuntimeError, "leaves historical"):
            audit.require_retired_surface_absent(historical, historical - {"case/sweep.py"}, set())

    def test_restored_non_entrypoint_is_rejected(self) -> None:
        historical = {"case/README.md", "case/main.rs", "case/pin_mu0.3_s1.toml"}
        with self.assertRaisesRegex(RuntimeError, "restored"):
            audit.require_retired_surface_absent(historical, historical, {"case/pin_mu0.3_s1.toml"})

    def test_empty_predecessor_is_rejected(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "contains no"):
            audit.require_retired_surface_absent(set(), set(), set())

    def test_replacement_candidate_in_manifest_is_reported(self) -> None:
        self.assertEqual(
            audit.replacement_candidates(
                {"examples/angle_of_repose/main.rs", "src/constitutive.rs"},
                {"src/constitutive.rs": "generic friction"},
            ),
            ["examples/angle_of_repose/main.rs"],
        )

    def test_generic_friction_is_not_mislabelled_as_replacement(self) -> None:
        self.assertEqual(
            audit.replacement_candidates(
                {"src/constitutive.rs"}, {"src/constitutive.rs": "rolling friction coefficient"}
            ),
            [],
        )

    def test_documentation_mention_is_not_a_runnable_replacement(self) -> None:
        self.assertEqual(
            audit.replacement_candidates(
                {"docs/angle_of_repose.md"}, {"docs/angle_of_repose.md": "angle of repose"}
            ),
            [],
        )


if __name__ == "__main__":
    unittest.main()
