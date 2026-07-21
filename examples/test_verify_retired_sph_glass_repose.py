#!/usr/bin/env python3
"""Regression tests for the retired-SPH scope boundary, not physics validation."""

from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("verify_retired_sph_glass_repose.py")
SPEC = importlib.util.spec_from_file_location("retired_sph_audit", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


class RetiredSphScopeTests(unittest.TestCase):
    def test_exact_historical_blob_at_a_new_path_is_rejected(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "relocated"):
            audit.require_no_relocated_historical_blobs(
                {"retired/main.rs": "abc123"}, {"examples/new_name/main.rs": "abc123"}
            )

    def test_unrelated_new_implementation_is_not_called_a_relocation(self) -> None:
        audit.require_no_relocated_historical_blobs(
            {"retired/main.rs": "abc123"}, {"examples/new_name/main.rs": "def456"}
        )


if __name__ == "__main__":
    unittest.main()
