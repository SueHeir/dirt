#!/usr/bin/env python3
"""Keep retired-SPH provenance checks out of a green validation summary.

This is a CI-semantics guard, not a scientific validation test.  Its purpose is
to prevent a successful absence/citation audit from being reported beside DEM
benchmark PASS lines and then misread as calibration evidence.
"""

import importlib.util
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("ci_validation.py")
SPEC = importlib.util.spec_from_file_location("ci_validation", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
ci = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ci)


class ValidationScopeBoundaryTests(unittest.TestCase):
    def test_retired_sph_audits_are_not_green_validation_sweeps(self) -> None:
        prohibited = {
            "examples/verify_retired_sph_glass_repose.py",
            "examples/audit_retired_sph_claim.py",
        }
        configured = ci.sweep_paths(ci.SMOKE_SWEEPS) | ci.sweep_paths(ci.FULL_SWEEPS)
        self.assertTrue(prohibited.isdisjoint(configured))


if __name__ == "__main__":
    unittest.main()
