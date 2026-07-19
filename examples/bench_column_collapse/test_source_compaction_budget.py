#!/usr/bin/env python3
"""Regression tests for a source that can satisfy the unchanged release gate."""

import importlib.util
import pathlib
import unittest


ROOT = pathlib.Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("collapse_sweep", ROOT / "sweep.py")
SWEEP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SWEEP)


class SourceCompactionBudgetTests(unittest.TestCase):
    def test_source_height_compensates_deliberate_fcc_clearance(self):
        # This validates only the prescribed initial geometry; every dynamic
        # witness still has to pass its measured 2% release-aspect gate.
        for aspect in (0.5, 3.0, 5.0):
            target = aspect * SWEEP.L0 * SWEEP.SOURCE_DILATION
            layer = (2.0 / 3.0) ** 0.5 * SWEEP.SOURCE_DILATION * 2.0 * SWEEP.RADIUS
            for seed in SWEEP.SEEDS:
                points = SWEEP.active_column_positions(SWEEP.n_particles(aspect), aspect, seed)
                self.assertLessEqual(abs(SWEEP.source_preparation_height(points) - target), layer)

    def test_all_boundary_sources_remain_exact_and_nonoverlapping(self):
        for aspect in (0.5, 3.0, 5.0):
            count = SWEEP.n_particles(aspect)
            for seed in SWEEP.SEEDS:
                points = SWEEP.active_column_positions(count, aspect, seed)
                SWEEP.audit_active_source(points, count, aspect)
                self.assertGreaterEqual(
                    SWEEP.source_min_separation(points),
                    2.0 * SWEEP.RADIUS * (1.0 - 1.0e-10),
                )


if __name__ == "__main__":
    unittest.main()
