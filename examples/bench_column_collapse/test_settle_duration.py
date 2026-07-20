#!/usr/bin/env python3
"""Keep the release-rest preparation interval identical in both solvers."""

import importlib.util
import pathlib
import unittest


HERE = pathlib.Path(__file__).parent
spec = importlib.util.spec_from_file_location("collapse_sweep", HERE / "sweep.py")
SWEEP = importlib.util.module_from_spec(spec)
spec.loader.exec_module(SWEEP)


class SettleDurationTests(unittest.TestCase):
    def test_preparation_is_longer_than_the_old_transient(self):
        # A release witness is an admission condition, so do not silently
        # regress to the 0.32 s insertion-transient-only preparation.
        self.assertEqual(SWEEP.SETTLE_STEPS, 800_000)
        self.assertAlmostEqual(SWEEP.SETTLE_STEPS * SWEEP.SETTLE_DT, 0.8)

    def test_release_preserves_four_seconds_at_the_resolved_timestep(self):
        self.assertEqual(SWEEP.COLLAPSE_STEPS, 4_000_000)
        self.assertAlmostEqual(SWEEP.COLLAPSE_STEPS * SWEEP.COLLAPSE_DT, 4.0)
        self.assertEqual(SWEEP.COLLAPSE_DT, SWEEP.SETTLE_DT)

    def test_lammps_uses_the_same_settle_duration(self):
        # Rendering receives the same ``settle_steps`` argument in
        # ``write_lammps_input``; retain the placeholder so a future deck
        # cannot drift to a separately hard-coded preparation length.
        self.assertIn("run             {settle_steps}", SWEEP.LMP_TEMPLATE)
        self.assertIn("timestep        {settle_dt}", SWEEP.LMP_TEMPLATE)
        self.assertIn("timestep        {collapse_dt}", SWEEP.LMP_TEMPLATE)


if __name__ == "__main__":
    unittest.main()
