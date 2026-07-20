#!/usr/bin/env python3
"""Focused tests for the observer's support exclusion and toe geometry."""

import importlib.util
import pathlib
import tempfile
import unittest


MODULE = pathlib.Path(__file__).with_name("independent_observer.py")
SPEC = importlib.util.spec_from_file_location("independent_observer", MODULE)
OBSERVER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(OBSERVER)


class IndependentObserverGeometryTests(unittest.TestCase):
    def test_observer_uses_declared_raw_arrest_cadence(self):
        """A mismatched cadence would make every genuine DIRT witness fail."""
        self.assertEqual(OBSERVER.ARREST_INTERVAL, 100_000)

    def test_frozen_bed_is_not_measured_as_runout(self):
        # The long support would create a false 0.60 m toe if it were selected
        # by z alone.  It must be removed by its release-frame identity.
        release = [
            (0.0, 0.0, 0.003, 0.0015),
            (0.60, 0.0, 0.003, 0.0015),
            (0.095, 0.0, 0.006, 0.0015),
            (0.098, 0.0, 0.006, 0.0015),
        ]
        final = [
            (0.0, 0.0, 0.003, 0.0015),
            (0.60, 0.0, 0.003, 0.0015),
            (0.095, 0.0, 0.006, 0.0015),
            (0.098, 0.0, 0.006, 0.0015),
        ]
        _, mobile = OBSERVER.split_frozen_support(release, final)
        self.assertEqual(len(mobile), 2)
        self.assertAlmostEqual(OBSERVER.interval_toe(mobile), 0.0995)

    def test_toe_requires_the_component_anchored_to_release_footprint(self):
        deposit = [
            (-0.20, 0.0, 0.006, 0.0015),
            (0.095, 0.0, 0.006, 0.0015),
            (0.098, 0.0, 0.006, 0.0015),
        ]
        self.assertAlmostEqual(OBSERVER.interval_toe(deposit), 0.0995)

    def test_moved_support_fails_closed(self):
        release = [(0.0, 0.0, 0.003, 0.0015), (0.01, 0.0, 0.006, 0.0015)]
        final = [(0.00001, 0.0, 0.003, 0.0015), (0.01, 0.0, 0.006, 0.0015)]
        with self.assertRaisesRegex(ValueError, "frozen support"):
            OBSERVER.split_frozen_support(release, final)

    def test_source_population_rejects_underfilled_release(self):
        """Equal release/final counts must not admit an underfilled source."""
        with tempfile.TemporaryDirectory() as directory:
            sweep = pathlib.Path(directory) / "sweep"
            case = sweep / "a0p5"
            case.mkdir(parents=True)
            (sweep / "rough_base.csv").write_text("0,0,0.0015\n")
            (case / "active_column.csv").write_text("0.0015,0.0015,0.0045\n0.0045,0.0015,0.0045\n")
            release = [(0.0, 0.0, 0.0015, 0.0015),
                       (0.0015, 0.0015, 0.0045, 0.0015)]
            frozen = [release[0]]
            active = [release[1]]
            old = OBSERVER.SWEEP
            OBSERVER.SWEEP = str(sweep)
            try:
                with self.assertRaisesRegex(ValueError, "population"):
                    OBSERVER.verify_source_population(0.5, 0, release, active, frozen)
            finally:
                OBSERVER.SWEEP = old

    def test_release_geometry_rejects_grains_past_closed_gate(self):
        active = [
            (0.0015, 0.0015, 0.006, 0.0015),
            (0.1010, 0.0015, 0.006, 0.0015),
        ]
        with self.assertRaisesRegex(ValueError, "still-active gate"):
            OBSERVER.release_aspect(active)


if __name__ == "__main__":
    unittest.main()
