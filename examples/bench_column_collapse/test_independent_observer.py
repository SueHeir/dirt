#!/usr/bin/env python3
"""Focused tests for the observer's support exclusion and toe geometry."""

import importlib.util
import pathlib
import unittest


MODULE = pathlib.Path(__file__).with_name("independent_observer.py")
SPEC = importlib.util.spec_from_file_location("independent_observer", MODULE)
OBSERVER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(OBSERVER)


class IndependentObserverGeometryTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
