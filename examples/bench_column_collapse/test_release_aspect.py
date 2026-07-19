#!/usr/bin/env python3
"""Regression tests for measured-aspect admission in both analyses."""

import importlib.util
import pathlib
import unittest


HERE = pathlib.Path(__file__).parent


def load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SWEEP = load("sweep")
OBSERVER = load("independent_observer")


class ReleaseAspectAdmissionTests(unittest.TestCase):
    def test_driver_accepts_exact_aspect_and_rejects_off_schedule_release(self):
        self.assertEqual(SWEEP.checked_release_aspect(2.0 * SWEEP.L0, 2.0), 2.0)
        with self.assertRaisesRegex(ValueError, "differs from scheduled"):
            SWEEP.checked_release_aspect(2.05 * SWEEP.L0, 2.0)

    def test_independent_observer_has_the_same_strict_aspect_boundary(self):
        self.assertEqual(OBSERVER.checked_release_aspect(3.0, 3.0), 3.0)
        with self.assertRaisesRegex(ValueError, "differs from scheduled"):
            OBSERVER.checked_release_aspect(3.07, 3.0)

    def test_driver_rejects_a_height_correct_but_compacted_release_width(self):
        """The external law uses the released H_i/L_i, not the gate location."""
        with self.assertRaisesRegex(ValueError, "release width"):
            SWEEP.checked_release_dimensions(2.0 * SWEEP.L0, 0.95 * SWEEP.L0, 2.0)

    def test_observer_normalizes_to_the_measured_release_front(self):
        active = [
            (0.0015, 0.0015, 0.006, 0.0015),
            (0.0945, 0.0015, 0.006, 0.0015),
        ]
        height, width, left, right = OBSERVER.release_dimensions(active)
        self.assertAlmostEqual(width, OBSERVER.L0)
        self.assertAlmostEqual(right - left, OBSERVER.L0)
        self.assertAlmostEqual(height / width, OBSERVER.release_aspect(active))


if __name__ == "__main__":
    unittest.main()
