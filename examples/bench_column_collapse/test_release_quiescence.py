#!/usr/bin/env python3
"""Adversarial admission tests for the still-gated preparation witness."""

import csv
import importlib.util
import pathlib
import tempfile
import unittest


HERE = pathlib.Path(__file__).parent


def load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SWEEP = load("sweep")
OBSERVER = load("independent_observer")


class ReleaseQuiescenceTests(unittest.TestCase):
    def preparation(self, directory, speed, count=7, final_step=None):
        path = pathlib.Path(directory) / "preparation.csv"
        with path.open("w", newline="") as source:
            writer = csv.DictWriter(source, fieldnames=["settle_step", "particle_count", "max_speed_m_s"])
            writer.writeheader()
            final_step = SWEEP.SETTLE_STEPS if final_step is None else final_step
            for step in range(final_step - 3 * SWEEP.PREPARATION_SAMPLE_INTERVAL,
                              final_step + 1, SWEEP.PREPARATION_SAMPLE_INTERVAL):
                writer.writerow({"settle_step": step, "particle_count": count, "max_speed_m_s": speed})
        return path

    def test_driver_rejects_moving_still_gated_window(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "preparation-rest Fr"):
                SWEEP.checked_preparation_window(self.preparation(directory, 0.05), 7)

    def test_observer_independently_rejects_moving_still_gated_window(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "preparation-rest"):
                OBSERVER.prepared_at_rest(self.preparation(directory, 0.05), 7)

    def test_final_still_gated_sample_precedes_gate_deactivation(self):
        source = (HERE / "main.rs").read_text()
        self.assertIn("record_preparation_quiescence.run_if(in_stage(\"settle\"))", source)
        self.assertIn("begin_collapse.run_if(in_stage(\"collapse\"))", source)

    def test_driver_rejects_a_preparation_window_that_does_not_end_at_release(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.preparation(directory, 0.0, final_step=SWEEP.SETTLE_STEPS - SWEEP.PREPARATION_SAMPLE_INTERVAL)
            with self.assertRaisesRegex(ValueError, "invalid preparation-rest"):
                SWEEP.checked_preparation_window(path, 7)

    def test_driver_rejects_preparation_population_drift(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "invalid preparation-rest"):
                SWEEP.checked_preparation_window(self.preparation(directory, 0.0, count=6), 7)


if __name__ == "__main__":
    unittest.main()
