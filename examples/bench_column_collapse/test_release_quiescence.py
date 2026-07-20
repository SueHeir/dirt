#!/usr/bin/env python3
"""Adversarial admission tests for the pre-release kinetic witness."""

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
    def state(self, directory, speed):
        path = pathlib.Path(directory) / "release_state.csv"
        with path.open("w", newline="") as source:
            writer = csv.DictWriter(source, fieldnames=["particle_count", "max_speed_m_s"])
            writer.writeheader()
            writer.writerow({"particle_count": 7, "max_speed_m_s": speed})
        return path

    def test_driver_rejects_moving_pre_release_frame(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.state(directory, 0.05)
            vmax = SWEEP.checked_release_state(path, 7)
            with self.assertRaisesRegex(ValueError, "pre-release Fr"):
                SWEEP.checked_release_quiescence(vmax)

    def test_observer_independently_rejects_moving_pre_release_frame(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.state(directory, 0.05)
            with self.assertRaisesRegex(ValueError, "pre-release kinetic"):
                OBSERVER.released_at_rest(path, 7)

    def test_recorder_writes_release_state_before_gate_deactivation(self):
        source = (HERE / "main.rs").read_text()
        self.assertLess(source.index("column_collapse_release_state.csv"),
                        source.index("walls.deactivate_by_name(GATE_NAME)"))


if __name__ == "__main__":
    unittest.main()
