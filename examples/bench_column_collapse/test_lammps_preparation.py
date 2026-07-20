#!/usr/bin/env python3
"""Admission tests for the LAMMPS still-gated preparation witness."""

import pathlib
import sys
import tempfile
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import sweep


class LammpsPreparationTests(unittest.TestCase):
    def write_window(self, directory, rows):
        path = pathlib.Path(directory) / "preparation.txt"
        path.write_text("# Time-averaged data for fix preparation\n" +
                        "\n".join(f"{step} {speed}" for step, speed in rows) + "\n")
        return path

    def test_accepts_the_final_still_gated_window(self):
        with tempfile.TemporaryDirectory() as directory:
            start = sweep.SETTLE_STEPS - (sweep.PREPARATION_WINDOW_SAMPLES - 1) * sweep.PREPARATION_SAMPLE_INTERVAL
            path = self.write_window(directory, [
                (start + i * sweep.PREPARATION_SAMPLE_INTERVAL, 0.0)
                for i in range(sweep.PREPARATION_WINDOW_SAMPLES)
            ])
            self.assertEqual(sweep.lammps_preparation_window(path), [0.0] * sweep.PREPARATION_WINDOW_SAMPLES)

    def test_rejects_a_window_not_ending_at_gate_removal(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.write_window(directory, [
                ((i + 1) * sweep.PREPARATION_SAMPLE_INTERVAL, 0.0)
                for i in range(sweep.PREPARATION_WINDOW_SAMPLES)
            ])
            with self.assertRaisesRegex(ValueError, "invalid LAMMPS preparation"):
                sweep.lammps_preparation_window(path)

    def test_rendered_input_separates_preparation_from_released_arrest(self):
        self.assertIn("fix             preparation", sweep.LMP_TEMPLATE)
        self.assertIn("unfix           preparation", sweep.LMP_TEMPLATE)
        self.assertIn("fix             arrest", sweep.LMP_TEMPLATE)


if __name__ == "__main__":
    unittest.main()
