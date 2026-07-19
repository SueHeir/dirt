#!/usr/bin/env python3
"""Adversarial admission tests for the optional LAMMPS witness."""

import csv
import os
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep


class LammpsReleaseSupportTests(unittest.TestCase):
    BASE = [(0.0015, 0.0015, 0.0015), (0.0045, 0.0015, 0.0015)]

    def write_release(self, directory, points):
        path = os.path.join(directory, "release.csv")
        with open(path, "w", newline="") as f:
            out = csv.DictWriter(f, fieldnames=["x", "y", "z", "radius"])
            out.writeheader()
            for x, y, z in points:
                out.writerow({"x": x, "y": y, "z": z, "radius": 0.0015})
        return path

    def test_accepts_text_roundtrip_of_every_frozen_coordinate(self):
        with tempfile.TemporaryDirectory() as directory, \
             mock.patch.object(sweep, "rough_base_positions", return_value=self.BASE):
            path = self.write_release(directory, [
                (0.001500000001, 0.0015, 0.0015),
                (0.0045, 0.001499999999, 0.0015),
                (0.010, 0.0015, 0.0045),
            ])
            sweep.checked_lammps_release_support(path)

    def test_rejects_same_population_when_one_support_coordinate_moved(self):
        with tempfile.TemporaryDirectory() as directory, \
             mock.patch.object(sweep, "rough_base_positions", return_value=self.BASE):
            # The third row preserves the total population but cannot replace
            # the immutable base point at x=4.5 mm.
            path = self.write_release(directory, [
                self.BASE[0], (0.010, 0.0015, 0.0045), (0.012, 0.0015, 0.0045),
            ])
            with self.assertRaisesRegex(ValueError, "does not retain"):
                sweep.checked_lammps_release_support(path)


if __name__ == "__main__":
    unittest.main()
