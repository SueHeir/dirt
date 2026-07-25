#!/usr/bin/env python3
"""Keep the LAMMPS cache on the same measured-geometry contract as its audit."""

import os
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep


class LammpsSummaryGeometryTests(unittest.TestCase):
    def test_cached_overlay_uses_each_witness_release_front_and_width(self):
        """A gate coordinate must not replace the measured released footprint."""
        aspect = sweep.ASPECTS[0]
        with tempfile.TemporaryDirectory() as directory, \
             mock.patch.object(sweep, "ASPECTS", [aspect]), \
             mock.patch.object(sweep, "SEEDS", [0]), \
             mock.patch.object(sweep, "case_dir_seed", return_value=directory), \
             mock.patch.object(sweep, "n_particles", return_value=1), \
             mock.patch.object(sweep, "total_particles", return_value=1), \
             mock.patch.object(sweep, "write_lammps_input"), \
             mock.patch.object(sweep, "lammps_preflight", return_value=True), \
             mock.patch.object(sweep, "lammps_dump_to_csv", return_value=True), \
             mock.patch.object(sweep, "csv_particle_count", return_value=1), \
             mock.patch.object(sweep, "checked_lammps_release_support"), \
             mock.patch.object(sweep, "release_geometry", return_value=(0.050, 0.094, 0.0, 0.094)), \
             mock.patch.object(sweep, "checked_release_dimensions", return_value=0.050 / 0.094), \
             mock.patch.object(sweep, "lammps_preparation_window", return_value=[0.0]), \
             mock.patch.object(sweep, "lammps_arrest_window", return_value=[0.0]), \
             mock.patch.object(sweep, "lammps_max_speed", return_value=0.0), \
             mock.patch.object(sweep, "measure_column", return_value=(0.0, 0.188)), \
             mock.patch.object(sweep, "write_lammps_case_receipt"), \
             mock.patch.object(sweep.os.path, "isfile", return_value=True), \
             mock.patch.object(sweep.os, "remove"), \
             mock.patch.object(sweep.subprocess, "run", return_value=mock.Mock(returncode=0)):
            rows = sweep.run_lammps_sweep("lmp")

        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertAlmostEqual(row["L0"], 0.094)
        self.assertAlmostEqual(row["release_front"], 0.094)
        self.assertAlmostEqual(row["aspect"], 0.050 / 0.094)
        self.assertAlmostEqual(row["runout_norm"], (0.188 - 0.094) / 0.094)


if __name__ == "__main__":
    unittest.main()
