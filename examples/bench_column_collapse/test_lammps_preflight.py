#!/usr/bin/env python3
"""Execution-contract coverage for the optional LAMMPS overlay preflight."""

import os
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep


class LammpsPreflightTests(unittest.TestCase):
    def test_preflight_runs_both_rendered_stages_at_zero_steps(self):
        with tempfile.TemporaryDirectory() as directory, \
             mock.patch.object(sweep.tempfile, "TemporaryDirectory") as temporary, \
             mock.patch.object(sweep, "write_lammps_input") as render, \
             mock.patch.object(sweep.subprocess, "run", return_value=mock.Mock(returncode=0)) as run:
            temporary.return_value.__enter__.return_value = directory
            self.assertTrue(sweep.lammps_preflight("lmp"))

        render.assert_called_once()
        self.assertEqual(render.call_args.kwargs["settle_steps"], 0)
        self.assertEqual(render.call_args.kwargs["collapse_steps"], 0)
        self.assertEqual(render.call_args.kwargs["output_dir"], directory)
        self.assertEqual(run.call_args.args[0][:2], ["lmp", "-in"])

    def test_preflight_refuses_to_start_overlay_after_lammps_failure(self):
        with tempfile.TemporaryDirectory() as directory, \
             mock.patch.object(sweep.tempfile, "TemporaryDirectory") as temporary, \
             mock.patch.object(sweep, "write_lammps_input"), \
             mock.patch.object(sweep.subprocess, "run", return_value=mock.Mock(returncode=1)):
            temporary.return_value.__enter__.return_value = directory
            self.assertFalse(sweep.lammps_preflight("lmp"))


if __name__ == "__main__":
    unittest.main()
