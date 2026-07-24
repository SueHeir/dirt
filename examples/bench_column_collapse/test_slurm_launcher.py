#!/usr/bin/env python3
"""Regression coverage for the batch-worker execution boundary.

The launcher is intentionally a shell program so Slurm can execute one
immutable manifest entry directly.  These checks pin the safety properties
that matter before a multi-hour trajectory starts: syntax stays valid and the
worker must load an explicit reproducible build environment rather than rely on
the submit host's interactive shell state.
"""

import pathlib
import subprocess
import unittest


LAUNCHER = pathlib.Path(__file__).with_name("run_slurm_array.sh")


class SlurmLauncherTests(unittest.TestCase):
    def test_launcher_has_valid_bash_syntax(self):
        subprocess.run(["bash", "-n", str(LAUNCHER)], check=True)

    def test_launcher_requires_a_declared_build_environment(self):
        text = LAUNCHER.read_text()
        self.assertIn('DIRT_BUILD_ENV:-"$HOME/projects/.build-env"', text)
        self.assertIn('source "$build_env"', text)
        self.assertIn('missing DIRT build environment', text)


if __name__ == "__main__":
    unittest.main()
