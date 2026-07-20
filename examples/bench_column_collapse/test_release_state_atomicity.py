#!/usr/bin/env python3
"""The pre-release kinetic witness must share the case-receipt lifecycle."""

import os
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep


class ReleaseStateAtomicityTests(unittest.TestCase):
    def test_case_paths_include_the_release_kinetic_witness(self):
        """Ensemble derivation and receipt hashing use one canonical path map."""
        paths = sweep.dirt_case_paths(sweep.ASPECTS[0], sweep.SEEDS[0])
        self.assertIn("release_state", paths)
        self.assertTrue(paths["release_state"].endswith("column_collapse_release_state.csv"))

    def test_rerun_invalidates_release_state_with_the_other_witnesses(self):
        """A fresh trajectory must not inherit an old pre-release rest sample."""
        with tempfile.TemporaryDirectory() as directory, \
             mock.patch.object(sweep, "case_dir_seed", return_value=directory):
            data = os.path.join(directory, "data")
            os.mkdir(data)
            release_state = os.path.join(data, "column_collapse_release_state.csv")
            untouched = os.path.join(data, "unrelated.csv")
            open(release_state, "w").close()
            open(untouched, "w").close()
            sweep._clear_case_evidence(sweep.ASPECTS[0], sweep.SEEDS[0])
            self.assertFalse(os.path.exists(release_state))
            self.assertTrue(os.path.exists(untouched))


if __name__ == "__main__":
    unittest.main()
