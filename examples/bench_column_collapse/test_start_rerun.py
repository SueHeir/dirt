#!/usr/bin/env python3
"""Regression coverage for explicit column-collapse witness refreshes.

The expensive DIRT executable is replaced with a recorder here: this test
checks the dispatch contract, not a fabricated physical result.  A real case
is only considered admissible after its raw files and receipt pass the driver's
normal checks.
"""

import os
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep


class StartRerunTest(unittest.TestCase):
    def test_rerun_dispatches_an_admitted_case(self):
        """--rerun must execute, whereas normal resume may reuse evidence."""
        case = (sweep.ASPECTS[0], sweep.SEEDS[0])
        calls = []

        def record_run(a, seed, binary, env):
            calls.append((a, seed))
            return a, seed

        with tempfile.TemporaryDirectory() as directory, \
             mock.patch.object(sweep, "DATA_DIR", directory), \
             mock.patch.object(sweep, "checked_protocol_manifest", return_value={}), \
             mock.patch.object(sweep.subprocess, "run"), \
             mock.patch.object(sweep.os.path, "isfile", return_value=True), \
             mock.patch.object(sweep, "_case_evidence_error", return_value=None), \
             mock.patch.object(sweep, "_run_dirt_case", side_effect=record_run):
            sweep.start(jobs=1, rerun=True, selected_cases=[case])

        self.assertEqual(calls, [case])

    def test_case_lock_rejects_a_duplicate_writer_but_not_another_case(self):
        """Only duplicate witnesses conflict; independent scheduler work proceeds."""
        with tempfile.TemporaryDirectory() as directory, \
             mock.patch.object(sweep, "SWEEP_DIR", directory):
            with sweep.exclusive_case_lock(0.5, 0, "written"):
                with self.assertRaisesRegex(ValueError, "a=0.5 seed=0"):
                    with sweep.exclusive_case_lock(0.5, 0, "written"):
                        pass
                with self.assertRaisesRegex(ValueError, "cannot derive"):
                    with sweep.shared_case_lock(0.5, 0, "derive the ensemble"):
                        pass
                with sweep.exclusive_case_lock(0.5, 1, "written"):
                    pass


if __name__ == "__main__":
    unittest.main()
