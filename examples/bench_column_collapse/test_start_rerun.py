#!/usr/bin/env python3
"""Regression coverage for explicit column-collapse witness refreshes.

The expensive DIRT executable is replaced with a recorder here: this test
checks the dispatch contract, not a fabricated physical result.  A real case
is only considered admissible after its raw files and receipt pass the driver's
normal checks.
"""

import contextlib
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
             mock.patch.object(sweep.os, "remove"), \
             mock.patch.object(sweep, "_case_evidence_error", return_value=None), \
             mock.patch.object(sweep, "_run_dirt_case", side_effect=record_run):
            sweep.start(jobs=1, rerun=True, selected_cases=[case])

        self.assertEqual(calls, [case])

    def test_lammps_is_not_implicitly_launched_after_a_dirt_campaign(self):
        """An installed external binary must not make an optional overlay mandatory."""
        rows = [{"nominal_aspect": aspect} for aspect in sweep.ASPECTS]
        cases = [(aspect, seed) for aspect in sweep.ASPECTS for seed in sweep.SEEDS]
        with tempfile.TemporaryDirectory() as directory, \
             mock.patch.object(sweep, "DATA_DIR", directory), \
             mock.patch.object(sweep, "checked_protocol_manifest", return_value={}), \
             mock.patch.object(sweep.subprocess, "run"), \
             mock.patch.object(sweep.os.path, "isfile", return_value=True), \
             mock.patch.object(sweep, "_case_evidence_error", return_value=None), \
             mock.patch.object(sweep, "derive_dirt_ensemble", return_value=rows), \
             mock.patch.object(sweep, "_write_runout"), \
             mock.patch.object(sweep, "find_lammps", return_value="lmp") as find_lammps, \
             mock.patch.object(sweep, "run_lammps_sweep") as run_lammps:
            sweep._start(jobs=1, selected_cases=cases)

        find_lammps.assert_not_called()
        run_lammps.assert_not_called()

    def test_lammps_overlay_requires_explicit_opt_in(self):
        """The external dynamics leg remains available when deliberately requested."""
        rows = [{"nominal_aspect": aspect} for aspect in sweep.ASPECTS]
        cases = [(aspect, seed) for aspect in sweep.ASPECTS for seed in sweep.SEEDS]
        with tempfile.TemporaryDirectory() as directory, \
             mock.patch.object(sweep, "DATA_DIR", directory), \
             mock.patch.object(sweep, "checked_protocol_manifest", return_value={}), \
             mock.patch.object(sweep.subprocess, "run"), \
             mock.patch.object(sweep.os.path, "isfile", return_value=True), \
             mock.patch.object(sweep.os, "remove"), \
             mock.patch.object(sweep, "_case_evidence_error", return_value=None), \
             mock.patch.object(sweep, "derive_dirt_ensemble", return_value=rows), \
             mock.patch.object(sweep, "_write_runout"), \
             mock.patch.object(sweep, "find_lammps", return_value="lmp"), \
             mock.patch.object(sweep, "run_lammps_sweep", return_value=rows) as run_lammps:
            sweep._start(jobs=1, selected_cases=cases, with_lammps=True)

        run_lammps.assert_called_once_with("lmp")

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

    def test_rerun_clears_preparation_witness_with_the_other_raw_evidence(self):
        """A new trajectory must not inherit an old quiet-preparation tail."""
        case = (sweep.ASPECTS[0], sweep.SEEDS[0])
        with tempfile.TemporaryDirectory() as directory, \
             mock.patch.object(sweep, "SWEEP_DIR", directory):
            data = os.path.join(sweep.case_dir_seed(*case), "data")
            os.makedirs(data)
            names = (
                "column_collapse_results.csv",
                "column_collapse_release.csv",
                "column_collapse_final_state.csv",
                "column_collapse_arrest.csv",
                "column_collapse_preparation.csv",
                sweep.CASE_RECEIPT_NAME,
            )
            for name in names:
                with open(os.path.join(data, name), "w") as witness:
                    witness.write("stale\n")

            sweep._clear_case_evidence(*case)

            self.assertEqual(os.listdir(data), [])

    def test_worker_rejects_a_terminated_but_unqualified_witness(self):
        """Scheduler success must mean physical admission, not only exit zero."""
        case = (sweep.ASPECTS[0], sweep.SEEDS[0])
        with tempfile.TemporaryDirectory() as directory, \
             mock.patch.object(sweep, "SWEEP_DIR", directory), \
             mock.patch.object(sweep, "exclusive_case_lock", lambda *args: contextlib.nullcontext()), \
             mock.patch.object(sweep.subprocess, "run"), \
             mock.patch.object(sweep, "write_case_receipt"), \
             mock.patch.object(sweep, "_case_evidence_error", return_value="terminal Fr=0.2 > 0.01"):
            cdir = sweep.case_dir_seed(*case)
            os.makedirs(cdir)
            with open(os.path.join(cdir, "config.toml"), "w"):
                pass
            with self.assertRaisesRegex(RuntimeError, "not admitted: terminal Fr"):
                sweep._run_dirt_case(*case, binary="dirt", env={})


if __name__ == "__main__":
    unittest.main()
