#!/usr/bin/env python3
"""Fail-closed coverage for the scheduler-to-witness boundary."""

import csv
import os
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep


class JobManifestValidationTests(unittest.TestCase):
    def _rows(self):
        cases = []
        rows = []
        for index, (aspect, seed) in enumerate(
                ((a, s) for a in sweep.ASPECTS for s in sweep.SEEDS), start=1):
            cases.append({"nominal_aspect": aspect, "seed": seed,
                          "active_count": index, "active_source_sha256": f"source-{index}"})
            rows.append({"index": str(index), "aspect": f"{aspect:g}", "seed": str(seed),
                         "active_count": str(index), "active_source_sha256": f"source-{index}",
                         "protocol_sha256": "protocol",
                         "command": f"python3 examples/bench_column_collapse/sweep.py start --case {aspect:g},{seed}"})
        return cases, rows

    def test_rejects_duplicate_case_even_when_row_count_is_33(self):
        cases, rows = self._rows()
        rows[-1]["aspect"] = rows[-2]["aspect"]
        rows[-1]["seed"] = rows[-2]["seed"]
        with tempfile.TemporaryDirectory() as directory, \
             mock.patch.object(sweep, "SWEEP_DIR", directory), \
             mock.patch.object(sweep, "checked_protocol_manifest",
                               return_value={"protocol_sha256": "protocol", "cases": cases}):
            path = os.path.join(directory, sweep.JOB_MANIFEST_NAME)
            with open(path, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaisesRegex(ValueError, "disagrees"):
                sweep.validate_job_manifest()


if __name__ == "__main__":
    unittest.main()
