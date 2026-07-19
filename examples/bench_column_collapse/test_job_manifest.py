#!/usr/bin/env python3
"""Regression coverage for distributed column-collapse witness dispatch."""

import csv
import os
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep


class JobManifestTest(unittest.TestCase):
    def test_manifest_has_one_digest_bound_command_per_declared_case(self):
        """A scheduler map must cover exactly the required 11 x 3 witnesses."""
        cases = []
        for aspect in sweep.ASPECTS:
            for seed in sweep.SEEDS:
                cases.append({
                    "nominal_aspect": aspect,
                    "seed": seed,
                    "active_count": 100 + seed,
                    "active_source_sha256": f"source-{aspect}-{seed}",
                })
        manifest = {"protocol_sha256": "protocol-digest", "cases": cases}
        with tempfile.TemporaryDirectory() as directory, \
             mock.patch.object(sweep, "SWEEP_DIR", directory), \
             mock.patch.object(sweep, "checked_protocol_manifest", return_value=manifest):
            path = sweep.emit_jobs()
            with open(path, newline="") as f:
                rows = list(csv.DictReader(f, delimiter="\t"))

        self.assertEqual(len(rows), len(sweep.ASPECTS) * len(sweep.SEEDS))
        self.assertEqual([int(row["index"]) for row in rows], list(range(1, 34)))
        self.assertEqual({row["protocol_sha256"] for row in rows}, {"protocol-digest"})
        self.assertEqual(len({(row["aspect"], row["seed"]) for row in rows}), 33)
        for row in rows:
            self.assertIn(f"--case {row['aspect']},{row['seed']}", row["command"])
            self.assertTrue(row["active_source_sha256"].startswith("source-"))


if __name__ == "__main__":
    unittest.main()
