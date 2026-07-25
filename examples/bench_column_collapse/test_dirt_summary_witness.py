#!/usr/bin/env python3
"""The aggregate DIRT CSV must never be the evidence used by ``graph``."""

import csv
import os
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep


def witnessed_rows():
    return [{
        "nominal_aspect": aspect,
        "aspect": aspect,
        "L0": 0.096,
        "H": 0.096 * aspect,
        "L_f": 0.096 * (1.0 + aspect),
        "release_front": 0.096,
        "runout_norm": aspect,
        "runout_std": 0.0,
        "n_seeds": len(sweep.SEEDS),
        "protocol_sha256": "raw-witness-contract",
    } for aspect in sweep.ASPECTS]


class DirtSummaryWitnessTests(unittest.TestCase):
    def test_graph_refuses_an_edited_summary_even_when_it_has_all_aspects(self):
        rows = witnessed_rows()
        with tempfile.TemporaryDirectory() as directory:
            cache = os.path.join(directory, "runout.csv")
            sweep._write_runout(cache, rows)
            with open(cache, newline="") as source:
                cached = list(csv.DictReader(source))
            cached[0]["runout_norm"] = "99.0"
            with open(cache, "w", newline="") as target:
                writer = csv.DictWriter(target, fieldnames=cached[0].keys())
                writer.writeheader()
                writer.writerows(cached)
            with mock.patch.object(sweep, "RUNOUT_CSV", cache), \
                 mock.patch.object(sweep, "derive_dirt_ensemble", return_value=rows):
                with self.assertRaisesRegex(ValueError, "disagrees with raw witnesses"):
                    sweep.load_verified_dirt()

    def test_complete_raw_witnesses_are_used_when_the_cache_is_absent(self):
        rows = witnessed_rows()
        with tempfile.TemporaryDirectory() as directory, \
             mock.patch.object(sweep, "RUNOUT_CSV", os.path.join(directory, "runout.csv")), \
             mock.patch.object(sweep, "derive_dirt_ensemble", return_value=rows):
            self.assertEqual(sweep.load_verified_dirt(), rows)


if __name__ == "__main__":
    unittest.main()
