#!/usr/bin/env python3
"""Regression coverage for non-numerical column-collapse campaign status."""

import os
import sys
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep


class CampaignStatusTests(unittest.TestCase):
    def test_partial_raw_evidence_is_incomplete_not_a_physical_failure(self):
        admitted = (sweep.ASPECTS[0], sweep.SEEDS[0])

        def evidence(aspect, seed):
            return None if (aspect, seed) == admitted else "missing release witness"

        with mock.patch.object(sweep, "checked_protocol_manifest",
                               return_value={"protocol_sha256": "canonical"}), \
             mock.patch.object(sweep, "_case_evidence_error", side_effect=evidence):
            status = sweep.campaign_status()

        self.assertEqual(status["state"], "INCOMPLETE")
        self.assertEqual(status["protocol_sha256"], "canonical")
        self.assertEqual(status["admitted_cases"], [{"aspect": 0.5, "seed": 0}])
        self.assertEqual(len(status["inadmissible_cases"]), 32)

    def test_complete_raw_ensemble_is_only_ready_for_graph(self):
        with mock.patch.object(sweep, "checked_protocol_manifest",
                               return_value={"protocol_sha256": "canonical"}), \
             mock.patch.object(sweep, "_case_evidence_error", return_value=None):
            status = sweep.campaign_status()

        self.assertEqual(status["state"], "READY_FOR_GRAPH")
        self.assertEqual(len(status["admitted_cases"]), 33)
        self.assertEqual(status["inadmissible_cases"], [])

    def test_bad_source_provenance_is_unprepared(self):
        with mock.patch.object(sweep, "checked_protocol_manifest",
                               side_effect=ValueError("rough-base source does not match protocol")):
            status = sweep.campaign_status()

        self.assertEqual(status["state"], "UNPREPARED")
        self.assertIn("rough-base", status["reason"])


if __name__ == "__main__":
    unittest.main()
