#!/usr/bin/env python3
"""Regression coverage for DIRT executable provenance in raw witnesses."""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep


class DirtBinaryReceiptTests(unittest.TestCase):
    def test_identity_is_content_bound_not_just_a_target_path(self):
        """Replacing a rebuilt target at the same path changes its identity."""
        with tempfile.TemporaryDirectory() as directory:
            binary = os.path.join(directory, "bench_column_collapse")
            with open(binary, "wb") as output:
                output.write(b"first recorder")
            first = sweep.dirt_binary_identity(binary)

            with open(binary, "wb") as output:
                output.write(b"different recorder")
            second = sweep.dirt_binary_identity(binary)

        self.assertEqual(first["path"], second["path"])
        self.assertNotEqual(first["sha256"], second["sha256"])

    def test_identity_rejects_missing_recorder(self):
        with self.assertRaisesRegex(ValueError, "DIRT executable"):
            sweep.dirt_binary_identity("/nonexistent/bench_column_collapse")


if __name__ == "__main__":
    unittest.main()
