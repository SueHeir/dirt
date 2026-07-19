import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from protocol_admission import admission_failures


class ProtocolAdmissionTests(unittest.TestCase):
    def test_current_surrogate_is_not_eligible_for_a_replication_score(self):
        failures = admission_failures()
        self.assertIn("dimension", failures)
        self.assertIn("lateral_boundary", failures)
        self.assertIn("stress_observable", failures)
        self.assertIn("fabric_observable", failures)


if __name__ == "__main__":
    unittest.main()
