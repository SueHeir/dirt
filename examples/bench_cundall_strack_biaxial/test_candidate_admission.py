import pathlib
import unittest

from candidate_admission import REQUIRED_PROTOCOL_FIELDS, decide


LEDGER = pathlib.Path(__file__).with_name("data") / "lammps_22jul2025_periodic_candidate.csv"


class CandidateAdmissionTests(unittest.TestCase):
    def test_periodic_lammps_candidate_is_rejected_before_response_scoring(self):
        decision = decide(LEDGER)
        self.assertFalse(decision.eligible)
        self.assertEqual(set(decision.failures), set(REQUIRED_PROTOCOL_FIELDS))

    def test_admission_ledger_is_complete_not_a_selective_protocol_check(self):
        decision = decide(LEDGER)
        self.assertIn("lateral_boundary_and_resultant", decision.failures)
        self.assertIn("contact_or_fabric_evolution", decision.failures)


if __name__ == "__main__":
    unittest.main()
