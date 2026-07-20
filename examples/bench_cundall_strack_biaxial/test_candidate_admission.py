import pathlib
import tempfile
import unittest

from candidate_admission import (
    REQUIRED_PROTOCOL_FIELDS,
    audit_archived_lammps_deck,
    audit_archived_lammps_ledger,
    decide,
)


LEDGER = pathlib.Path(__file__).with_name("data") / "lammps_22jul2025_periodic_candidate.csv"
DECK = pathlib.Path(__file__).with_name("data") / "lammps_22jul2025_periodic_candidate.lmp"


class CandidateAdmissionTests(unittest.TestCase):
    def test_periodic_lammps_candidate_is_rejected_before_response_scoring(self):
        decision = decide(LEDGER)
        self.assertFalse(decision.eligible)
        self.assertEqual(set(decision.failures), set(REQUIRED_PROTOCOL_FIELDS))

    def test_admission_ledger_is_complete_not_a_selective_protocol_check(self):
        decision = decide(LEDGER)
        self.assertIn("lateral_boundary_and_resultant", decision.failures)
        self.assertIn("contact_or_fabric_evolution", decision.failures)

    def test_archived_deck_authenticates_the_negative_control_protocol(self):
        self.assertIsNone(audit_archived_lammps_deck(DECK))

    def test_archived_ledger_authenticates_the_complete_rejection_record(self):
        self.assertIsNone(audit_archived_lammps_ledger(LEDGER))

    def test_tampered_deck_is_rejected_before_ledger_interpretation(self):
        with tempfile.TemporaryDirectory() as directory:
            tampered = pathlib.Path(directory) / "candidate.lmp"
            tampered.write_text(DECK.read_text() + "\n# tampered\n")
            with self.assertRaisesRegex(RuntimeError, "checksum mismatch"):
                audit_archived_lammps_deck(tampered)

    def test_tampered_ledger_is_rejected_before_protocol_interpretation(self):
        with tempfile.TemporaryDirectory() as directory:
            tampered = pathlib.Path(directory) / "candidate.csv"
            tampered.write_text(LEDGER.read_text() + "\n# tampered\n")
            with self.assertRaisesRegex(RuntimeError, "ledger checksum mismatch"):
                audit_archived_lammps_ledger(tampered)


if __name__ == "__main__":
    unittest.main()
