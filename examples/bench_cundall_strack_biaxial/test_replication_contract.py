import unittest

from replication_contract import REQUIRED_SERIES, decide
from sweep import audit_external_evidence


class ReplicationContractTests(unittest.TestCase):
    def test_primary_source_is_not_promoted_from_snapshots_to_a_response_series(self):
        missing = set(audit_external_evidence())
        decision = decide("Cundall-Strack 1979", {
            series: series not in missing for series in REQUIRED_SERIES
        })
        self.assertFalse(decision.eligible)
        self.assertEqual(missing, set(REQUIRED_SERIES))

    def test_complete_equivalent_data_is_admissible_without_a_tolerance(self):
        decision = decide("future independent data", {
            series: True for series in REQUIRED_SERIES
        })
        self.assertTrue(decision.eligible)


if __name__ == "__main__":
    unittest.main()
