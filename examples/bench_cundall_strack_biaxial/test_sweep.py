"""Regression tests for the wall-reaction measurement contract."""
import importlib.util
import pathlib
import unittest


SPEC = importlib.util.spec_from_file_location(
    "cundall_sweep", pathlib.Path(__file__).with_name("sweep.py")
)
SWEEP = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(SWEEP)


def row(strain, horizontal, vertical, ratio):
    return {"axial_strain": strain, "f_h_mean": horizontal, "f_v_mean": vertical,
            "wall_force_ratio": ratio, "contacts": 1.0, "lateral_strain": -strain,
            "phi": 0.55, "coordination": 4.5}


class MeasurementContractTests(unittest.TestCase):
    def test_primary_source_rows_are_audited(self):
        self.assertEqual(SWEEP.read_reference(), {"A": 0.39, "B": 0.33})

    def test_primary_source_lacks_all_goal_required_trajectories(self):
        self.assertEqual(SWEEP.audit_external_evidence(), [
            "state_registration", "stress_or_deviatoric_path",
            "volumetric_strain_or_dilatancy_path", "contact_or_fabric_evolution",
        ])

    def test_recomputed_wall_ratio_is_accepted(self):
        passed, _ = SWEEP.evaluate([row(0.0, 2.0, 4.0, 0.5), row(0.1, 3.0, 6.0, 0.5)])
        self.assertTrue(passed)

    def test_inconsistent_precomputed_ratio_is_rejected(self):
        passed, checks = SWEEP.evaluate([row(0.0, 2.0, 4.0, 0.5), row(0.1, 3.0, 6.0, 0.6)])
        self.assertFalse(passed)
        self.assertFalse(checks["ratio_recomputed"])

    def test_lateral_wall_motion_is_required_for_a_moving_wall_cell(self):
        rows = [row(0.0, 2.0, 4.0, 0.5), row(0.1, 3.0, 6.0, 0.5)]
        rows[-1]["lateral_strain"] = 0.0
        passed, checks = SWEEP.evaluate(rows)
        self.assertFalse(passed)
        self.assertFalse(checks["lateral_compression"])

    def test_loose_insertion_transient_is_rejected(self):
        rows = [row(0.0, 2.0, 4.0, 0.5), row(0.1, 3.0, 6.0, 0.5)]
        rows[0]["phi"] = 0.49
        rows[1]["coordination"] = 3.9
        passed, checks = SWEEP.evaluate(rows)
        self.assertFalse(passed)
        self.assertFalse(checks["dense_solid_fraction"])
        self.assertFalse(checks["dense_coordination"])

    def test_only_the_primary_source_is_considered_for_replication_admission(self):
        decision = SWEEP.replication_evidence_decision()
        self.assertFalse(decision.eligible)
        self.assertEqual(decision.candidate, "Cundall--Strack 1979 primary source")


if __name__ == "__main__":
    unittest.main()
