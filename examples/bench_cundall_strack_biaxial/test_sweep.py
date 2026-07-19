"""Regression tests for the wall-reaction measurement contract."""
import importlib.util
import pathlib
import unittest
from unittest import mock


SPEC = importlib.util.spec_from_file_location(
    "cundall_sweep", pathlib.Path(__file__).with_name("sweep.py")
)
SWEEP = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(SWEEP)


def row(strain, horizontal, vertical, ratio):
    return {"axial_strain": strain, "f_h_mean": horizontal, "f_v_mean": vertical,
            "wall_force_ratio": ratio, "contacts": 1.0}


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

    def test_independent_comparison_reports_raw_disagreement(self):
        rows = [row(0.0, 1.0, 1.0, 1.0), row(0.07, 1.0, 2.0, 0.5)]
        reference = [
            {"axial_strain": 0.0, "syy": 1.0},
            {"axial_strain": 0.07, "syy": 3.0},
        ]
        with mock.patch.object(SWEEP, "read", return_value=reference):
            comparison = SWEEP.compare_independent_lammps(rows)
        self.assertEqual(comparison["samples"], 12)
        self.assertGreater(comparison["normalized_axial_rmse"], 0.0)


if __name__ == "__main__":
    unittest.main()
