#!/usr/bin/env python3
"""Independent checks for the analytical static-sliding preflight."""

import importlib.util
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).parent
SPEC = importlib.util.spec_from_file_location("physical_feasibility", HERE / "physical_feasibility.py")
physical_feasibility = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(physical_feasibility)


class StaticSlidingFeasibilityTests(unittest.TestCase):
    def test_coulomb_limit_matches_independent_hand_values(self):
        self.assertAlmostEqual(physical_feasibility.single_contact_sliding_coefficient(22.0),
                               0.4040262258, places=9)
        self.assertAlmostEqual(physical_feasibility.single_contact_sliding_coefficient(26.0),
                               0.4877325886, places=9)

    def test_declared_glass_sliding_value_cannot_support_low_band_edge(self):
        mu_p = physical_feasibility.declared_sliding_friction(HERE / "config.toml")
        self.assertEqual(mu_p, 0.16)
        self.assertFalse(physical_feasibility.single_contact_supports_band(mu_p, (22.0, 26.0)))

    def test_positive_control_is_not_wired_to_the_failing_material(self):
        self.assertTrue(physical_feasibility.single_contact_supports_band(0.5, (22.0, 26.0)))

    def test_config_reader_rejects_non_sds_model(self):
        source = (HERE / "config.toml").read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            config = Path(directory) / "bad.toml"
            config.write_text(source.replace('rolling_model = "sds"',
                                             'rolling_model = "constant"'), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "SDS"):
                physical_feasibility.declared_sliding_friction(config)


if __name__ == "__main__":
    unittest.main()
