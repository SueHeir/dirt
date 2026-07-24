#!/usr/bin/env python3
"""Regression guard for the safe, non-forced Verlet rebuild policy."""

import importlib.util
import os
import tempfile
import tomllib
import unittest


HERE = os.path.dirname(os.path.abspath(__file__))
SPEC = importlib.util.spec_from_file_location("column_collapse_sweep", os.path.join(HERE, "sweep.py"))
sweep = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(sweep)


class NeighborScheduleTests(unittest.TestCase):
    def test_generated_case_uses_displacement_driven_rebuild(self):
        with tempfile.TemporaryDirectory() as tmp:
            rendered = sweep.TOML_TEMPLATE.format(
                aspect=0.5, count=1, seed=0, youngs=sweep.YOUNGS_MOD,
                poisson=sweep.POISSON, restitution=sweep.RESTITUTION,
                friction=sweep.FRICTION, radius=sweep.RADIUS, density=sweep.DENSITY,
                l0=f"{sweep.L0:.4f}", w=f"{sweep.W:.4f}", y_high=f"{sweep.W + 0.003:.4f}",
                rough_base="rough.csv", active_z_low=f"{sweep.BASE_Z + sweep.RADIUS:.4f}",
            base_select_z=f"{sweep.BASE_SELECT_Z:.4f}", gate_z_low=f"{sweep.GATE_Z_LOW:.4f}",
            gate_x_high=f"{sweep.L0 + 5.0 * sweep.RADIUS:.4f}",
            gate_y_high=f"{sweep.W + 0.0060:.4f}",
                insert_top="0.0100", z_high="0.2000", active_column="active.csv", output_dir=tmp,
                settle_dt=f"{sweep.SETTLE_DT:.3e}", collapse_dt=f"{sweep.COLLAPSE_DT:.3e}",
                preparation_max_displacement=f"{sweep.PREPARATION_MAX_DISPLACEMENT:.3e}",
                settle_steps=sweep.SETTLE_STEPS, collapse_steps=sweep.COLLAPSE_STEPS,
            )
        neighbor = tomllib.loads(rendered)["neighbor"]
        self.assertEqual(neighbor["skin_fraction"], 1.1)
        self.assertEqual(neighbor["every"], 0)
        self.assertNotIn("check", neighbor)


if __name__ == "__main__":
    unittest.main()
