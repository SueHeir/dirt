#!/usr/bin/env python3
"""Regression checks for the finite gate above the rough support layer."""

import importlib.util
import os
import tempfile
import tomllib
import unittest


HERE = os.path.dirname(os.path.abspath(__file__))
SPEC = importlib.util.spec_from_file_location("column_collapse_sweep", os.path.join(HERE, "sweep.py"))
sweep = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(sweep)


class GateSupportSeparationTests(unittest.TestCase):
    def test_gate_starts_above_frozen_base_and_below_mobile_source(self):
        base = sweep.rough_base_positions()
        active = sweep.active_column_positions(sweep.n_particles(0.5), 0.5, 0)
        self.assertTrue(all(z < sweep.GATE_Z_LOW for _, _, z in base))
        self.assertTrue(all(z > sweep.GATE_Z_LOW for _, _, z in active))

    def test_rendered_gate_has_same_lower_bound(self):
        # Keep this small rendering check independent of generated campaign files.
        with tempfile.TemporaryDirectory() as tmp:
            active = os.path.join(tmp, "active.csv")
            base = os.path.join(tmp, "base.csv")
            open(active, "w").close()
            open(base, "w").close()
            rendered = sweep.TOML_TEMPLATE.format(
                aspect=0.5, count=1, seed=0, youngs=sweep.YOUNGS_MOD,
                poisson=sweep.POISSON, restitution=sweep.RESTITUTION,
                friction=sweep.FRICTION, radius=sweep.RADIUS, density=sweep.DENSITY,
                l0=f"{sweep.L0:.4f}", w=f"{sweep.W:.4f}", y_high=f"{sweep.W + 0.003:.4f}",
                rough_base=base, active_z_low=f"{sweep.BASE_Z + sweep.RADIUS:.4f}",
                base_select_z=f"{sweep.BASE_SELECT_Z:.4f}",
                gate_z_low=f"{sweep.GATE_Z_LOW:.4f}", insert_top="0.0100", z_high="0.2000",
                gate_x_high=f"{sweep.L0 + 5.0 * sweep.RADIUS:.4f}",
                gate_y_high=f"{sweep.W + 0.0060:.4f}",
                active_column=active, output_dir=tmp,
                settle_dt=f"{sweep.SETTLE_DT:.3e}", collapse_dt=f"{sweep.COLLAPSE_DT:.3e}",
                preparation_max_displacement=f"{sweep.PREPARATION_MAX_DISPLACEMENT:.3e}",
                settle_steps=sweep.SETTLE_STEPS, collapse_steps=sweep.COLLAPSE_STEPS,
            )
        gate = tomllib.loads(rendered)["wall"][-1]
        self.assertEqual(gate["type"], "region")
        self.assertFalse(gate["inside"])
        self.assertEqual(gate["region"]["type"], "block")
        self.assertEqual(gate["region"]["min"][0], sweep.L0)
        self.assertEqual(gate["region"]["min"][2], sweep.GATE_Z_LOW)
        self.assertGreater(gate["region"]["max"][0], sweep.L0)
        self.assertGreater(gate["region"]["max"][1], sweep.W)

    def test_rendered_type_zero_is_the_frozen_rough_base(self):
        """The dynamic type group must immobilize support, never mobile grains."""
        rendered = sweep.TOML_TEMPLATE.format(
            aspect=0.5, count=1, seed=0, youngs=sweep.YOUNGS_MOD,
            poisson=sweep.POISSON, restitution=sweep.RESTITUTION,
            friction=sweep.FRICTION, radius=sweep.RADIUS, density=sweep.DENSITY,
            l0=f"{sweep.L0:.4f}", w=f"{sweep.W:.4f}", y_high=f"{sweep.W + 0.003:.4f}",
            rough_base="rough.csv", active_z_low=f"{sweep.BASE_Z + sweep.RADIUS:.4f}",
            base_select_z=f"{sweep.BASE_SELECT_Z:.4f}",
            gate_z_low=f"{sweep.GATE_Z_LOW:.4f}", insert_top="0.0100", z_high="0.2000",
            gate_x_high=f"{sweep.L0 + 5.0 * sweep.RADIUS:.4f}",
            gate_y_high=f"{sweep.W + 0.0060:.4f}",
            active_column="active.csv", output_dir="output",
            settle_dt=f"{sweep.SETTLE_DT:.3e}", collapse_dt=f"{sweep.COLLAPSE_DT:.3e}",
            preparation_max_displacement=f"{sweep.PREPARATION_MAX_DISPLACEMENT:.3e}",
            settle_steps=sweep.SETTLE_STEPS, collapse_steps=sweep.COLLAPSE_STEPS,
        )
        config = tomllib.loads(rendered)
        self.assertEqual(config["dem"]["materials"][0]["name"], "rough_glass")
        self.assertEqual(config["group"][0]["name"], "rough_base")
        self.assertEqual(config["group"][0]["type"], [0])
        self.assertEqual(config["particles"]["insert"][0]["material"], "rough_glass")
        self.assertEqual(config["particles"]["insert"][1]["material"], "glass")


if __name__ == "__main__":
    unittest.main()
