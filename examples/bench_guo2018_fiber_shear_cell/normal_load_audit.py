#!/usr/bin/env python3
"""Independently derive Guo's gravity load from reported stress and cell area.

Guo et al. (doi:10.1002/aic.16397, pp. 5--8) specifies a periodic
64-mm-by-36-mm control cell and applies normal stress by assigning a weight to
the vertically free upper wall.  For a prescribed stress, the *assigned load*
is consequently determined by sigma*A.  This is not evidence for the missing
wall-sphere diameter or layout, and it must not be mistaken for a runnable DIRT
wall assembly.
"""
import argparse
import math
import unittest

G = 9.81
LX = 0.064
LZ = 0.036
PRINTED_CASE_PA = 1735.0


def mass_for_normal_stress(stress_pa: float, lx: float = LX, lz: float = LZ) -> float:
    if stress_pa <= 0 or lx <= 0 or lz <= 0:
        raise ValueError("stress and periodic cell dimensions must be positive")
    return stress_pa * lx * lz / G


def normal_stress_for_mass(mass_kg: float, lx: float = LX, lz: float = LZ) -> float:
    if mass_kg <= 0 or lx <= 0 or lz <= 0:
        raise ValueError("mass and periodic cell dimensions must be positive")
    return mass_kg * G / (lx * lz)


class NormalLoadTests(unittest.TestCase):
    def test_printed_case_round_trips_without_a_wall_geometry_assumption(self):
        mass = mass_for_normal_stress(PRINTED_CASE_PA)
        self.assertAlmostEqual(mass, 0.4074862385321101, places=15)
        self.assertAlmostEqual(normal_stress_for_mass(mass), PRINTED_CASE_PA, places=10)

    def test_invalid_inputs_are_not_silently_accepted(self):
        with self.assertRaises(ValueError):
            mass_for_normal_stress(0)
        with self.assertRaises(ValueError):
            normal_stress_for_mass(-1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--stress-pa", type=float, default=PRINTED_CASE_PA)
    args = parser.parse_args()
    if args.self_test:
        raise SystemExit(not unittest.main(argv=["normal_load_audit"], exit=False).result.wasSuccessful())
    mass = mass_for_normal_stress(args.stress_pa)
    recovered = normal_stress_for_mass(mass)
    if not math.isclose(recovered, args.stress_pa, rel_tol=0.0, abs_tol=1e-10):
        raise RuntimeError("normal-load round trip failed")
    print(f"NORMAL LOAD: {args.stress_pa:g} Pa over {LX * 1000:g} mm x {LZ * 1000:g} mm -> {mass:.15g} kg at g={G:g} m/s^2")
    print("LIMIT: this derives only the assigned gravity load; wall-sphere diameter/layout remain unreported.")


if __name__ == "__main__":
    main()
