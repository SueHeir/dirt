#!/usr/bin/env python3
"""Fail-closed validator for a Guo rubber-cord shear-cell campaign.

This program deliberately consumes solver-written histories only.  It has no
DIRT-result fixture and it never substitutes a requested servo load for the
measured lid reaction.
"""
import argparse
import csv
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
REFERENCE = HERE / "data" / "guo_2019_rubber_cord.csv"
LOADS = (651.0, 1735.0, 3470.0)
SHEAR_STRAIN = 0.50
WINDOW = 40
LOAD_TOL = 0.15
SHEAR_TOL = 0.30
PHI_TOL = 0.05


def reference():
    values = {}
    with REFERENCE.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if row["source"] != "experiment" or row["material"] != "rubber_cord":
                continue
            values.setdefault(float(row["normal_stress_pa"]), {})[row["observable"]] = float(row["value"])
    if set(values) != set(LOADS) or any(set(v) != {"shear_stress_pa", "solid_fraction"} for v in values.values()):
        raise ValueError("incomplete external reference series")
    return values


def mean(rows, field):
    return sum(float(r[field]) for r in rows) / len(rows)


def measure(path, pressure, width):
    """Return post-drive solver observables after checking the full protocol."""
    with Path(path).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    required = {"stage", "shear_strain", "normal_stress_pa", "shear_stress_pa", "solid_fraction", "n_atoms"}
    if not rows or not required <= set(rows[0]):
        raise ValueError(f"{path}: not a complete solver cell history")
    expected_atoms = 8500 if width == 64 else 12750 if width == 96 else None
    if expected_atoms is None or any(int(r["n_atoms"]) != expected_atoms for r in rows):
        raise ValueError(f"{path}: unexpected or changing fibre population")
    stages = {r["stage"] for r in rows}
    if not {"settle", "normal_load", "shear"} <= stages:
        raise ValueError(f"{path}: missing settle, measured-load, or driven-shear stage")
    loading = [r for r in rows if r["stage"] == "normal_load"]
    if len(loading) < WINDOW or abs(mean(loading[-WINDOW:], "normal_stress_pa") - pressure) / pressure > LOAD_TOL:
        raise ValueError(f"{path}: measured lid load does not qualify at {pressure:g} Pa")
    steady = [r for r in rows if r["stage"] == "shear" and float(r["shear_strain"]) >= SHEAR_STRAIN]
    if len(steady) < WINDOW:
        raise ValueError(f"{path}: no declared steady-strain window")
    final = steady[-WINDOW:]
    if abs(mean(final, "normal_stress_pa") - pressure) / pressure > LOAD_TOL:
        raise ValueError(f"{path}: lid load drifted outside tolerance during shear")
    return {name: mean(final, name) for name in ("normal_stress_pa", "shear_stress_pa", "solid_fraction")}


def compare(cases):
    """Compare complete campaign observations to external figures and size check."""
    if set(cases) != {(p, w) for p in LOADS for w in (64, 96)}:
        raise ValueError("campaign must contain all three loads at 64 and 96 mm")
    refs = reference()
    report = []
    for pressure in LOADS:
        observed, ref = cases[pressure, 64], refs[pressure]
        shear_error = abs(observed["shear_stress_pa"] - ref["shear_stress_pa"]) / ref["shear_stress_pa"]
        phi_error = abs(observed["solid_fraction"] - ref["solid_fraction"])
        if shear_error > SHEAR_TOL or phi_error > PHI_TOL:
            raise ValueError(f"{pressure:g} Pa: external Fig. 6/7 gate failed (tau={shear_error:.3f}, phi={phi_error:.3f})")
        larger = cases[pressure, 96]
        if abs(observed["shear_stress_pa"] - larger["shear_stress_pa"]) / max(abs(observed["shear_stress_pa"]), 1.0) > SHEAR_TOL or abs(observed["solid_fraction"] - larger["solid_fraction"]) > PHI_TOL:
            raise ValueError(f"{pressure:g} Pa: 64/96-mm sensitivity gate failed")
        report.append((pressure, observed, ref, shear_error, phi_error))
    return report


def parse_case(value):
    pressure, width, path = value.split(":", 2)
    return float(pressure), int(width), Path(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", default=[], metavar="PA:MM:HISTORY", help="one solver history; provide all six cases")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(ValidatorTests)
        raise SystemExit(not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful())
    cases = {}
    for encoded in args.case:
        pressure, width, path = parse_case(encoded)
        cases[pressure, width] = measure(path, pressure, width)
    for pressure, observed, ref, shear_error, phi_error in compare(cases):
        print(f"{pressure:.0f} Pa: tau={observed['shear_stress_pa']:.1f}/{ref['shear_stress_pa']:.1f} (rel={shear_error:.3f}); phi={observed['solid_fraction']:.3f}/{ref['solid_fraction']:.3f} (abs={phi_error:.3f})")
    print("PASS: solver histories satisfy measured-load, external-reference, and size-sensitivity gates")


class ValidatorTests(unittest.TestCase):
    def history(self, normal=651.0, stages=True):
        rows = []
        if stages:
            rows.append(["settle", 0, normal, 0, .30, 8500])
            rows.extend([["normal_load", 0, normal, 0, .33, 8500]] * WINDOW)
        rows.extend([["shear", .51, normal, 480, .33, 8500]] * WINDOW)
        fd = tempfile.NamedTemporaryFile("w", newline="", delete=False)
        with fd:
            writer = csv.writer(fd); writer.writerow(["stage", "shear_strain", "normal_stress_pa", "shear_stress_pa", "solid_fraction", "n_atoms"]); writer.writerows(rows)
        return Path(fd.name)

    def test_measures_real_history(self):
        path = self.history()
        self.assertAlmostEqual(measure(path, 651, 64)["shear_stress_pa"], 480)
        path.unlink()

    def test_rejects_missing_protocol_stage(self):
        path = self.history(stages=False)
        with self.assertRaisesRegex(ValueError, "missing"):
            measure(path, 651, 64)
        path.unlink()


if __name__ == "__main__":
    main()
