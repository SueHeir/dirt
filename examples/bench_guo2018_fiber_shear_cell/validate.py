#!/usr/bin/env python3
"""Fail-closed validator for a Guo rubber-cord shear-cell campaign.

This program deliberately consumes solver-written histories only.  It has no
DIRT-result fixture and it never substitutes a requested servo load for the
measured lid reaction.
"""
import argparse
import csv
import hashlib
import json
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
# Guo et al. Fig. 6 prints the rubber-cord experimental fit
# tau = 0.83 sigma_yy - 60 Pa.  This independently constrains the stored
# figure digitization; it is not derived from any DIRT observation.
FIG6_RUBBER_CORD_SLOPE = 0.83
FIG6_RUBBER_CORD_INTERCEPT_PA = -60.0
FIG6_DIGITIZATION_REL_TOL = 0.03


def reference():
    values = {}
    with REFERENCE.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if row["source"] != "experiment" or row["material"] != "rubber_cord":
                continue
            values.setdefault(float(row["normal_stress_pa"]), {})[row["observable"]] = float(row["value"])
    if set(values) != set(LOADS) or any(set(v) != {"shear_stress_pa", "solid_fraction"} for v in values.values()):
        raise ValueError("incomplete external reference series")
    for pressure, observables in values.items():
        figure_fit = FIG6_RUBBER_CORD_SLOPE * pressure + FIG6_RUBBER_CORD_INTERCEPT_PA
        if abs(observables["shear_stress_pa"] - figure_fit) / figure_fit > FIG6_DIGITIZATION_REL_TOL:
            raise ValueError(f"Fig. 6 digitization at {pressure:g} Pa disagrees with its printed fit")
    return values


def mean(rows, field):
    return sum(float(r[field]) for r in rows) / len(rows)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def measure(path, pressure, width, expected_atoms=None):
    """Return post-drive solver observables after checking the full protocol."""
    with Path(path).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    required = {"stage", "shear_strain", "normal_stress_pa", "shear_stress_pa", "solid_fraction", "n_atoms"}
    if not rows or not required <= set(rows[0]):
        raise ValueError(f"{path}: not a complete solver cell history")
    if expected_atoms is None:
        # 64-mm population is the Table-2 500-fibre case.  The 96-mm
        # sensitivity input preserves areal loading (750 fibres); it is not
        # presented as a primary-source population.
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


def measure_case(case_dir, pressure, width):
    """Measure a runner-receipted case, never a free-standing CSV."""
    case_dir = Path(case_dir)
    manifest_path, receipt_path = case_dir / "case_manifest.json", case_dir / "solver_receipt.json"
    if not manifest_path.is_file() or not receipt_path.is_file():
        raise ValueError(f"{case_dir}: missing run_case manifest or successful solver receipt")
    manifest, receipt = json.loads(manifest_path.read_text()), json.loads(receipt_path.read_text())
    if manifest.get("pressure_pa") != pressure or manifest.get("width_mm") != width:
        raise ValueError(f"{case_dir}: manifest does not identify the requested case")
    history = case_dir / manifest.get("solver_history", "")
    if not history.is_file() or receipt.get("history_sha256") != sha256(history):
        raise ValueError(f"{case_dir}: history differs from the successful solver receipt")
    for name, digest in manifest.get("input_sha256", {}).items():
        if not (case_dir / name).is_file() or sha256(case_dir / name) != digest:
            raise ValueError(f"{case_dir}: solver input {name} differs from its manifest")
    if receipt.get("input_sha256") != manifest.get("input_sha256"):
        raise ValueError(f"{case_dir}: receipt is not bound to this solver input")
    return measure(history, pressure, width, manifest.get("expected_global_atoms"))


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
    # A manifest/receipt can establish provenance only after the input matches
    # the paper's stated periodic planar DEM reduction.  This is distinct from
    # the experimental ring apparatus and prevents the current cylindrical
    # cup surrogate from being compared to Figs. 6/7.
    if args.case:
        from source_contract import require_published_control_cell
        require_published_control_cell(HERE / "config.toml")
    cases = {}
    for encoded in args.case:
        pressure, width, path = parse_case(encoded)
        cases[pressure, width] = measure_case(path, pressure, width)
    for pressure, observed, ref, shear_error, phi_error in compare(cases):
        print(f"{pressure:.0f} Pa: tau={observed['shear_stress_pa']:.1f}/{ref['shear_stress_pa']:.1f} (rel={shear_error:.3f}); phi={observed['solid_fraction']:.3f}/{ref['solid_fraction']:.3f} (abs={phi_error:.3f})")
    print("PASS: solver histories satisfy measured-load, external-reference, and size-sensitivity gates")


class ValidatorTests(unittest.TestCase):
    def test_external_digitization_agrees_with_printed_figure_fit(self):
        self.assertEqual(set(reference()), set(LOADS))

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

    def test_rejects_missing_protocol_stage(self):
        path = self.history(stages=False)
        with self.assertRaisesRegex(ValueError, "missing"):
            measure(path, 651, 64)
        path.unlink()

    def test_rejects_requested_load_when_lid_did_not_measure_it(self):
        path = self.history(normal=0.0)
        with self.assertRaisesRegex(ValueError, "measured lid load"):
            measure(path, 651, 64)
        path.unlink()


if __name__ == "__main__":
    main()
