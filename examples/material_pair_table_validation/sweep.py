#!/usr/bin/env python3
"""Plot typed Material pair-table values against a pre-redesign golden table."""
import csv
import subprocess
from pathlib import Path
import matplotlib.pyplot as plt

ROOT = Path(__file__).parent
TOL = 1e-12

def read_file(name):
    with (ROOT / "data" / name).open() as f:
        return {r["property"]: float(r["value"]) for r in csv.DictReader(f)}

def measure_typed():
    result = subprocess.run(
        ["cargo", "test", "--quiet", "-p", "dirt_atom", "--no-default-features",
         "--features", "precision-double", "tests::typed_materials_match_independent_legacy_mixing_rules",
         "--", "--exact", "--nocapture"],
        cwd=ROOT.parent.parent, check=True, text=True, capture_output=True,
    )
    rows = [line.removeprefix("PAIR_TABLE,") for line in result.stdout.splitlines()
            if line.startswith("PAIR_TABLE,")]
    return {r["property"]: float(r["value"]) for r in csv.DictReader(rows)}

legacy, typed = read_file("legacy_pair_table.csv"), measure_typed()
assert legacy.keys() == typed.keys()
errors = {k: abs(typed[k] - legacy[k]) / max(abs(legacy[k]), 1.0) for k in legacy}
failed = [k for k, e in errors.items() if e > TOL]
print(f"PAIR TABLE: {'PASS' if not failed else 'FAIL'} — {len(legacy) - len(failed)}/{len(legacy)} properties within {TOL:g} relative tolerance")
labels = list(errors)
fig, ax = plt.subplots(figsize=(12, 4.8))
ax.semilogy(range(len(labels)), [max(errors[k], 1e-18) for k in labels], "o", label="typed vs origin/main golden")
ax.axhline(TOL, color="crimson", linestyle="--", label=f"pass limit ({TOL:g})")
ax.set_xticks(range(len(labels)), labels, rotation=65, ha="right", fontsize=8)
ax.set_ylabel("relative error (floor 1e-18)")
ax.set_title("Typed Material pair table against pre-redesign golden values")
ax.legend()
fig.tight_layout()
(ROOT / "plots").mkdir(exist_ok=True)
fig.savefig(ROOT / "plots" / "typed_vs_legacy_pair_table.png", dpi=160)
if failed:
    raise SystemExit("failed: " + ", ".join(failed))
