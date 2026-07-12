#!/usr/bin/env python3
"""Validate public SimulationFixture measurements against a fixed contract."""
import csv
import math
import subprocess
from pathlib import Path
import matplotlib.pyplot as plt

ROOT = Path(__file__).parent
REPO = ROOT.parent.parent
TOL = 1e-15
FIELDS = ("atom_rows", "nlocal", "natoms", "dem_rows", "materials", "pair_rows", "pair_columns", "dt")

def read_rows(text):
    return {row["case"]: row for row in csv.DictReader(text.splitlines())}

expected = read_rows((ROOT / "data" / "fixture_contract.csv").read_text())
result = subprocess.run(
    ["cargo", "run", "--quiet", "--example", "simulation_fixture_validation", "--no-default-features", "--features", "precision-double"],
    cwd=REPO, check=True, text=True, capture_output=True,
)
measured = read_rows(result.stdout)
if expected.keys() != measured.keys():
    raise SystemExit(f"case mismatch: expected {expected.keys()}, measured {measured.keys()}")

errors, labels = [], []
for case in expected:
    for field in FIELDS:
        want, got = float(expected[case][field]), float(measured[case][field])
        errors.append(abs(got - want) / max(abs(want), 1.0))
        labels.append(f"{case}\n{field}")
    for field in ("csr_offsets", "csr_indices"):
        errors.append(0.0 if measured[case][field] == expected[case][field] else math.inf)
        labels.append(f"{case}\n{field}")

failed = [label for label, error in zip(labels, errors) if error > TOL]
fig, ax = plt.subplots(figsize=(12, 5))
ax.semilogy(range(len(errors)), [max(error, 1e-18) for error in errors], "o", color="#2f6f73", label="public fixture measurement vs contract")
ax.axhline(TOL, color="crimson", linestyle="--", label=f"PASS limit ({TOL:g} relative error)")
ax.set_xticks(range(len(labels)), labels, rotation=65, ha="right", fontsize=8)
ax.set_ylabel("relative error (exact matches at 1e-18 floor)")
ax.set_title("SimulationFixture structural contract through its public API")
ax.legend()
fig.tight_layout()
(ROOT / "plots").mkdir(exist_ok=True)
fig.savefig(ROOT / "plots" / "fixture_contract.png", dpi=160)
print(f"SIMULATION FIXTURE: {'PASS' if not failed else 'FAIL'} — {len(errors) - len(failed)}/{len(errors)} structural measurements match the committed contract")
if failed:
    raise SystemExit("failed: " + ", ".join(failed))
