#!/usr/bin/env python3
"""
Validation driver for `Walls::activate_by_name`.

The Rust example holds a single sphere at fixed overlap with a named plane wall
and records the normal force over three phases:

  active_before -> deactivate_by_name -> activate_by_name

The check is intentionally narrow: with unchanged geometry and material state,
the deactivated window must be exactly force-free while the reactivated window
must recover the same nonzero force as the initial active window.
"""

import csv
import os
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_wall_activate_by_name"
CONFIG = os.path.join(SCRIPT_DIR, "config.toml")
CSV = os.path.join(SCRIPT_DIR, "data", "wall_activate_by_name_results.csv")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
PLOT = os.path.join(PLOT_DIR, "wall_activate_by_name_force.png")


def build_and_run():
    subprocess.run(
        [
            "cargo",
            "run",
            "--release",
            "--example",
            EXAMPLE,
            "--no-default-features",
            "--features",
            "precision-double",
            "--",
            CONFIG,
        ],
        cwd=REPO_ROOT,
        check=True,
    )


def load_rows():
    with open(CSV, newline="") as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        row["sample"] = int(row["sample"])
        row["wall_active"] = row["wall_active"] == "true"
        row["particle_fz"] = float(row["particle_fz"])
        row["wall_force"] = float(row["wall_force"])
    return rows


def validate(rows):
    active_before = [r for r in rows if r["phase"] == "active_before"]
    deactivated = [r for r in rows if r["phase"] == "deactivated"]
    reactivated = [r for r in rows if r["phase"] == "reactivated"]
    checks = []

    def check(name, ok, detail):
        checks.append(ok)
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")

    active_force = sum(r["particle_fz"] for r in active_before) / len(active_before)
    reactivated_force = sum(r["particle_fz"] for r in reactivated) / len(reactivated)
    max_inactive = max(abs(r["particle_fz"]) for r in deactivated)
    max_wall_inactive = max(abs(r["wall_force"]) for r in deactivated)
    rel_recovery = abs(reactivated_force - active_force) / active_force

    check(
        "active wall produces a nonzero repulsive particle force",
        active_force > 0.0 and all(r["wall_active"] for r in active_before),
        f"mean active Fz = {active_force:.6e} N",
    )
    check(
        "deactivate_by_name makes the named wall force-free",
        max_inactive < 1.0e-14
        and max_wall_inactive < 1.0e-14
        and not any(r["wall_active"] for r in deactivated),
        f"max inactive particle |Fz| = {max_inactive:.3e} N, wall |F| = {max_wall_inactive:.3e} N",
    )
    check(
        "activate_by_name restores the original force response",
        rel_recovery < 1.0e-12 and all(r["wall_active"] for r in reactivated),
        f"mean reactivated Fz = {reactivated_force:.6e} N, rel recovery error = {rel_recovery:.3e}",
    )

    ok = all(checks)
    print(f"\n{sum(checks)}/{len(checks)} checks passed")
    print("ALL CHECKS PASSED" if ok else "CHECKS FAILED")
    return ok


def plot(rows):
    os.makedirs(PLOT_DIR, exist_ok=True)
    x = [r["sample"] for r in rows]
    fz = [r["particle_fz"] for r in rows]
    expected = [fz[0] if r["expected_response"] == "nonzero" else 0.0 for r in rows]

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(x, expected, color="0.55", lw=4, alpha=0.35, label="expected on/off/on response")
    ax.plot(x, fz, marker="o", color="#1f77b4", lw=1.8, label="DIRT particle Fz")
    for label in ("deactivated", "reactivated"):
        xs = [r["sample"] for r in rows if r["phase"] == label]
        if xs:
            ax.axvline(min(xs) - 0.5, color="0.2", lw=0.9, ls="--")
            ax.text(min(xs), max(fz) * 0.92, label, fontsize=9)
    ax.set_xlabel("sample")
    ax.set_ylabel("normal wall force on particle (N)")
    ax.set_title("Named wall deactivation and activate_by_name reactivation")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(PLOT, dpi=160)
    print(f"[plot] {PLOT}")


def main():
    if len(sys.argv) < 2 or sys.argv[1] != "graph":
        build_and_run()
    rows = load_rows()
    ok = validate(rows)
    plot(rows)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
