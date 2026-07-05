#!/usr/bin/env python3
"""
Deterministic clump Monte Carlo inertia sampler validation.

Commands:
    python3 examples/bench_clump_inertia_sampler/sweep.py
    python3 examples/bench_clump_inertia_sampler/sweep.py start
    python3 examples/bench_clump_inertia_sampler/sweep.py graph
"""

import csv
import os
import subprocess
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
CONFIG = os.path.join(SCRIPT_DIR, "config.toml")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
CSV = os.path.join(DATA_DIR, "inertia_sampler.csv")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
PLOT = os.path.join(PLOT_DIR, "inertia_sampler_determinism.png")


def start():
    os.makedirs(DATA_DIR, exist_ok=True)
    subprocess.run(
        [
            "cargo",
            "run",
            "--release",
            "--example",
            "bench_clump_inertia_sampler",
            "--no-default-features",
            "--features",
            "precision-double",
            "--",
            CONFIG,
        ],
        cwd=REPO_ROOT,
        check=True,
    )


def read_rows():
    rows = []
    with open(CSV, newline="") as f:
        for row in csv.DictReader(f):
            row["sample_count"] = int(row["sample_count"])
            row["mass_rel_err"] = float(row["mass_rel_err"])
            row["max_diag_rel_err"] = float(row["max_diag_rel_err"])
            row["bitwise_repeat"] = row["bitwise_repeat"] == "true"
            rows.append(row)
    return rows


def graph():
    if not os.path.exists(CSV):
        start()

    rows = read_rows()
    repeat_rows = [r for r in rows if r["mode"] in ("default_repeat", "explicit_seed_repeat")]
    spread_rows = [r for r in rows if r["mode"] == "seed_spread"]

    repeat_ok = all(r["bitwise_repeat"] for r in repeat_rows)
    counts = sorted({r["sample_count"] for r in spread_rows})
    tolerance = 0.05
    final_count = max(counts)
    final_errors = [r["max_diag_rel_err"] for r in spread_rows if r["sample_count"] == final_count]
    spread_ok = max(final_errors) <= tolerance

    os.makedirs(PLOT_DIR, exist_ok=True)
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"matplotlib unavailable, skipped plot: {e}")
        return 1

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)

    modes = ["default_repeat", "explicit_seed_repeat"]
    labels = ["default seed", "explicit seed"]
    failures = [
        sum(1 for r in repeat_rows if r["mode"] == mode and not r["bitwise_repeat"])
        for mode in modes
    ]
    ax0.bar(labels, failures, color=["#4c78a8", "#f58518"])
    ax0.axhline(0.0, color="black", linewidth=0.8)
    ax0.set_ylabel("bitwise repeat failures")
    ax0.set_title("Repeatability gate")
    ax0.set_ylim(0, max(1, max(failures) + 1))
    for i, v in enumerate(failures):
        ax0.text(i, v + 0.03, str(v), ha="center", va="bottom")

    for count in counts:
        ys = [r["max_diag_rel_err"] for r in spread_rows if r["sample_count"] == count]
        xs = [count] * len(ys)
        ax1.scatter(xs, ys, s=20, alpha=0.65, edgecolors="none")
        mean = sum(ys) / len(ys)
        ax1.plot([count * 0.88, count * 1.12], [mean, mean], color="black", linewidth=1.0)
    ax1.axhline(tolerance, color="#d62728", linestyle="--", linewidth=1.2, label="5% tolerance")
    ax1.set_xscale("log")
    ax1.set_xlabel("Monte Carlo samples")
    ax1.set_ylabel("max diagonal inertia rel. error")
    ax1.set_title("Seed-to-seed spread vs analytic sphere")
    ax1.legend(loc="upper right")

    status = "PASS" if repeat_ok and spread_ok else "FAIL"
    fig.suptitle(f"Clump inertia sampler determinism: {status}")
    fig.savefig(PLOT, dpi=180)
    print(f"wrote {PLOT}")
    print(f"repeatability: {'PASS' if repeat_ok else 'FAIL'}")
    print(
        f"seed spread at n={final_count}: max_rel_err={max(final_errors):.4%} "
        f"<= {tolerance:.1%} -> {'PASS' if spread_ok else 'FAIL'}"
    )
    print("ALL CHECKS PASSED" if repeat_ok and spread_ok else "CHECKS FAILED")
    return 0 if repeat_ok and spread_ok else 1


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    if cmd == "start":
        start()
        return 0
    if cmd == "graph":
        return graph()
    if cmd == "all":
        start()
        return graph()
    print("usage: sweep.py [all|start|graph]", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
