#!/usr/bin/env python3
import csv
import os
import subprocess
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
CSV_PATH = os.path.join(SCRIPT_DIR, "data", "wall_twisting_parity.csv")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")


def start():
    cmd = [
        "cargo",
        "run",
        "--release",
        "--example",
        "bench_wall_twisting_parity",
        "--no-default-features",
        "--features",
        "precision-double",
        "--",
        os.path.join(SCRIPT_DIR, "config.toml"),
    ]
    subprocess.run(cmd, cwd=REPO_ROOT, check=True)


def graph():
    rows = []
    with open(CSV_PATH) as f:
        for row in csv.DictReader(f):
            row["torque_x"] = float(row["torque_x"])
            row["expected_tau_x"] = float(row["expected_tau_x"])
            row["rel_err"] = float(row["rel_err"])
            rows.append(row)

    max_err = max(r["rel_err"] for r in rows)
    passed = max_err < 1.0e-12

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(PLOT_DIR, exist_ok=True)
    labels = [r["geometry"] for r in rows]
    measured = [r["torque_x"] for r in rows]
    expected = [r["expected_tau_x"] for r in rows]
    rel_err = [r["rel_err"] for r in rows]
    gate = 1.0e-12

    fig, (ax, err_ax) = plt.subplots(
        2,
        1,
        figsize=(7, 5.5),
        sharex=True,
        gridspec_kw={"height_ratios": [2.0, 1.15]},
    )
    xs = range(len(rows))
    ax.bar([x - 0.18 for x in xs], measured, width=0.36, label="DIRT torque")
    ax.bar([x + 0.18 for x in xs], expected, width=0.36, label="plane-law reference")
    ax.set_ylabel("twisting torque x (N m)")
    ax.set_title("Wall twisting torque parity")
    ax.legend()
    ax.grid(axis="y", alpha=0.25)

    err_ax.bar(xs, rel_err, width=0.5, color="#4c78a8", label="relative error")
    err_ax.axhline(gate, color="#c43c39", linestyle="--", linewidth=1.5, label="PASS gate: 1e-12")
    err_ax.set_xticks(list(xs), labels)
    err_ax.set_ylabel("relative error")
    err_ax.set_ylim(0.0, gate * 1.25)
    err_ax.grid(axis="y", alpha=0.25)
    err_ax.legend(loc="upper right")
    err_ax.text(
        0.01,
        0.88,
        f"{'PASS' if passed else 'FAIL'}: max rel. error = {max_err:.2e}",
        transform=err_ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "wall_twisting_parity.png"), dpi=160)

    print("=== Wall twisting parity ===")
    for r in rows:
        print(
            f"  {r['geometry']:8s} torque_x={r['torque_x']:.12e} "
            f"ref={r['expected_tau_x']:.12e} rel_err={r['rel_err']:.2e}"
        )
    print(f"RESULT: {'PASS' if passed else 'FAIL'} (max_rel_err={max_err:.2e})")
    if not passed:
        sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        start()
        graph()
    elif sys.argv[1] == "start":
        start()
    elif sys.argv[1] == "graph":
        graph()
    else:
        raise SystemExit(f"unknown command: {sys.argv[1]}")
