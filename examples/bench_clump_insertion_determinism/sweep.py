#!/usr/bin/env python3
"""Config-level clump insertion determinism gate."""

import csv
import math
import os
import shutil
import subprocess
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_clump_insertion_determinism"
CONFIG = os.path.join(SCRIPT_DIR, "config.toml")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
PLOT = os.path.join(PLOT_DIR, "clump_insertion_determinism.png")


def config_for(seed):
    with open(CONFIG) as f:
        text = f.read()
    return text.replace("seed = 20260705", f"seed = {seed}")


def run_case(label, seed):
    os.makedirs(DATA_DIR, exist_ok=True)
    cfg = os.path.join(DATA_DIR, f"{label}.toml")
    out = os.path.join(DATA_DIR, f"{label}.csv")
    with open(cfg, "w") as f:
        f.write(config_for(seed))
    cmd = [
        "cargo", "run", "--release", "--example", EXAMPLE,
        "--no-default-features", "--features", "precision-double", "--", cfg, out,
    ]
    subprocess.run(cmd, cwd=REPO_ROOT, check=True)
    return out


def numeric_values(path):
    vals = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            for key in ("x", "y", "z", "w"):
                if row[key]:
                    vals.append(float(int(row[key])))
    return vals


def max_delta(path_a, path_b):
    a = numeric_values(path_a)
    b = numeric_values(path_b)
    if len(a) != len(b):
        return math.inf
    return max(abs(x - y) for x, y in zip(a, b))


def graph(same_delta, changed_delta):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(PLOT_DIR, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    labels = ["same seed\nrun A vs B", "changed seed\nA vs C"]
    values = [same_delta, changed_delta]
    colors = ["#2f6f73", "#b45f3c"]
    ax.bar(labels, values, color=colors)
    ax.axhline(0.0, color="black", linewidth=1.0)
    ax.set_ylabel("max absolute difference in fingerprint bits")
    ax.set_title("Config-level clump insertion determinism")
    ax.text(0, same_delta, "byte-identical", ha="center", va="bottom")
    ax.text(1, changed_delta, "diverges", ha="center", va="bottom")
    ax.set_ylim(0, max(1.0, changed_delta * 1.15))
    fig.tight_layout()
    fig.savefig(PLOT, dpi=160)


def main():
    if len(sys.argv) == 1 or sys.argv[1] == "all":
        shutil.rmtree(DATA_DIR, ignore_errors=True)
        a = run_case("same_a", 20260705)
        b = run_case("same_b", 20260705)
        c = run_case("changed", 20260706)
    elif sys.argv[1] == "graph":
        a = os.path.join(DATA_DIR, "same_a.csv")
        b = os.path.join(DATA_DIR, "same_b.csv")
        c = os.path.join(DATA_DIR, "changed.csv")
    else:
        raise SystemExit("usage: sweep.py [all|graph]")

    same_bytes = open(a, "rb").read() == open(b, "rb").read()
    same_delta = max_delta(a, b)
    changed_delta = max_delta(a, c)
    ok = same_bytes and same_delta == 0.0 and math.isfinite(changed_delta) and changed_delta > 0.0
    graph(same_delta, changed_delta)
    print(f"same_seed_byte_identical={same_bytes}")
    print(f"same_seed_max_delta={same_delta:.0f}")
    print(f"changed_seed_max_delta={changed_delta:.0f}")
    print("ALL CHECKS PASSED" if ok else "VALIDATION FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
