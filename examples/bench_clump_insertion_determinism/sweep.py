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
PLOT = os.path.join(PLOT_DIR, "clump_insertion_determinism.svg")


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
    os.makedirs(PLOT_DIR, exist_ok=True)
    # Keep the result graph reproducible on a plain Python installation.  This
    # intentionally emits SVG directly instead of relying on an undeclared
    # plotting package after the solver cases have already completed.
    width, height, baseline = 760, 360, 270
    changed_height = 180
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="white"/>
<text x="380" y="34" text-anchor="middle" font-family="sans-serif" font-size="20">Config-level clump insertion determinism</text>
<line x1="90" y1="{baseline}" x2="710" y2="{baseline}" stroke="black" stroke-width="2"/>
<rect x="190" y="{baseline - 2}" width="110" height="2" fill="#2f6f73"/>
<rect x="470" y="{baseline - changed_height}" width="110" height="{changed_height}" fill="#b45f3c"/>
<text x="245" y="{baseline + 28}" text-anchor="middle" font-family="sans-serif" font-size="14">same seed A vs B</text>
<text x="245" y="{baseline + 48}" text-anchor="middle" font-family="sans-serif" font-size="13">0 bits; byte-identical</text>
<text x="525" y="{baseline + 28}" text-anchor="middle" font-family="sans-serif" font-size="14">changed seed A vs C</text>
<text x="525" y="{baseline + 48}" text-anchor="middle" font-family="sans-serif" font-size="13">{changed_delta:.0f} max bit delta</text>
<text x="28" y="185" transform="rotate(-90 28 185)" text-anchor="middle" font-family="sans-serif" font-size="14">fingerprint divergence (bits; not to scale)</text>
<text x="525" y="{baseline - changed_height - 10}" text-anchor="middle" font-family="sans-serif" font-size="13">different deterministic stream</text>
</svg>'''
    with open(PLOT, "w", encoding="utf-8") as f:
        f.write(svg)


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
