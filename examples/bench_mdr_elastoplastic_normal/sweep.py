#!/usr/bin/env python3
"""
MDR elastic-plastic normal-contact benchmark.

Runs DIRT's selectable `contact_model = "mdr"` for a prescribed loading and
unloading path, then gates the measured normal force against an independent
Python implementation of the same documented one-dimensional MDR reference used
for this DIRT subset. The plot shows the force-displacement hysteresis and the
tolerance band used by the gate.
"""

import csv
import math
import os
import subprocess
import sys

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
CONFIG = os.path.join(SCRIPT_DIR, "config.toml")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
TRACE = os.path.join(DATA_DIR, "mdr_trace.csv")
REF = os.path.join(DATA_DIR, "mdr_reference.csv")
PLOT = os.path.join(PLOT_DIR, "mdr_force_trace.png")

REL_TOL = 2.0e-10
ABS_TOL = 1.0e-8


def load_config():
    with open(CONFIG, "rb") as f:
        return tomllib.load(f)


def nonadhesive_force(delta, a_shape, e_eff, b_shape):
    if delta <= 0.0 or a_shape <= 0.0 or b_shape <= 0.0:
        return 0.0
    x = max(0.0, min(1.0, delta / a_shape))
    return 0.25 * e_eff * a_shape * b_shape * (
        math.acos(1.0 - 2.0 * x) - (2.0 - 4.0 * x) * math.sqrt(max(0.0, x - x * x))
    )


def reference_trace(cfg):
    r_eff = cfg["radius"] / 2.0
    e_eff = 1.0 / (2.0 * (1.0 - cfg["poisson_ratio"] ** 2) / cfg["youngs_mod"])
    gamma = cfg["surface_energy"]
    y = cfg["yield_stress"]
    n = max(2, int(cfg["points_per_leg"]))
    max_delta = cfg["max_overlap"]
    delta_max = 0.0
    yielded = False
    delta_y = 0.0
    a_adh = 0.0
    rows = []
    for phase, values in (
        ("loading", [i * max_delta / (n - 1) for i in range(n)]),
        ("unloading", [(1.0 - i / (n - 1)) * max_delta for i in range(n)]),
    ):
        for delta in values:
            delta_max = max(delta_max, delta)
            if (not yielded) and y > 0.0 and delta > 0.0:
                p_hertz = 4.0 * e_eff * math.sqrt(delta) / (3.0 * math.pi * math.sqrt(r_eff))
                p_y = y * (1.75 * math.exp(-4.4 * delta_max / r_eff) + 1.0)
                if p_hertz > p_y:
                    yielded = True
                    delta_y = delta
            if yielded:
                c_a = math.pi * (delta_y * delta_y - delta_y * r_eff)
                a_max = math.sqrt(max(1.0e-30, 2.0 * delta_max * r_eff - delta_max * delta_max + c_a / math.pi))
                p_y = y * (1.75 * math.exp(-4.4 * delta_max / r_eff) + 1.0)
                a_shape = max(1.0e-30, 4.0 * p_y / e_eff * a_max)
                b_shape = max(1.0e-30, 2.0 * a_max)
                delta_e_max = 0.5 * a_shape
                plastic_offset = max(0.0, delta_max - delta_e_max)
                delta_e = max(0.0, min(delta - plastic_offset, delta_e_max))
            else:
                a_shape = 4.0 * r_eff
                b_shape = 2.0 * r_eff
                delta_e = delta
            force = nonadhesive_force(delta_e, a_shape, e_eff, b_shape)
            if delta_e > 0.0:
                a_contact = b_shape * math.sqrt(max(0.0, a_shape - delta_e)) * math.sqrt(delta_e) / a_shape
            else:
                a_contact = 0.0
            if gamma > 0.0:
                if delta >= delta_max or a_contact >= a_adh:
                    a_adh = 0.99 * a_contact
                elif a_adh > 0.0:
                    force -= 2.0 * math.pi * gamma * a_adh
            rows.append({"phase": phase, "delta": delta, "force": force})
    return rows


def run():
    os.makedirs(DATA_DIR, exist_ok=True)
    cmd = [
        "cargo",
        "run",
        "--release",
        "--example",
        "bench_mdr_elastoplastic_normal",
        "--no-default-features",
        "--features",
        "precision-double",
        "--",
        CONFIG,
    ]
    subprocess.run(cmd, cwd=REPO_ROOT, check=True)


def graph():
    cfg = load_config()
    ref = reference_trace(cfg)
    with open(REF, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["phase", "delta", "force"], lineterminator="\n")
        w.writeheader()
        w.writerows(ref)
    with open(TRACE, newline="") as f:
        actual = list(csv.DictReader(f))
    if len(actual) != len(ref):
        raise SystemExit(f"FAIL: trace length {len(actual)} != reference {len(ref)}")
    max_err = 0.0
    for got, exp in zip(actual, ref):
        fg = float(got["force"])
        fe = exp["force"]
        err = abs(fg - fe) / max(abs(fe), ABS_TOL)
        max_err = max(max_err, err)
        if abs(fg - fe) > ABS_TOL and err > REL_TOL:
            raise SystemExit(
                f"FAIL: {got['phase']} delta={float(got['delta']):.6e} "
                f"force={fg:.12e} ref={fe:.12e} rel_err={err:.3e}"
            )

    os.makedirs(PLOT_DIR, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for phase, marker in [("loading", "o"), ("unloading", "s")]:
        xs = [float(r["delta"]) for r in actual if r["phase"] == phase]
        ys = [float(r["force"]) for r in actual if r["phase"] == phase]
        yr = [r["force"] for r in ref if r["phase"] == phase]
        ax.plot(xs, yr, "-", label=f"{phase} reference")
        ax.plot(xs, ys, marker, ms=3, linestyle="none", label=f"{phase} DIRT")
        band = [ABS_TOL + REL_TOL * max(abs(v), ABS_TOL) for v in yr]
        ax.fill_between(xs, [v - b for v, b in zip(yr, band)], [v + b for v, b in zip(yr, band)], alpha=0.16)
    ax.set_xlabel("overlap delta [m]")
    ax.set_ylabel("normal force [N]")
    ax.set_title("MDR elastic-plastic normal response: DIRT vs analytic reference")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(PLOT, dpi=180)
    print(f"PASS: MDR trace max relative error {max_err:.3e} (tol {REL_TOL:.1e}); plot {PLOT}")


def main():
    action = sys.argv[1] if len(sys.argv) > 1 else "all"
    if action in ("run", "all"):
        run()
    if action in ("graph", "all"):
        graph()


if __name__ == "__main__":
    main()
