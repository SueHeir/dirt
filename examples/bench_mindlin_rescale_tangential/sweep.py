#!/usr/bin/env python3
"""
Mindlin unloading-rescale tangential benchmark driver.

The Rust recorder loads a tangential contact at fixed peak normal overlap, then
unloads the normal overlap with v_t = 0. With e = 1 and a high Coulomb cap, the
documented LAMMPS recurrences are exact:

  history:                 xi_n = xi_{n-1}
  mindlin_rescale:         xi_n = xi_{n-1} * a_n/a_{n-1}  on unloading
  mindlin_rescale/force:   Fte_n = Fte_{n-1} * a_n/a_{n-1} on unloading
  linear_nohistory:        no stored tangential history

This script checks DIRT against those independent recurrences and plots the
unloading force against contact radius.
"""

import os
import sys
import tomllib
import subprocess
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_mindlin_rescale_tangential"
CONFIG = os.path.join(SCRIPT_DIR, "config.toml")
CSV = os.path.join(SCRIPT_DIR, "data", "mindlin_rescale_tangential_results.csv")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
PLOT = os.path.join(PLOT_DIR, "mindlin_rescale_unload.png")


def build_and_run():
    print("[build] cargo build --release --example", EXAMPLE)
    subprocess.run(
        [
            "cargo",
            "build",
            "--release",
            "--example",
            EXAMPLE,
            "--no-default-features",
            "--features",
            "precision-double",
        ],
        cwd=REPO_ROOT,
        check=True,
    )
    binary = os.path.join(REPO_ROOT, "target", "release", "examples", EXAMPLE)
    print("[run]", binary)
    subprocess.run([binary, CONFIG], cwd=REPO_ROOT, check=True)


def load():
    rows = np.genfromtxt(CSV, delimiter=",", names=True, dtype=None, encoding="utf-8")
    with open(CONFIG, "rb") as f:
        cfg = tomllib.load(f)
    return rows, cfg


def split(rows, model):
    m = rows["model"] == model
    return rows[m]


def effective_shear_modulus_identical(e, nu):
    return 1.0 / (4.0 * (2.0 - nu) * (1.0 + nu) / e)


def expected_force(model, cfg):
    sc = cfg["scenario"]
    mat = cfg["material"]
    r_eff = float(sc["radius"]) / 2.0
    g = effective_shear_modulus_identical(float(mat["youngs_mod"]), float(mat["poisson_ratio"]))
    dt = float(sc["dt"])
    vload = float(sc["tangential_velocity"])
    load_steps = int(sc["load_steps"])
    unload_steps = int(sc["unload_steps"])
    peak = float(sc["peak_overlap"])
    final = float(sc["final_overlap"])

    force_history = model == "mindlin_rescale/force"
    displacement_history = model in ("history", "mindlin_rescale")
    rescale = model in ("mindlin_rescale", "mindlin_rescale/force")

    h = 0.0
    prev_a = 0.0
    values = []
    for step in range(load_steps + unload_steps):
        if step < load_steps:
            overlap = peak
            vt = vload
            phase = "load"
        else:
            k = step - load_steps
            frac = k / (unload_steps - 1) if unload_steps > 1 else 1.0
            overlap = peak + frac * (final - peak)
            vt = 0.0
            phase = "unload"

        a = (r_eff * overlap) ** 0.5
        kt = 8.0 * g * a
        if model == "linear_nohistory":
            h = 0.0
            ft = 0.0
        else:
            if rescale and prev_a > 0.0 and a < prev_a:
                h *= a / prev_a
            if force_history:
                h += kt * vt * dt
                ft = h
            elif displacement_history:
                h += vt * dt
                ft = kt * h
            else:
                raise ValueError(model)
        prev_a = a
        values.append((phase, overlap, a, ft, abs(h)))
    return np.array(
        values,
        dtype=[("phase", "U8"), ("overlap", "f8"), ("a", "f8"), ("ft", "f8"), ("hmag", "f8")],
    )


def validate():
    rows, cfg = load()
    checks = []

    def check(name, ok, detail):
        checks.append(bool(ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")

    for model in ("history", "mindlin_rescale", "mindlin_rescale/force", "linear_nohistory"):
        got = split(rows, model)
        exp = expected_force(model, cfg)
        max_ft = max(1.0, float(np.max(np.abs(exp["ft"]))))
        err = float(np.max(np.abs(got["ft"].astype(float) - exp["ft"])))
        herr = float(np.max(np.abs(got["hmag"].astype(float) - exp["hmag"])))
        hscale = max(1.0, float(np.max(np.abs(exp["hmag"]))))
        check(
            f"{model} matches documented recurrence",
            err < 1e-9 * max_ft and herr < 1e-9 * hscale,
            f"max |Ft-ref|={err:.3e} N, max |history-ref|={herr:.3e}",
        )

    hist = split(rows, "history")
    rescale = split(rows, "mindlin_rescale")
    nohist = split(rows, "linear_nohistory")
    force = split(rows, "mindlin_rescale/force")
    unload = hist["phase"] == "unload"
    first_unload = np.flatnonzero(unload)[0]
    last_unload = np.flatnonzero(unload)[-1]

    a_peak = float(hist["a"][first_unload])
    a_final = float(hist["a"][last_unload])
    expected_ratio = a_final / a_peak
    measured_ratio = abs(float(rescale["ft"][last_unload] / hist["ft"][last_unload]))
    check(
        "mindlin_rescale unloads quadratically relative to history",
        abs(measured_ratio - expected_ratio) < 2e-10,
        f"Ft_rescale/Ft_history={measured_ratio:.6f}, a_final/a_peak={expected_ratio:.6f}",
    )

    check(
        "linear_nohistory remains zero during v_t=0 unload",
        float(np.max(np.abs(nohist["ft"][unload]))) < 1e-12,
        f"max unload |Ft|={float(np.max(np.abs(nohist['ft'][unload]))):.3e} N",
    )

    force_drop = abs(float(force["hmag"][last_unload] / force["hmag"][first_unload]))
    check(
        "mindlin_rescale/force gates elastic-force history on unload",
        abs(force_drop - expected_ratio) < 2e-10,
        f"|Fte_final|/|Fte_peak|={force_drop:.6f}, a_final/a_peak={expected_ratio:.6f}",
    )

    print(f"\n{sum(checks)}/{len(checks)} checks passed")
    ok = all(checks)
    print("ALL CHECKS PASSED" if ok else "CHECKS FAILED")
    return ok


def plot():
    rows, cfg = load()
    os.makedirs(PLOT_DIR, exist_ok=True)
    fig, (ax, ax_resid) = plt.subplots(
        2,
        1,
        figsize=(8.0, 7.0),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )
    colors = {
        "history": "#4c78a8",
        "mindlin_rescale": "#d62728",
        "mindlin_rescale/force": "#2ca02c",
        "linear_nohistory": "#666666",
    }
    for model, color in colors.items():
        got = split(rows, model)
        exp = expected_force(model, cfg)
        unload = got["phase"] == "unload"
        ax.plot(exp["a"][unload], exp["ft"][unload], color=color, linewidth=2.0, label=f"{model} ref")
        ax.scatter(got["a"][unload], got["ft"][unload], s=12, color=color, alpha=0.7, edgecolors="none")
        resid = got["ft"][unload].astype(float) - exp["ft"][unload]
        ax_resid.scatter(got["a"][unload], resid, s=10, color=color, alpha=0.7, edgecolors="none", label=model)

    ax.set_ylabel("Tangential force $F_t$ [N]")
    ax.set_title("Mindlin unloading rescale against documented recurrence")
    ax.grid(True, color="#dddddd", linewidth=0.8)
    ax.legend(loc="best", frameon=False, fontsize=8)

    ax_resid.axhline(0.0, color="#777777", linewidth=0.8)
    ax_resid.set_xlabel("Contact radius $a = \\sqrt{R^*\\delta}$ [m]")
    ax_resid.set_ylabel("DIRT - ref [N]")
    ax_resid.grid(True, color="#dddddd", linewidth=0.8)
    fig.tight_layout()
    fig.savefig(PLOT, dpi=180)
    plt.close(fig)
    print(f"[plot] {PLOT}")


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    if cmd in ("all", "start"):
        build_and_run()
    ok = validate()
    if ok:
        plot()
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
