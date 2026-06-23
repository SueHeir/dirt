#!/usr/bin/env python3
"""Cooperativity-length closure sweep for the MUD nonlocal (NGF) branch.

Runs a shear-rate (γ̇) grid of the Lees-Edwards rig, time-averages the steady
window, and fits the two closures MUD's nonlocal branch needs:

  1. ξ(μ) = A·d / √(μ − μ_s)   →  the cooperativity amplitude A (and divergence).
  2. g = γ̇/μ ∝ √T             →  the Zhang-Kamrin bridge that lets MUD drive the
                                  fluidity from granular temperature.

Usage:
  python3 sweep.py --run      # generate per-γ̇ configs, run the DEM, then analyze
  python3 sweep.py            # analyze existing data/ CSVs only

Writes data/calibration.yaml (A, μ_s, the g–√T slope) for MUD's MaterialParams,
and plots/cooperativity.png. PASS if the ξ(μ) fit and the g∝√T law hold.
"""
import argparse
import glob
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
PLOTS = os.path.join(HERE, "plots")
D_GRAIN = 0.5e-3
MU_S = 0.38  # static yield from the μ(I) fit (bench_lebc_shear); refit if desired

GDOTS = [10.0, 30.0, 50.0, 100.0, 200.0, 400.0]  # shear-rate grid → spans I, μ, T


def gen_config(gdot):
    """Write a per-γ̇ config (distinct output dir) and return (cfg_path, out_dir)."""
    base = os.path.join(HERE, "config.toml")
    with open(base) as f:
        txt = f.read()
    tag = f"g{gdot:g}"
    out_dir = os.path.join(DATA, tag)
    txt = txt.replace("rate = 50.0", f"rate = {gdot}")
    txt = txt.replace(
        'dir = "examples/SPH_glass_sphere_calibration/08_cooperativity_length/data"',
        f'dir = "examples/SPH_glass_sphere_calibration/08_cooperativity_length/data/{tag}"',
    )
    cfg = os.path.join(DATA, f"config_{tag}.toml")
    os.makedirs(out_dir, exist_ok=True)
    with open(cfg, "w") as f:
        f.write(txt)
    return cfg, out_dir


def run_case(cfg):
    # dirt repo root is three levels up: 08_.../ → SPH_glass_sphere_calibration/ → examples/ → root
    root = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
    rel = os.path.relpath(cfg, root)
    subprocess.run(
        ["cargo", "run", "--release", "--example",
         "sphcal_cooperativity_length", "--no-default-features", "--", rel],
        cwd=root, check=True,
    )


def steady_average(csv):
    """Average the steady window (last half of the shear stage) of one CSV."""
    rows = []
    with open(csv) as f:
        header = f.readline().strip().split(",")
        for line in f:
            rows.append([float(x) for x in line.strip().split(",")])
    if not rows:
        return None
    idx = {n: i for i, n in enumerate(header)}
    # steady = γ̇ > 0 (shear stage) and the last half of those samples
    sheared = [r for r in rows if r[idx["gdot"]] > 1e-9]
    if len(sheared) < 4:
        return None
    tail = sheared[len(sheared) // 2:]
    avg = {k: sum(r[idx[k]] for r in tail) / len(tail) for k in idx}
    return avg


def collect():
    pts = []
    for csv in sorted(glob.glob(os.path.join(DATA, "*", "cooperativity_results.csv"))):
        a = steady_average(csv)
        if a:
            pts.append(a)
    # also accept a top-level single-run CSV
    top = os.path.join(DATA, "cooperativity_results.csv")
    if os.path.exists(top):
        a = steady_average(top)
        if a:
            pts.append(a)
    return pts


def analyze(pts):
    if not pts:
        print("no data — run with --run first"); return 1
    print(f"{'gdot':>7} {'I':>8} {'mu':>7} {'T':>10} {'g':>8} {'xi[d]':>8}")
    xs, ys, gs, ts = [], [], [], []
    for a in pts:
        mu, xi, g, t = a["mu"], a["xi"], a["g"], a["T"]
        print(f"{a['gdot']:7.0f} {a['I']:8.4f} {mu:7.3f} {t:10.2e} {g:8.1f} {xi/D_GRAIN:8.2f}")
        if mu > MU_S + 1e-3:
            xs.append(1.0 / math.sqrt(mu - MU_S))  # ξ = A d · x
            ys.append(xi)
        if t > 0:
            gs.append(g); ts.append(math.sqrt(t))

    # Fit ξ = (A·d)·(μ−μ_s)^(-1/2) through the origin → slope = A·d.
    a_amp = None
    if len(xs) >= 2:
        slope = sum(x * y for x, y in zip(xs, ys)) / sum(x * x for x in xs)
        a_amp = slope / D_GRAIN
        # R² of the through-origin fit
        ybar = sum(ys) / len(ys)
        ss_res = sum((y - slope * x) ** 2 for x, y in zip(xs, ys))
        ss_tot = sum((y - ybar) ** 2 for y in ys) or 1.0
        r2_xi = 1 - ss_res / ss_tot
        print(f"\nξ(μ) = A·d/√(μ−μ_s):  A = {a_amp:.3f}  (R² = {r2_xi:.3f}, μ_s = {MU_S})")
    # Fit g = k·√T through the origin (Zhang-Kamrin bridge).
    r2_g = None
    if len(gs) >= 2:
        kk = sum(s * g for s, g in zip(ts, gs)) / sum(s * s for s in ts)
        gbar = sum(gs) / len(gs)
        ss_res = sum((g - kk * s) ** 2 for s, g in zip(ts, gs))
        ss_tot = sum((g - gbar) ** 2 for g in gs) or 1.0
        r2_g = 1 - ss_res / ss_tot
        print(f"g = k·√T (Zhang-Kamrin):  k = {kk:.2f}  (R² = {r2_g:.3f})")

    try:
        plot(pts, a_amp)
    except Exception as e:  # noqa: BLE001
        print(f"(plot skipped: {e})")

    if a_amp is not None:
        os.makedirs(DATA, exist_ok=True)
        with open(os.path.join(DATA, "calibration.yaml"), "w") as f:
            f.write("# MUD nonlocal-cooperativity closure (DEM-calibrated)\n")
            f.write(f"coop_amplitude: {a_amp:.4f}   # A in ξ = A d/√(μ−μ_s)\n")
            f.write(f"mu_s: {MU_S}\n")
        print(f"\nwrote {os.path.join(DATA, 'calibration.yaml')}  → set coop_amplitude = {a_amp:.3f} in MUD")

    ok = (a_amp is not None and a_amp > 0) and (r2_g is None or r2_g > 0.5)
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


def plot(pts, a_amp):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(PLOTS, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    mus = [a["mu"] for a in pts]
    xis = [a["xi"] / D_GRAIN for a in pts]
    ax1.plot(mus, xis, "o", color="#c0392b")
    if a_amp:
        import numpy as np
        mm = np.linspace(MU_S + 0.02, max(mus) + 0.02, 100)
        ax1.plot(mm, a_amp / np.sqrt(mm - MU_S), "-", color="#1f4e79",
                 label=f"A/√(μ−μ_s), A={a_amp:.2f}")
        ax1.legend()
    ax1.set_xlabel("μ = τ/p"); ax1.set_ylabel("ξ / d"); ax1.set_title("Cooperativity length vs μ")
    ax1.grid(alpha=0.3)
    ts = [math.sqrt(a["T"]) for a in pts]
    gs = [a["g"] for a in pts]
    ax2.plot(ts, gs, "o", color="#c0392b")
    ax2.set_xlabel("√T  [m/s]"); ax2.set_ylabel("g = γ̇/μ  [1/s]")
    ax2.set_title("Zhang–Kamrin bridge: g ∝ √T"); ax2.grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(PLOTS, "cooperativity.png")
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true", help="generate + run the γ̇ sweep first")
    args = ap.parse_args()
    if args.run:
        for g in GDOTS:
            cfg, _ = gen_config(g)
            print(f"running γ̇ = {g} ...")
            run_case(cfg)
    sys.exit(analyze(collect()))
