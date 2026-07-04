#!/usr/bin/env python3
"""
History-free (`linear_nohistory`) tangential contact model benchmark driver.

Two glass spheres are held at a FIXED normal overlap (so F_n is constant) and
driven with a prescribed, reversing relative TANGENTIAL velocity — a triangle
"load -> unload -> reverse" path (0 -> +V -> -V -> +V -> -V). The Rust recorder
(`main.rs`) evaluates the DIRT pair contact force each step for BOTH tangential
models and writes (model, step, v_t, F_n, F_t, |xi|).

We validate that DIRT's new `linear_nohistory` model reproduces the documented
LAMMPS `pair_granular` velocity-Coulomb law (doc/src/pair_granular.rst,
"tangential linear_nohistory"):

    F_t = -min(mu * F_n, eta_t * |v_t|) * t_hat ,   t_hat = v_t / |v_t|

with ZERO accumulated displacement, and that it is genuinely DISTINCT from the
history-based Mindlin path (which accumulates xi and retains an elastic force at
the v_t = 0 crossings).

Independent reference anchors (not self-consistency):
  * mu is the sliding-friction coefficient from the input deck (config.toml).
  * The Coulomb cap |F_t| = mu * |F_n| is the documented LAMMPS critical force.
  * The velocity-Coulomb shape min(mu*F_n, eta_t*|v_t|) is the documented form;
    eta_t is identified from the sub-cap viscous branch (a straight line through
    the origin — the signature of a history-FREE law) and the full min() shape is
    then verified against every recorded point.

Commands:
    python3 examples/bench_nohistory_tangential/sweep.py            # build + run + validate
    python3 examples/bench_nohistory_tangential/sweep.py graph      # validate existing CSV only

Exit code 0 = PASS, 1 = FAIL.
"""

import os
import sys
import tomllib
import subprocess
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_nohistory_tangential"
CONFIG = os.path.join(SCRIPT_DIR, "config.toml")
CSV = os.path.join(SCRIPT_DIR, "data", "nohistory_tangential_results.csv")


def build_and_run():
    print("[build] cargo build --release --example", EXAMPLE)
    subprocess.run(
        ["cargo", "build", "--release", "--example", EXAMPLE,
         "--no-default-features", "--features", "precision-double"],
        cwd=REPO_ROOT, check=True,
    )
    binary = os.path.join(REPO_ROOT, "target", "release", "examples", EXAMPLE)
    print("[run]", binary)
    subprocess.run([binary, CONFIG], cwd=REPO_ROOT, check=True)


def load():
    rows = np.genfromtxt(CSV, delimiter=",", names=True, dtype=None, encoding="utf-8")
    with open(CONFIG, "rb") as f:
        cfg = tomllib.load(f)
    mu = float(cfg["material"]["friction"])
    return rows, mu


def split(rows, model):
    m = rows["model"] == model
    return (rows["vt"][m].astype(float), rows["fn"][m].astype(float),
            rows["ft"][m].astype(float), rows["xi"][m].astype(float))


def validate():
    rows, mu = load()
    vt_nh, fn_nh, ft_nh, xi_nh = split(rows, "linear_nohistory")
    vt_h, fn_h, ft_h, xi_h = split(rows, "history")

    if len(vt_nh) == 0 or len(vt_h) == 0:
        print("FATAL: missing rows for one of the models"); return False

    mu_fn_nh = mu * np.abs(fn_nh)          # documented Coulomb cap, per row
    cap = np.median(mu_fn_nh)              # F_n is held constant by construction
    zero_v = np.abs(vt_nh) < 1e-9          # v_t = 0 crossings (exact leg boundaries)
    capped = np.abs(ft_nh) > 0.999 * mu_fn_nh
    subcap = (~capped) & (np.abs(vt_nh) > 1e-6)

    checks = []

    def check(name, ok, detail):
        checks.append(ok)
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")

    # ── linear_nohistory: the documented velocity-Coulomb law ──────────────────
    # C1 — zero accumulated displacement.
    max_xi = float(np.max(np.abs(xi_nh)))
    check("nohistory zero displacement accumulation (|xi| == 0)",
          max_xi < 1e-12, f"max|xi| = {max_xi:.3e} (tol 1e-12)")

    # C2 — history-free memory: F_t == 0 whenever v_t == 0 (no elastic term).
    max_ft_at_zero = float(np.max(np.abs(ft_nh[zero_v]))) if zero_v.any() else np.inf
    check("nohistory F_t == 0 at every v_t = 0 crossing",
          max_ft_at_zero < 1e-9, f"max|F_t(v_t=0)| = {max_ft_at_zero:.3e} N (tol 1e-9)")

    # C3 — viscous branch is linear THROUGH THE ORIGIN (no elastic offset).
    #      Identify eta_t (= slope) from the sub-cap points.
    A = np.vstack([vt_nh[subcap], np.ones(subcap.sum())]).T
    (slope, intercept), *_ = np.linalg.lstsq(A, ft_nh[subcap], rcond=None)
    eta_t = slope
    resid = ft_nh[subcap] - slope * vt_nh[subcap]
    rms = float(np.sqrt(np.mean(resid**2)))
    check("nohistory viscous branch linear through origin (F_t = eta_t v_t)",
          abs(intercept) < 1e-4 * cap and rms < 1e-4 * cap,
          f"eta_t = {eta_t:.4f} N.s/m, intercept = {intercept:.3e} N, rms = {rms:.3e} N")

    # C4 — full documented shape F_t = sign(v_t) * min(mu F_n, eta_t |v_t|).
    ft_pred = np.sign(vt_nh) * np.minimum(mu_fn_nh, eta_t * np.abs(vt_nh))
    max_err = float(np.max(np.abs(ft_nh - ft_pred)))
    check("nohistory matches documented min(mu F_n, eta_t |v_t|) at every step",
          max_err < 1e-3 * cap, f"max|F_t - F_t_lammps| = {max_err:.3e} N ({max_err/cap:.2e}*cap, tol 1e-3)")

    # C5 — anti-trivial: BOTH regimes (Coulomb-capped sliding AND viscous) present.
    check("nohistory exercises both sliding (capped) and viscous regimes",
          capped.sum() > 0 and subcap.sum() > 10,
          f"{int(capped.sum())} capped rows, {int(subcap.sum())} viscous rows")

    # C6 — sliding Coulomb magnitude equals mu*|F_n| exactly (independent mu).
    #      "Truly sliding" = the documented min() selects the Coulomb branch,
    #      i.e. eta_t|v_t| exceeds the cap (exclude the thin switch-over band).
    truly_capped = eta_t * np.abs(vt_nh) > mu_fn_nh * 1.001
    if truly_capped.any():
        cap_err = float(np.max(np.abs(np.abs(ft_nh[truly_capped]) - mu_fn_nh[truly_capped])))
    else:
        cap_err = np.inf
    check("nohistory sliding force equals mu*|F_n| (documented Coulomb cap)",
          cap_err < 1e-6 * cap, f"max||F_t| - mu|F_n|| = {cap_err:.3e} N over "
          f"{int(truly_capped.sum())} sliding rows (tol 1e-6)")

    # ── history model: genuinely DISTINCT (this is the point of the goal) ──────
    # C7 — history accumulates a nonzero tangential displacement.
    max_xi_h = float(np.max(np.abs(xi_h)))
    check("history model accumulates displacement (|xi| > 0)",
          max_xi_h > 1e-8, f"max|xi_history| = {max_xi_h:.3e} m")

    # C8 — the distinguishing behavior: at v_t = 0 the history model retains a
    #      large elastic force while the history-free model is exactly zero.
    zero_v_h = np.abs(vt_h) < 1e-9
    ft_h_at_zero = float(np.max(np.abs(ft_h[zero_v_h]))) if zero_v_h.any() else 0.0
    check("history retains elastic force at v_t = 0 (nohistory does not)",
          ft_h_at_zero > 0.5 * cap and max_ft_at_zero < 1e-9,
          f"|F_t^history(v_t=0)| = {ft_h_at_zero:.3f} N  vs  "
          f"|F_t^nohistory(v_t=0)| = {max_ft_at_zero:.3e} N")

    npass = sum(checks)
    ntot = len(checks)
    print(f"\n{npass}/{ntot} checks passed")
    ok = npass == ntot
    print("ALL CHECKS PASSED" if ok else "CHECKS FAILED")
    return ok


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    if cmd in ("all", "start"):
        build_and_run()
    ok = validate()
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
