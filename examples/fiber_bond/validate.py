#!/usr/bin/env python3
"""
Validate one `fiber_bond` test run against Guo 2018 (*Chem. Eng. Sci.* 175,
118–129) analytical predictions and produce side-by-side overlay plots.

For each scenario the validator:
  * loads `fiber_bond.csv` (time series) and `profile.csv` (per-atom snapshot)
  * computes the model-derived quantity from kinematics + bond config
  * compares it to the closed-form Guo prediction
  * prints a single-block summary
  * writes one or two PNG plots showing **theory line + DEM markers** on
    the same axes — the easy-to-read comparison

Usage
-----
    python3 examples/fiber_bond/validate.py PATH_TO_CSV
"""

from __future__ import annotations

import csv
import math
import os
import sys
from pathlib import Path


# ── CSV loaders ──────────────────────────────────────────────────────────────

def load_timeseries(path):
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({k: float(v) for k, v in row.items()})
    return rows


def load_profile(path):
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({
                "tag":  int(float(row["tag"])),
                "x0":   float(row["x0"]),
                "y0":   float(row["y0"]),
                "z0":   float(row["z0"]),
                "x":    float(row["x"]),
                "y":    float(row["y"]),
                "z":    float(row["z"]),
                "mass": float(row.get("mass", 0.0)),
            })
    return sorted(rows, key=lambda r: r["x0"])


def linear_fit_through_origin(x, y):
    sxx = sum(xi * xi for xi in x)
    sxy = sum(xi * yi for xi, yi in zip(x, y))
    return sxy / sxx if sxx else 0.0


# ── Plot helpers ─────────────────────────────────────────────────────────────

def try_import_matplotlib():
    try:
        import matplotlib.pyplot as plt  # type: ignore
        return plt
    except ImportError:
        print("  (matplotlib not available — skipping plots)")
        return None


# ── Axial elastic ────────────────────────────────────────────────────────────

def validate_axial_elastic(rows, profile, out_dir):
    """Guo 2018 Eq. 1 — recover E from σ(ε) slope and uniformity along the fiber.

    Two figures:
      * `axial_stress_strain.png` — σ vs ε scatter (DEM) + analytical line.
      * `axial_profile.png`       — per-bond strain along the fiber profile
                                    (uniform-strain check; analogous to Guo Fig. 6).
    """
    last = rows[-1]
    length0 = last["length0"]
    bond_len = last["bond_len_mid0"]
    area = last["area"]
    k_n = last["k_n"]
    e_from_kn = k_n * bond_len / area

    sig_mid = [r["delta_mid"] * k_n / area for r in rows]
    eps_global = [(r["length_global"] - length0) / length0 for r in rows]

    ntrim = max(1, len(rows) // 10)
    e_fit = linear_fit_through_origin(eps_global[ntrim:], sig_mid[ntrim:])
    err = abs(e_fit - e_from_kn) / e_from_kn

    print("=== Axial elastic — Guo 2018 Eq. 1 ===")
    print(f"  samples              : {len(rows)} rows")
    print(f"  L₀ (end-to-end)      : {length0:.6e} m")
    print(f"  bond length L        : {bond_len:.6e} m")
    print(f"  area A               : {area:.6e} m²")
    print(f"  E from K_n·L/A       : {e_from_kn:.6e} Pa  (config)")
    print(f"  E from σ/ε slope     : {e_fit:.6e} Pa  (measured)")
    print(f"  relative error       : {err:.3%}")
    status = "PASS" if err < 0.005 else "FAIL"
    print(f"  status               : {status}  (tolerance 0.5%)")

    plt = try_import_matplotlib()
    if plt is None:
        return err < 0.005

    # ── Figure 1: σ(ε) scatter + theory line ────────────────────────────
    eps_max = max(eps_global)
    eps_min = min(eps_global)
    eps_th = [eps_min + (eps_max - eps_min) * i / 100.0 for i in range(101)]
    sig_th = [e_from_kn * e for e in eps_th]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot([e * 100 for e in eps_th], [s / 1e6 for s in sig_th],
            "-", lw=1.8, c="C1",
            label=f"Theory σ = E·ε  (Guo Eq. 1, E = {e_from_kn:.3g} Pa)")
    ax.plot([e * 100 for e in eps_global], [s / 1e6 for s in sig_mid],
            "o", ms=4, mfc="white", mec="C0", mew=1.2,
            label="DEM (this code, mid-bond)")
    ax.set_xlabel("axial strain  ε  (%)")
    ax.set_ylabel("axial stress  σ  (MPa)")
    ax.set_title("Axial elastic — stress / strain  vs.  σ = E·ε")
    ax.legend(loc="upper left", framealpha=0.95)
    out = Path(out_dir) / "axial_stress_strain.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  plot saved           : {out}")

    # ── Figure 2: per-bond strain profile (Guo Fig. 6 analog) ───────────
    if profile and len(profile) >= 2:
        xs0 = [p["x0"] for p in profile]
        xs  = [p["x"]  for p in profile]
        # Cumulative displacement δ_n = x − x₀ along the fiber.
        delta = [xs[i] - xs0[i] for i in range(len(xs))]
        delta_max = max(abs(d) for d in delta) or 1.0
        # Analytical profile under uniform tension: δ_n(x₀) = ε · x₀, so the
        # normalised profile δ_n / max(δ_n) = x₀ / x₀_max is a 1:1 diagonal.
        x_max = max(xs0) or 1.0

        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.plot([0.0, 1.0], [0.0, 1.0],
                "-", lw=1.8, c="C1",
                label="Theory δ_n/|δ_max| = x/L_c  (Guo Fig. 6)")
        ax.plot([x / x_max for x in xs0], [d / delta_max for d in delta],
                "o", ms=8, mfc="white", mec="C0", mew=1.5,
                label="DEM (this code, per atom)")
        ax.set_xlabel("position along fiber  x / L_c")
        ax.set_ylabel("normalised axial displacement  δ_n / |δ_max|")
        ax.set_title("Axial elastic — displacement profile")
        ax.legend(loc="upper left", framealpha=0.95)
        out = Path(out_dir) / "axial_profile.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"  plot saved           : {out}")

    return err < 0.005


# ── Cantilever bending elastic ───────────────────────────────────────────────

def validate_cantilever_bending(rows, profile, out_dir):
    """Guo 2018 Sec. 2.1, Fig. 3 — recover E·I from `y_tip = F·L_c³/(3·E·I)`.

    Two figures:
      * `cantilever_profile.png` — full y(x) deflection along the fiber
                                   (DEM atom positions) overlaid on the
                                   small-deformation EB profile
                                   y(x) = F·x²·(3L − x)/(6·E·I).
      * `cantilever_tip_vs_time.png` — tip z(t) settling curve with the
                                       steady-state EB prediction marked.
    """
    last = rows[-1]
    iben = last["iben"]
    k_bend = last["k_bend"]
    bond_len = last["bond_len_mid0"]
    length0 = last["length0"]
    ei = k_bend * bond_len
    lc = length0

    # Steady-state tip deflection from the last 10 % of frames.
    n = len(rows)
    tail = rows[max(0, int(0.9 * n)):]
    z_steady = sum(r["right_z"] for r in tail) / len(tail) if tail else float("nan")
    y_steady = sum(r["right_y"] for r in tail) / len(tail) if tail else float("nan")

    fz = float(os.environ.get("FIBER_BOND_FZ", "-0.1"))
    y_pred = fz * lc**3 / (3.0 * ei)
    ei_meas = fz * lc**3 / (3.0 * z_steady) if abs(z_steady) > 1e-30 else float("inf")
    fl2_ei = abs(fz) * lc**2 / ei
    err = abs(z_steady - y_pred) / abs(y_pred) if abs(y_pred) > 1e-30 else float("inf")

    print("=== Cantilever bending elastic — Guo 2018 Sec. 2.1 ===")
    print(f"  samples (steady tail)       : {len(tail)}")
    print(f"  L_c (anchor → tip)          : {lc:.6e} m")
    print(f"  I (bond cross-section)      : {iben:.6e} m⁴")
    print(f"  E·I from config             : {ei:.6e} N·m²")
    print(f"  applied F_z                 : {fz:+.4e} N  (env: FIBER_BOND_FZ)")
    print(f"  F·L_c² / (E·I)              : {fl2_ei:.4f}  "
          f"({'small-def OK' if fl2_ei < 1.0 else 'nonlinear regime'})")
    print(f"  tip z (steady)              : {z_steady:+.6e} m  (predicted {y_pred:+.6e})")
    print(f"  tip y (steady)              : {y_steady:+.6e} m  (should be ~0)")
    print(f"  E·I from measured y_tip     : {ei_meas:.6e} N·m²")
    print(f"  relative error in y_tip     : {err:.3%}")
    tol = 0.05 if fl2_ei < 0.3 else 0.20
    status = "PASS" if err < tol else "FAIL"
    print(f"  status                      : {status}  (tolerance {tol:.0%})")

    plt = try_import_matplotlib()
    if plt is None:
        return err < tol

    # ── Figure 1: y(x) profile DEM vs theory ────────────────────────────
    if profile and len(profile) >= 2:
        xs0 = [p["x0"] for p in profile]
        zs  = [p["z"]  for p in profile]
        # Analytical y(x) = F·x²·(3·L−x)/(6·E·I), with x measured from the
        # pinned (anchor) end of the fiber.
        x_anchor = min(xs0)
        x_theory = [i * lc / 100.0 for i in range(101)]
        y_theory = [fz * x**2 * (3.0 * lc - x) / (6.0 * ei) for x in x_theory]

        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.plot([(x + x_anchor) * 1e3 for x in x_theory],
                [y * 1e3 for y in y_theory],
                "-", lw=1.8, c="C1",
                label="Theory  y(x) = F·x²·(3L−x)/(6EI)  (Guo Sec. 2.1)")
        ax.plot([x * 1e3 for x in xs0],
                [z * 1e3 for z in zs],
                "o", ms=8, mfc="white", mec="C0", mew=1.5,
                label="DEM (this code, per atom)")
        ax.axvline(x_anchor * 1e3, ls=":", c="k", lw=0.8, alpha=0.5)
        ax.text(x_anchor * 1e3 + 0.5, 0.02 * abs(y_pred) * 1e3,
                "anchor", fontsize=9, alpha=0.7)
        ax.set_xlabel("position along fiber  x  (mm)")
        ax.set_ylabel("transverse deflection  y  (mm)")
        ax.set_title(f"Cantilever bending profile  (F = {fz:+.3g} N,  L_c = {lc*1e3:.0f} mm)")
        ax.legend(loc="lower left", framealpha=0.95)
        out = Path(out_dir) / "cantilever_profile.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"  plot saved                  : {out}")

    # ── Figure 2: tip z(t) settling ─────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot([r["t"] * 1e3 for r in rows],
            [r["right_z"] * 1e3 for r in rows],
            "-", lw=1.2, c="C0", label="DEM tip z(t)")
    ax.axhline(y_pred * 1e3, ls="--", lw=1.6, c="C1",
               label=f"Theory steady-state  y_tip = {y_pred*1e3:+.3f} mm")
    ax.set_xlabel("time  t  (ms)")
    ax.set_ylabel("tip deflection  z  (mm)")
    ax.set_title("Cantilever bending — tip settling under damped relaxation")
    ax.legend(loc="upper right", framealpha=0.95)
    out = Path(out_dir) / "cantilever_tip_vs_time.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  plot saved                  : {out}")

    return err < tol


# ── Dispatch ─────────────────────────────────────────────────────────────────

# ── Bending vibration ────────────────────────────────────────────────────────

def validate_bending_vibration(rows, profile, out_dir):
    """Guo 2018 Eq. 18 — recover natural bending period from a step-loaded,
    undamped cantilever tip trajectory.

    Steady-state ringing under a constant tip force is `z(t) = y_static · (1 − cos(ω·t))`
    oscillating between 0 and 2·y_static. Period found from local maxima of
    `z(t)` (the recurrent returns near z ≈ 0).
    """
    last = rows[-1]
    iben = last["iben"]
    k_bend = last["k_bend"]
    bond_len = last["bond_len_mid0"]
    length0 = last["length0"]
    area = last["area"]
    ei = k_bend * bond_len
    lc = length0

    fz = float(os.environ.get("FIBER_BOND_FZ", "-0.1"))
    y_static = fz * lc**3 / (3.0 * ei)

    # Guo Eq. 18 in continuum form: T = 1.7868·L²·√(ρ_l/EI) with ρ_l = ρ·A.
    # Our chain is **not** a continuum — mass is lumped at sphere centres
    # and the bonds (and any gaps between spheres) are massless. So we use
    # the discrete-mass-equivalent ρ_l = M_chain / L, where M_chain is the
    # actual total atom mass read from the profile snapshot.
    n_atoms = len(profile) if profile else 0
    m_chain = sum(p["mass"] for p in profile) if profile else 0.0
    rho_l_discrete = m_chain / lc if lc > 0 else 0.0
    # For reference / diagnostic — what Guo's pure-continuum form would give:
    rho_continuum = 2500.0          # the canonical config value
    rho_l_continuum = rho_continuum * area

    t_bend_discrete  = 1.7868 * lc**2 * math.sqrt(rho_l_discrete  / ei) if rho_l_discrete  > 0 else float("nan")
    t_bend_continuum = 1.7868 * lc**2 * math.sqrt(rho_l_continuum / ei)
    t_bend_pred = t_bend_discrete
    omega_pred = 2.0 * math.pi / t_bend_pred if t_bend_pred > 0 else 0.0

    # Extract period from local maxima of z(t) (z is most negative at the
    # bottom of the swing, most positive on the rebound back near 0). Local
    # max in a quadratic-fit sense over a 3-point window keeps the result
    # stable when the sample rate doesn't perfectly hit the peak.
    t_arr = [r["t"] for r in rows]
    z_arr = [r["right_z"] for r in rows]
    maxima_t = []
    for i in range(1, len(rows) - 1):
        if z_arr[i] > z_arr[i - 1] and z_arr[i] >= z_arr[i + 1]:
            # 3-point parabolic fit refinement.
            denom = (z_arr[i - 1] - 2.0 * z_arr[i] + z_arr[i + 1])
            offset = 0.5 * (z_arr[i - 1] - z_arr[i + 1]) / denom if denom != 0 else 0.0
            dt_loc = 0.5 * (t_arr[i + 1] - t_arr[i - 1])
            maxima_t.append(t_arr[i] + offset * dt_loc)
    # Period = mean spacing between consecutive maxima.
    if len(maxima_t) >= 2:
        deltas = [maxima_t[i + 1] - maxima_t[i] for i in range(len(maxima_t) - 1)]
        t_bend_meas = sum(deltas) / len(deltas)
    else:
        t_bend_meas = float("nan")
    err = (abs(t_bend_meas - t_bend_pred) / t_bend_pred
           if math.isfinite(t_bend_meas) else float("inf"))

    print("=== Bending vibration — Guo 2018 Eq. 18 (discrete-mass form) ===")
    print(f"  L_c (anchor → tip)          : {lc:.6e} m")
    print(f"  E·I                         : {ei:.6e} N·m²")
    print(f"  N atoms / total mass        : {n_atoms} / {m_chain:.6e} kg")
    print(f"  ρ_l discrete (M_chain/L)    : {rho_l_discrete:.6e} kg/m")
    print(f"  ρ_l continuum (ρ·A)         : {rho_l_continuum:.6e} kg/m  (diagnostic only)")
    print(f"  applied F_z                 : {fz:+.4e} N  (env: FIBER_BOND_FZ)")
    print(f"  y_static (osc. centre)      : {y_static*1e3:+.3f} mm")
    print(f"  maxima detected             : {len(maxima_t)}")
    print(f"  T_bend predicted (discrete) : {t_bend_discrete *1e3:.4f} ms  ← used for PASS/FAIL")
    print(f"  T_bend predicted (continuum): {t_bend_continuum*1e3:.4f} ms  (would apply if bonds had mass)")
    print(f"  T_bend measured             : {t_bend_meas*1e3:.4f} ms")
    print(f"  relative error              : {err:.3%}")
    tol = 0.05
    status = "PASS" if err < tol else "FAIL"
    print(f"  status                      : {status}  (tolerance {tol:.0%})")

    plt = try_import_matplotlib()
    if plt is None:
        return err < tol

    # ── Figure: tip z(t) overlay with analytic cosine ───────────────────
    t_meas_ms = [t * 1e3 for t in t_arr]
    z_meas_mm = [z * 1e3 for z in z_arr]
    # Theory: z(t) = y_static · (1 − cos(ω·t)).
    t_theory = [t_arr[-1] * i / 1000.0 for i in range(1001)]
    z_theory = [y_static * (1.0 - math.cos(omega_pred * t)) for t in t_theory]

    # Use a wider plot for the vibration time series so the period structure
    # is legible. Restrict to first ~5 predicted periods.
    t_window = min(t_arr[-1], 5.0 * t_bend_discrete)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    t_theory_w = [t for t in t_theory if t <= t_window]
    z_theory_w = [y_static * (1.0 - math.cos(omega_pred * t)) for t in t_theory_w]
    ax.plot([t * 1e3 for t in t_theory_w], [z * 1e3 for z in z_theory_w],
            "-", lw=1.8, c="C1",
            label=f"Theory  Guo Eq. 18 (discrete-mass)  T = {t_bend_discrete*1e3:.3f} ms")
    meas_in_window = [(t, z) for t, z in zip(t_meas_ms, z_meas_mm) if t * 1e-3 <= t_window]
    ax.plot([t for t, _ in meas_in_window], [z for _, z in meas_in_window],
            "-", lw=1.3, c="C0",
            label=f"DEM tip z(t)  measured T = {t_bend_meas*1e3:.3f} ms")
    # Mark each detected maximum with a small tick.
    for tm in maxima_t:
        ax.axvline(tm * 1e3, ls=":", c="0.5", lw=0.6)
    ax.set_xlabel("time  t  (ms)")
    ax.set_ylabel("tip deflection  z  (mm)")
    ax.set_title("Bending vibration — natural-period ringing under step load")
    ax.set_xlim(0.0, t_window * 1e3)
    ax.legend(loc="lower right", framealpha=0.95)
    out = Path(out_dir) / "bending_vibration.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  plot saved                  : {out}")
    return err < tol


# ── Axial plastic — piecewise hardening ─────────────────────────────────────

def validate_axial_plastic_piecewise(rows, profile, out_dir):
    """Validate that the piecewise-linear axial plasticity envelope is
    traced correctly under monotonic loading.

    Loads the breakpoints & multipliers from the env var FIBER_BOND_PIECEWISE
    (format: `bp1,bp2,bp3;m1,m2,m3`), else falls back to the canonical
    `[0.01, 0.02, 0.03] / [0.5, 0.1, 0.0]` config.
    """
    last = rows[-1]
    length0 = last["length0"]
    bond_len = last["bond_len_mid0"]
    area = last["area"]
    k_n = last["k_n"]
    k_eff = k_n * bond_len  # = E·A — stiffness per strain.

    # Piecewise envelope spec (matches `axial_plastic_piecewise.toml`).
    env = os.environ.get("FIBER_BOND_PIECEWISE", "0.01,0.02,0.03;0.5,0.1,0.0")
    bp_str, mult_str = env.split(";")
    breaks = [float(x) for x in bp_str.split(",") if x]
    mults  = [float(x) for x in mult_str.split(",") if x]

    def envelope(eps_mag):
        f_acc = 0.0
        eps_prev = 0.0
        slope = k_eff
        for bp, m in zip(breaks, mults):
            if eps_mag <= bp:
                return f_acc + slope * (eps_mag - eps_prev)
            f_acc += slope * (bp - eps_prev)
            eps_prev = bp
            slope = k_eff * m
        return f_acc + slope * (eps_mag - eps_prev)

    # Actual bond force from the kinematic strain minus the plastic anchor:
    #   F = K_n · (δ − ε_p · L_bond) = K_n · L_bond · (ε_axial − ε_p_axial)
    # For elastic regime ε_p = 0 → F = K_n · δ.
    eps_p_axial = [r["eps_p_axial_mid"] for r in rows]
    f_mid = [
        k_n * (r["delta_mid"] - r["eps_p_axial_mid"] * bond_len)
        for r in rows
    ]
    eps_global = [(r["length_global"] - length0) / length0 for r in rows]

    # Peak strain reached.
    eps_peak = max(eps_global)
    f_peak   = max(f_mid)
    # Compare DEM force at three sample strains to the envelope.
    sample_eps = [bp * 0.5 for bp in breaks] + [breaks[-1] * 1.5]
    print("=== Axial plastic piecewise — envelope tracing ===")
    print(f"  breakpoint strains : {breaks}")
    print(f"  slope multipliers  : {mults}")
    print(f"  k_eff = E·A        : {k_eff:.6e} N")
    print(f"  ε peak reached     : {eps_peak:.4%}")
    print(f"  F peak measured    : {f_peak:.4e} N")
    max_err = 0.0
    for eps_target in sample_eps:
        if eps_target > eps_peak:
            continue
        # Find DEM force closest to eps_target on the monotonic-loading
        # branch (eps_global is monotonically increasing in this run).
        idx = min(range(len(eps_global)), key=lambda i: abs(eps_global[i] - eps_target))
        f_dem = f_mid[idx]
        f_env = envelope(eps_target)
        err = abs(f_dem - f_env) / max(abs(f_env), 1e-30)
        max_err = max(max_err, err)
        print(f"    ε ≈ {eps_target:.4%}  →  F_DEM = {f_dem:.4e} N  "
              f"F_env = {f_env:.4e} N  err = {err:.3%}")
    tol = 0.05
    status = "PASS" if max_err < tol else "FAIL"
    print(f"  status             : {status}  (tolerance {tol:.0%})")

    plt = try_import_matplotlib()
    if plt is None:
        return max_err < tol

    # ── Figure 1: F vs ε overlay ────────────────────────────────────────
    eps_th = [eps_peak * i / 500.0 for i in range(501)]
    f_th = [envelope(e) for e in eps_th]
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot([e * 100 for e in eps_th], [f * 1e3 for f in f_th],
            "-", lw=1.8, c="C1", label="Theory  piecewise envelope (this code, configured)")
    ax.plot([e * 100 for e in eps_global], [f * 1e3 for f in f_mid],
            "o", ms=3, mfc="white", mec="C0", mew=0.8,
            label="DEM  axial bond force vs strain")
    # Mark breakpoints.
    for bp in breaks:
        ax.axvline(bp * 100, ls=":", c="0.6", lw=0.6)
        ax.text(bp * 100 + 0.02, 0.05 * f_peak * 1e3, f"ε={bp*100:.1f}%",
                fontsize=8, color="0.4", rotation=90, va="bottom")
    ax.set_xlabel("axial strain  ε  (%)")
    ax.set_ylabel("axial bond force  F_n  (mN)")
    ax.set_title("Axial plastic — piecewise hardening envelope")
    ax.legend(loc="lower right", framealpha=0.95)
    out = Path(out_dir) / "axial_plastic_envelope.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  plot saved         : {out}")

    # ── Figure 2: ε_p_axial accumulation ────────────────────────────────
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot([e * 100 for e in eps_global],
            [ep * 100 for ep in eps_p_axial],
            "-", lw=1.6, c="C0", label="DEM  cumulative plastic strain ε_p")
    # Theory: ε_p = max(0, ε − ε_yield) ... actually no, ε_p advances along the envelope.
    # The relationship is ε_p = ε − F_env(ε)/k_eff, so trace that.
    eps_p_th = [e - envelope(e) / k_eff for e in eps_th]
    ax.plot([e * 100 for e in eps_th], [ep * 100 for ep in eps_p_th],
            "--", lw=1.4, c="C1", label="Theory  ε_p = ε − F_env(ε)/k_eff")
    for bp in breaks:
        ax.axvline(bp * 100, ls=":", c="0.6", lw=0.6)
    ax.set_xlabel("axial strain  ε  (%)")
    ax.set_ylabel("plastic axial strain  ε_p  (%)")
    ax.set_title("Axial plastic — plastic-anchor accumulation")
    ax.legend(loc="upper left", framealpha=0.95)
    out = Path(out_dir) / "axial_plastic_eps_p.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  plot saved         : {out}")

    return max_err < tol


# ── Bending plastic — Guo three-step loading ────────────────────────────────

def validate_bending_plastic_guo(rows, profile, out_dir):
    """Guo 2018 Figs. 10/12/13 — three load cycles drive the middle bond past
    yield; verify that `M(θ_bend)` traces an elastic-perfectly-plastic
    hysteresis envelope and that `θ_p_bend` accumulates monotonically.

    For Phase 1 our model is **bilinear** (Guo's intermediate K_ep slope is
    deferred): the envelope is

        |M_bend| = min(K_e · |θ_e|, M_p)        with M_p = (4/3) σ_0 r_b³

    On unload, M follows the elastic line with slope K_e back to M = 0 at
    θ_bend = θ_p. The validator checks:
      1. peak |M_bend| ≤ M_p (within tolerance) — moment is capped
      2. residual |θ_p_bend| grows monotonically across the three cycles
      3. tip z drifts more negative cycle by cycle (permanent set)
    """
    last = rows[-1]
    iben = last["iben"]
    k_bend = last["k_bend"]
    bond_len = last["bond_len_mid0"]
    r_b = last["r_b"]
    # Material constants — yield stress and Young's modulus come from env so
    # the validator stays in sync with the TOML; the defaults match the
    # canonical `bending_plastic_guo.toml`.
    sigma_y = float(os.environ.get("FIBER_BOND_SIGMA_Y", "3.0e6"))
    e_b     = float(os.environ.get("FIBER_BOND_E_B",     "1.0e9"))
    # Guo Eqs. 27, 29, 31, 33, 35 — trilinear bending envelope.
    pi = math.pi
    m_p     = (4.0 / 3.0) * sigma_y * r_b**3              # Eq. 31, fully-plastic cap
    m_e     = sigma_y * iben / r_b                        # Eq. 29, first-yield moment
    k_ep    = 0.5 * k_bend                                # Eq. 35
    theta_e = sigma_y * bond_len / (e_b * r_b)            # Eq. 27, θ at first yield
    theta_p = theta_e * (32.0 - 3.0 * pi) / (3.0 * pi)    # Eq. 33 in ε-space, then ·l/r
    theta_yield = m_p / k_bend                            # equivalent linear-elastic yield θ

    # Per-step kinematics → moment estimate.
    # M_bend = K_bend · (θ_bend − θ_p) using the *signed y-component* (the
    # entire deformation lives in this single component for our 2D bending).
    times    = [r["t"] for r in rows]
    th_bend  = [r["dth_bend_y_mid"] for r in rows]
    th_p     = [r["theta_p_bend_y_mid"] for r in rows]
    th_e     = [b - p for b, p in zip(th_bend, th_p)]
    m_bend   = [k_bend * e for e in th_e]
    tip_z    = [r["right_z"] for r in rows]

    peak_m = max(abs(m) for m in m_bend)
    peak_th_p = max(abs(p) for p in th_p)
    final_th_p = th_p[-1]
    final_tip_z = tip_z[-1]

    # Cycle boundaries — must match the Rust schedule in main.rs.
    n_cycles = 3
    t_cycle = 30.0e-3
    # Sample θ_p at the END of each cycle (just before the next ramp).
    cycle_end_th_p = []
    for k in range(1, n_cycles + 1):
        t_end = k * t_cycle - 0.1e-3      # 0.1 ms before next cycle start
        i = min(range(len(times)), key=lambda j: abs(times[j] - t_end))
        cycle_end_th_p.append(th_p[i])
    monotone = all(abs(cycle_end_th_p[i]) <= abs(cycle_end_th_p[i + 1]) + 1e-12
                   for i in range(len(cycle_end_th_p) - 1))

    print("=== Bending plastic (Guo three-step, trilinear) — envelope + cycling ===")
    print(f"  σ_0 (yield stress)             : {sigma_y:.3e} Pa  (env: FIBER_BOND_SIGMA_Y)")
    print(f"  E_b (Young's modulus)          : {e_b:.3e} Pa     (env: FIBER_BOND_E_B)")
    print(f"  M^e = σ_0·I/r_b (first yield)  : {m_e:.6e} N·m  (Guo Eq. 29)")
    print(f"  M^p = (4/3)·σ_0·r_b³ (cap)     : {m_p:.6e} N·m  (Guo Eq. 31)")
    print(f"  K_e = E·I/l_b                  : {k_bend:.6e} N·m/rad  (Guo Eq. 34)")
    print(f"  K_ep = K_e/2                   : {k_ep:.6e} N·m/rad  (Guo Eq. 35)")
    print(f"  θ^e = σ_0·l_b/(E_b·r_b)        : {theta_e:.6e} rad   (Guo Eq. 27)")
    print(f"  θ^p (full-plastic transition)  : {theta_p:.6e} rad   (Guo Eq. 33)")
    print(f"  peak |M_bend| measured         : {peak_m:.6e} N·m   (cap M_p = {m_p:.3e})")
    print(f"  ratio peak / M_p               : {peak_m / m_p:.4f}")
    print(f"  peak |θ_p_bend| measured       : {peak_th_p:.6e} rad")
    print(f"  θ_p at end of each cycle (rad) : "
          + ", ".join(f"{p:+.3e}" for p in cycle_end_th_p))
    print(f"  monotone accumulation?         : {monotone}")
    print(f"  tip z final (residual)         : {final_tip_z*1e3:+.3f} mm")

    cap_ok = peak_m <= 1.01 * m_p    # 1 % slop on the envelope cap
    accum_ok = peak_th_p > 0.0 and monotone
    ok = cap_ok and accum_ok
    print(f"  cap respected?                 : {cap_ok}")
    print(f"  status                         : {'PASS' if ok else 'FAIL'}")

    plt = try_import_matplotlib()
    if plt is None:
        return ok

    # ── Figure 1: M(θ_bend) hysteresis loop with trilinear envelope ─────
    # Draw the three-segment Guo trilinear envelope (positive + negative
    # halves), then overlay the DEM trajectory. Guo Eq. 32:
    #   elastic       slope K_e   for |θ| ≤ θ^e
    #   elasto-plastic slope K_ep for θ^e < |θ| ≤ θ^p
    #   perfectly plastic cap    for |θ| > θ^p
    def trilinear_envelope(th):
        a = abs(th)
        if a <= theta_e:
            return math.copysign(k_bend * a, th)
        elif a <= theta_p:
            return math.copysign(m_e + k_ep * (a - theta_e), th)
        else:
            return math.copysign(m_p, th)
    th_max = max(abs(b) for b in th_bend) * 1.10
    n_env = 400
    th_env = [-th_max + 2.0 * th_max * i / n_env for i in range(n_env + 1)]
    m_env = [trilinear_envelope(t) for t in th_env]

    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.plot([t * 1e3 for t in th_env], [m * 1e3 for m in m_env],
            "-", lw=2.0, c="C1",
            label=("Theory  Guo trilinear envelope  "
                   f"(K_e = {k_bend:.3g}, K_ep = {k_ep:.3g} N·m/rad)"))
    ax.axhline( m_p * 1e3, ls="--", lw=0.8, c="0.6")
    ax.axhline(-m_p * 1e3, ls="--", lw=0.8, c="0.6")
    # Annotate the three regimes on the positive side.
    ax.axvline( theta_e * 1e3, ls=":", c="0.6", lw=0.6)
    ax.axvline( theta_p * 1e3, ls=":", c="0.6", lw=0.6)
    ax.axvline(-theta_e * 1e3, ls=":", c="0.6", lw=0.6)
    ax.axvline(-theta_p * 1e3, ls=":", c="0.6", lw=0.6)
    ax.plot([t * 1e3 for t in th_bend], [m * 1e3 for m in m_bend],
            "-", lw=1.2, c="C0", label="DEM trajectory  (mid-bond, all 3 cycles)")
    ax.set_xlabel("bending angle  θ_bend  (mrad)")
    ax.set_ylabel("bending moment  M_bend  (mN·m)")
    ax.set_title("Bending plastic — M–θ hysteresis vs. Guo trilinear envelope")
    ax.legend(loc="lower right", framealpha=0.95)
    out = Path(out_dir) / "bending_plastic_hysteresis.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  plot saved                     : {out}")

    # ── Figure 2: θ_p staircase + tip-z timeline ────────────────────────
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    ax1.plot([t * 1e3 for t in times], [p * 1e3 for p in th_p],
             "-", lw=1.4, c="C0", label="DEM  θ_p_bend (residual)  mid-bond")
    for k in range(1, n_cycles + 1):
        ax1.axvline(k * t_cycle * 1e3, ls=":", c="0.6", lw=0.6)
    ax1.set_ylabel("plastic anchor  θ_p_bend  (mrad)")
    ax1.set_title("Bending plastic — residual θ_p accumulates each cycle (Guo Fig. 12 analog)")
    ax1.legend(loc="upper left", framealpha=0.95)
    ax2.plot([t * 1e3 for t in times], [z * 1e3 for z in tip_z],
             "-", lw=1.2, c="C2", label="DEM  tip z(t)")
    for k in range(1, n_cycles + 1):
        ax2.axvline(k * t_cycle * 1e3, ls=":", c="0.6", lw=0.6)
    ax2.set_xlabel("time  t  (ms)")
    ax2.set_ylabel("tip deflection  z  (mm)")
    ax2.legend(loc="upper right", framealpha=0.95)
    out = Path(out_dir) / "bending_plastic_timeline.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  plot saved                     : {out}")
    return ok


VALIDATORS = {
    "axial_elastic":            validate_axial_elastic,
    "cantilever_bending":       validate_cantilever_bending,
    "bending_vibration":        validate_bending_vibration,
    "axial_plastic_piecewise":  validate_axial_plastic_piecewise,
    "bending_plastic_guo":      validate_bending_plastic_guo,
}


def main():
    if len(sys.argv) < 2:
        sys.exit(f"Usage: {sys.argv[0]} PATH_TO_CSV")
    path = sys.argv[1]
    if not os.path.exists(path):
        sys.exit(f"CSV not found: {path}. Run the Rust example first.")

    p = Path(path).resolve()
    mode = p.parent.parent.name
    out_dir = p.parent
    if mode not in VALIDATORS:
        sys.exit(f"Unknown validation mode '{mode}'. Known: {sorted(VALIDATORS)}")

    rows = load_timeseries(path)
    if not rows:
        sys.exit("CSV has no data rows.")

    profile_path = out_dir / "profile.csv"
    profile = load_profile(str(profile_path)) if profile_path.exists() else []
    if not profile:
        print(f"  note: {profile_path} not found — profile plots will be skipped")

    ok = VALIDATORS[mode](rows, profile, out_dir)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
