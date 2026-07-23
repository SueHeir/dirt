#!/usr/bin/env python3
"""
Validate one `fiber_bond_breakage` test run against analytical predictions
and produce side-by-side overlay plots.

The harness reuses the `fiber_bond` binary (one example registration; no
duplicate Rust code) and the same `fiber_bond.csv` / `profile.csv` /
`bond_thresholds.csv` schema. This script is the breakage-specific
validator counterpart to `examples/fiber_bond/validate.py`.

Mode dispatch picks the validator by the parent directory of the CSV —
e.g. `.../axial_stress_constant/data/fiber_bond.csv` runs
`validate_axial_stress_constant`.

Usage
-----
    python3 examples/fiber_bond_breakage/validate.py PATH_TO_CSV
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


def load_bond_thresholds(path):
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({
                "tag_a": int(float(row["tag_a"])),
                "tag_b": int(float(row["tag_b"])),
                "r0":    float(row["r0"]),
                "thr0":  float(row["thr0"]),
                "thr1":  float(row["thr1"]),
                "thr2":  float(row["thr2"]),
                "thr3":  float(row["thr3"]),
            })
    return rows


def first_break_event(rows):
    """Time and global strain at which `bond_count` first drops below its
    initial value. Returns `None` if no break happened during the run."""
    if not rows:
        return None
    initial_n = int(rows[0]["bond_count"])
    if initial_n == 0:
        return None
    length0 = rows[0]["length0"]
    for r in rows:
        if int(r["bond_count"]) < initial_n:
            eps = (r["length_global"] - length0) / length0
            return r["t"], eps
    return None


def try_import_matplotlib():
    try:
        import matplotlib.pyplot as plt  # type: ignore
        return plt
    except ImportError:
        print("  (matplotlib not available — skipping plots)")
        return None


def _common_axial_setup(rows):
    last = rows[-1]
    return {
        "length0":  last["length0"],
        "bond_len": last["bond_len_mid0"],
        "area":     last["area"],
        "k_n":      last["k_n"],
        "E":        last["k_n"] * last["bond_len_mid0"] / last["area"],
        # pull velocity in the canonical configs (matches `[[move_linear]] vx`)
        "vx":       float(os.environ.get("FIBER_BOND_VX", "0.1")),
    }


def _piecewise_envelope(eps_mag, e_mod, breakpoints, multipliers):
    sigma = 0.0
    eps_prev = 0.0
    slope = e_mod
    for bp, mult in zip(breakpoints, multipliers):
        if eps_mag <= bp:
            return sigma + slope * (eps_mag - eps_prev)
        sigma += slope * (bp - eps_prev)
        eps_prev = bp
        slope = e_mod * mult
    return sigma + slope * (eps_mag - eps_prev)


# ── Axial-stress + Constant threshold ────────────────────────────────────────

def validate_axial_stress_constant(rows, out_dir):
    """Predicted: a bond breaks at `ε = σ_max / E`, time `t = ε · L₀ / v`."""
    s = _common_axial_setup(rows)
    sigma_max = float(os.environ.get("FIBER_BOND_SIGMA_MAX", "5.0e6"))
    eps_pred = sigma_max / s["E"]
    t_pred = eps_pred * s["length0"] / s["vx"]
    ev = first_break_event(rows)
    print("=== Axial stress (constant threshold) — break prediction ===")
    print(f"  E (config)               : {s['E']:.4e} Pa")
    print(f"  σ_max (env)              : {sigma_max:.4e} Pa")
    print(f"  predicted ε_break        : {eps_pred:.6f}")
    print(f"  predicted t_break        : {t_pred*1e3:.4f} ms")
    if ev is None:
        print("  measured break event     : NONE — no bond broke during the run")
        print("  status                   : FAIL")
        return False
    t_meas, eps_meas = ev
    err_eps = abs(eps_meas - eps_pred) / eps_pred
    err_t   = abs(t_meas - t_pred) / t_pred
    print(f"  measured t / ε at break  : {t_meas*1e3:.4f} ms / ε = {eps_meas:.6f}")
    print(f"  relative error (ε)       : {err_eps:.3%}")
    print(f"  relative error (t)       : {err_t:.3%}")
    ok = err_eps < 0.05
    print(f"  status                   : {'PASS' if ok else 'FAIL'}  (tolerance 5%)")
    _plot_axial_break_event(rows, out_dir, eps_pred, t_pred,
                            title="Axial stress + Constant threshold",
                            mode_label=f"σ_max = {sigma_max:.2g} Pa",
                            filename="axial_stress_constant.png")
    return ok


# ── Axial-strain + Constant threshold ────────────────────────────────────────

def validate_axial_strain_constant(rows, out_dir):
    """Predicted: a bond breaks when ε_axial reaches ε_max."""
    s = _common_axial_setup(rows)
    eps_max = float(os.environ.get("FIBER_BOND_EPS_MAX", "0.005"))
    t_pred = eps_max * s["length0"] / s["vx"]
    ev = first_break_event(rows)
    print("=== Axial strain (constant threshold) — break prediction ===")
    print(f"  ε_max (env)              : {eps_max:.6f}")
    print(f"  predicted t_break        : {t_pred*1e3:.4f} ms")
    if ev is None:
        print("  measured break event     : NONE")
        print("  status                   : FAIL")
        return False
    t_meas, eps_meas = ev
    err_eps = abs(eps_meas - eps_max) / eps_max
    err_t = abs(t_meas - t_pred) / t_pred
    print(f"  measured t / ε at break  : {t_meas*1e3:.4f} ms / ε = {eps_meas:.6f}")
    print(f"  relative error (ε)       : {err_eps:.3%}")
    print(f"  relative error (t)       : {err_t:.3%}")
    ok = err_eps < 0.05
    print(f"  status                   : {'PASS' if ok else 'FAIL'}  (tolerance 5%)")
    _plot_axial_break_event(rows, out_dir, eps_max, t_pred,
                            title="Axial strain + Constant threshold",
                            mode_label=f"ε_max = {eps_max:.4f}",
                            filename="axial_strain_constant.png")
    return ok


# ── Axial-stress + Weibull thresholds ────────────────────────────────────────

def validate_axial_stress_weibull(rows, out_dir):
    """Predicted: weakest bond (lowest `thr0`) breaks first at
    `ε = thr0_min / E`. The validator reads per-bond thresholds from
    `bond_thresholds.csv` (written by the recorder at setup) — no
    re-implementation of the SplitMix64 / SmallRng / Weibull sampler in
    Python.
    """
    s = _common_axial_setup(rows)
    out = Path(out_dir)
    thr_path = out / "bond_thresholds.csv"
    if not thr_path.exists():
        print(f"  bond_thresholds.csv not found at {thr_path} — cannot validate")
        return False
    thrs = load_bond_thresholds(thr_path)
    if not thrs:
        print("  bond_thresholds.csv has no rows")
        return False
    tensile_vals = [t["thr0"] for t in thrs]
    sigma_min = min(tensile_vals)
    sigma_max = max(tensile_vals)
    sigma_mean = sum(tensile_vals) / len(tensile_vals)
    weakest = min(thrs, key=lambda t: t["thr0"])
    eps_pred = sigma_min / s["E"]
    t_pred = eps_pred * s["length0"] / s["vx"]
    ev = first_break_event(rows)
    print("=== Axial stress + Weibull thresholds — weakest-bond prediction ===")
    print(f"  E (config)               : {s['E']:.4e} Pa")
    print(f"  N bonds                  : {len(thrs)}")
    print(f"  σ thresholds (min/mean/max): "
          f"{sigma_min:.4e} / {sigma_mean:.4e} / {sigma_max:.4e} Pa")
    print(f"  weakest bond (tag pair)  : ({weakest['tag_a']}, {weakest['tag_b']})")
    print(f"  predicted ε_break        : {eps_pred:.6f}")
    print(f"  predicted t_break        : {t_pred*1e3:.4f} ms")
    if ev is None:
        print("  measured break event     : NONE")
        print("  status                   : FAIL")
        return False
    t_meas, eps_meas = ev
    err_eps = abs(eps_meas - eps_pred) / eps_pred
    err_t = abs(t_meas - t_pred) / t_pred
    print(f"  measured t / ε at break  : {t_meas*1e3:.4f} ms / ε = {eps_meas:.6f}")
    print(f"  relative error (ε)       : {err_eps:.3%}")
    print(f"  relative error (t)       : {err_t:.3%}")
    ok = err_eps < 0.05
    print(f"  status                   : {'PASS' if ok else 'FAIL'}  (tolerance 5%)")
    _plot_weibull_break(rows, thrs, out_dir, eps_pred, t_pred, s["E"])
    return ok


# ── Axial plasticity + AxialStress threshold ────────────────────────────────

def validate_axial_plastic_stress_constant(rows, out_dir):
    """Coupled axial piecewise plasticity and stress breakage.

    The reference is the analytical stress-strain envelope from the configured
    axial plastic return-map. The break threshold is deliberately in the
    post-yield hardening branch, so the predicted break strain depends on the
    plastic envelope:

        sigma = 10 MPa + 0.5 E (eps - 0.010) = 14 MPa
        eps_break = 0.018

    This exercises the coupling because an elastic-only prediction would break
    at eps = 0.014 instead.
    """
    s = _common_axial_setup(rows)
    sigma_break = float(os.environ.get("FIBER_BOND_SIGMA_BREAK", "14.0e6"))
    breakpoints = [0.01, 0.02, 0.03]
    multipliers = [0.5, 0.1, 0.0]
    e_mod = s["E"]
    area = s["area"]
    bond_len = s["bond_len"]
    k_n = s["k_n"]
    length0 = s["length0"]
    vx = s["vx"]

    eps_y = breakpoints[0]
    sigma_y = e_mod * eps_y
    # The canonical threshold is between first yield and the second
    # breakpoint. Keep the generic scan below as a guard if env overrides it.
    eps_pred = None
    eps_prev = 0.0
    sigma_prev = 0.0
    slope = e_mod
    for bp, mult in zip(breakpoints + [0.10], multipliers + [multipliers[-1]]):
        sigma_next = sigma_prev + slope * (bp - eps_prev)
        if sigma_break <= sigma_next or slope == 0.0:
            if slope <= 0.0:
                eps_pred = eps_prev if abs(sigma_break - sigma_prev) < 1e-9 else float("inf")
            else:
                eps_pred = eps_prev + (sigma_break - sigma_prev) / slope
            break
        sigma_prev = sigma_next
        eps_prev = bp
        slope = e_mod * mult
    if eps_pred is None:
        eps_pred = float("inf")
    t_pred = eps_pred * length0 / vx if math.isfinite(eps_pred) else float("inf")

    eps = [(r["length_global"] - length0) / length0 for r in rows]
    eps_mid = [r["strain_axial_mid"] for r in rows]
    # The recorder exposes the plastic anchor, so the measured load is the
    # actual coupled normal force after the return-map.
    stress = [
        k_n * (r["delta_mid"] - r["eps_p_axial_mid"] * bond_len) / area
        for r in rows
    ]
    eps_p = [r["eps_p_axial_mid"] for r in rows]
    ev = first_break_event(rows)

    # Envelope agreement up to the first break event. Use the mid-bond strain
    # for the force envelope itself; global strain is only for the first-break
    # event, matching the older axial breakage gates. Skip the tiny initial
    # transient by checking representative strains across elastic and plastic
    # branches rather than taking a max over near-zero stress.
    end = len(rows)
    if ev is not None:
        for i, r in enumerate(rows):
            if int(r["bond_count"]) < int(rows[0]["bond_count"]):
                end = max(1, i)
                break
    max_err = 0.0
    sample_eps = [0.005, 0.010, 0.014, 0.0175]
    for e_target in sample_eps:
        if e_target >= eps_pred:
            continue
        idx = min(range(end), key=lambda i: abs(eps_mid[i] - e_target))
        sig_i = stress[idx]
        sig_ref = _piecewise_envelope(e_target, e_mod, breakpoints, multipliers)
        if sig_ref > 1e-9:
            max_err = max(max_err, abs(sig_i - sig_ref) / sig_ref)

    print("=== Axial plastic + AxialStress breakage — coupled validation ===")
    print(f"  E (config)                    : {e_mod:.4e} Pa")
    print(f"  plastic breakpoints           : {breakpoints}")
    print(f"  plastic slope multipliers     : {multipliers}")
    print(f"  first-yield stress            : {sigma_y:.4e} Pa at eps = {eps_y:.4f}")
    print(f"  tensile break threshold       : {sigma_break:.4e} Pa")
    print(f"  predicted eps_break           : {eps_pred:.6f}")
    print(f"  elastic-only eps would be     : {sigma_break / e_mod:.6f}")
    print(f"  predicted t_break             : {t_pred*1e3:.4f} ms")
    print(f"  max envelope error pre-break  : {max_err:.3%}")
    if ev is None:
        print("  measured break event          : NONE")
        print("  status                        : FAIL")
        return False
    t_meas, eps_meas = ev
    err_eps = abs(eps_meas - eps_pred) / eps_pred
    err_t = abs(t_meas - t_pred) / t_pred
    had_plasticity = max(abs(p) for p in eps_p[:end]) > 1e-5
    print(f"  measured t / eps at break     : {t_meas*1e3:.4f} ms / eps = {eps_meas:.6f}")
    print(f"  relative error (eps)          : {err_eps:.3%}")
    print(f"  relative error (t)            : {err_t:.3%}")
    print(f"  plastic strain active?        : {had_plasticity}")
    ok = (max_err < 0.02) and (err_eps < 0.05) and had_plasticity
    print(f"  status                        : {'PASS' if ok else 'FAIL'}  "
          "(envelope 2%, break 5%)")
    _plot_axial_plastic_breakage(
        rows, out_dir, eps, stress, eps_p, breakpoints, multipliers,
        e_mod, sigma_break, eps_pred, t_pred, ok,
    )
    return ok


# ── Plot helpers ─────────────────────────────────────────────────────────────

def _plot_axial_break_event(rows, out_dir, eps_pred, t_pred, title, mode_label, filename):
    plt = try_import_matplotlib()
    if plt is None:
        return
    initial_n = int(rows[0]["bond_count"])
    length0 = rows[-1]["length0"]
    t_ms   = [r["t"] * 1e3 for r in rows]
    eps    = [(r["length_global"] - length0) / length0 for r in rows]
    nbonds = [r["bond_count"] for r in rows]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    ax1.plot(t_ms, eps, "-", lw=1.2, c="C0", label="DEM global strain ε(t)")
    ax1.axhline(eps_pred, ls="--", lw=1.6, c="C1",
                label=f"Theory  ε_break = {eps_pred:.4f}  ({mode_label})")
    ax1.axvline(t_pred * 1e3, ls=":", lw=1.0, c="C1",
                label=f"Theory  t_break = {t_pred*1e3:.2f} ms")
    ax1.set_ylabel("axial strain  ε")
    ax1.set_title(title)
    ax1.legend(loc="upper left", framealpha=0.95)
    ax2.plot(t_ms, nbonds, "-", lw=1.4, c="C2", label="DEM bond_count(t)")
    ax2.axvline(t_pred * 1e3, ls=":", lw=1.0, c="C1")
    ax2.set_xlabel("time  t  (ms)")
    ax2.set_ylabel("bonds remaining")
    ax2.set_ylim(-0.5, initial_n + 0.5)
    ax2.legend(loc="lower left", framealpha=0.95)
    fig.tight_layout()
    out = Path(out_dir) / filename
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  plot saved               : {out}")


def _plot_weibull_break(rows, thrs, out_dir, eps_pred, t_pred, e_mod):
    plt = try_import_matplotlib()
    if plt is None:
        return
    initial_n = int(rows[0]["bond_count"])
    length0 = rows[-1]["length0"]
    t_ms   = [r["t"] * 1e3 for r in rows]
    eps    = [(r["length_global"] - length0) / length0 for r in rows]
    nbonds = [r["bond_count"] for r in rows]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 6.5))

    # Panel 1: per-bond sampled thresholds (Weibull) → equivalent strain ε = σ/E.
    eps_thrs = sorted(t["thr0"] / e_mod for t in thrs)
    ranks = list(range(1, len(eps_thrs) + 1))
    ax1.plot(eps_thrs, ranks, "o", ms=8, mfc="white", mec="C0", mew=1.4,
             label="DEM  per-bond ε_break (sorted)  — Weibull-sampled")
    ax1.axvline(eps_pred, ls="--", lw=1.6, c="C1",
                label=f"Theory  weakest = min(ε_break) = {eps_pred:.4f}")
    ax1.set_xlabel("per-bond breaking strain  ε  = σ_break / E")
    ax1.set_ylabel("sorted rank")
    ax1.set_title("Per-bond Weibull-sampled tensile thresholds")
    ax1.legend(loc="lower right", framealpha=0.95)

    # Panel 2: trajectory + bond_count with prediction lines.
    ax2.plot(t_ms, eps, "-", lw=1.2, c="C0", label="DEM global strain ε(t)")
    ax2.axhline(eps_pred, ls="--", lw=1.6, c="C1",
                label=f"Predicted ε_break = {eps_pred:.4f}")
    ax2.axvline(t_pred * 1e3, ls=":", lw=1.0, c="C1",
                label=f"Predicted t_break = {t_pred*1e3:.2f} ms")
    ax2t = ax2.twinx()
    ax2t.plot(t_ms, nbonds, "-", lw=1.4, c="C2", label="DEM bond_count")
    ax2t.set_ylim(-0.5, initial_n + 0.5)
    ax2t.set_ylabel("bonds remaining", color="C2")
    ax2.set_xlabel("time  t  (ms)")
    ax2.set_ylabel("axial strain  ε", color="C0")
    ax2.set_title("Strain & bond-count trajectory")
    ax2.legend(loc="upper left", framealpha=0.95)
    ax2t.legend(loc="center right", framealpha=0.95)

    fig.tight_layout()
    out = Path(out_dir) / "axial_stress_weibull.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  plot saved               : {out}")


def _plot_axial_plastic_breakage(rows, out_dir, eps, stress, eps_p, breakpoints,
                                 multipliers, e_mod, sigma_break, eps_pred,
                                 t_pred, ok):
    plt = try_import_matplotlib()
    if plt is None:
        return
    initial_n = int(rows[0]["bond_count"])
    nbonds = [r["bond_count"] for r in rows]
    eps_max = max(max(eps), eps_pred) * 1.08
    eps_th = [eps_max * i / 500.0 for i in range(501)]
    area = rows[-1]["area"]
    force = [s * area for s in stress]
    f_th = [_piecewise_envelope(e, e_mod, breakpoints, multipliers) * area for e in eps_th]
    f_break = sigma_break * area

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 6.5), sharex=True)
    ax1.plot([e * 100 for e in eps_th], f_th,
             "-", lw=1.9, c="C1",
             label="Analytical elastic-plastic envelope")
    ax1.plot([e * 100 for e in eps], force,
             "o", ms=3, mfc="white", mec="C0", mew=0.8,
             label="DEM mid-bond axial force after plastic return")
    ax1.axhline(f_break, ls="--", lw=1.2, c="0.25",
                label=f"Break threshold = {f_break:.2f} N")
    ax1.axvline(eps_pred * 100, ls=":", lw=1.4, c="C3",
                label=f"Predicted break strain = {eps_pred*100:.2f}%")
    for bp in breakpoints:
        ax1.axvline(bp * 100, ls=":", c="0.7", lw=0.7)
    ax1.set_ylabel("axial force  F_n  (N)")
    ax1.set_title("Coupled axial plasticity + breakage force-strain validation")
    ax1.legend(loc="lower right", framealpha=0.95)

    ax2.plot([e * 100 for e in eps], [p * 100 for p in eps_p],
             "-", lw=1.3, c="C2", label="DEM plastic strain anchor")
    ax2b = ax2.twinx()
    ax2b.step([e * 100 for e in eps], nbonds, where="post",
              lw=1.4, c="C4", label="DEM bond_count")
    ax2.axvline(eps_pred * 100, ls=":", lw=1.4, c="C3")
    ax2.set_xlabel("global axial strain  eps  (%)")
    ax2.set_ylabel("plastic axial strain  eps_p  (%)", color="C2")
    ax2b.set_ylabel("bonds remaining", color="C4")
    ax2b.set_ylim(-0.5, initial_n + 0.5)
    ax2.legend(loc="upper left", framealpha=0.95)
    ax2b.legend(loc="center right", framealpha=0.95)
    fig.tight_layout()

    out = Path(out_dir) / "axial_plastic_stress_constant.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  plot saved                  : {out}")

    # Also write the committed summary plot beside the other validation figures.
    repo_plot_dir = Path(__file__).resolve().parent / "plots"
    repo_plot_dir.mkdir(parents=True, exist_ok=True)
    out2 = repo_plot_dir / "plastic_breakage_coupled_validation.png"
    fig.savefig(out2, dpi=150, bbox_inches="tight")
    print(f"  plot saved                  : {out2}")


# ── Cantilever-bend break event extraction ─────────────────────────────────

def _cantilever_break_event(rows):
    """First-break event for a cantilever-bend scenario: time, tip-z, and
    bending angle of the mid-bond at the moment `bond_count` drops below
    its initial value."""
    initial_n = int(rows[0]["bond_count"])
    for r in rows:
        if int(r["bond_count"]) < initial_n:
            return {
                "t":         r["t"],
                "tip_z":     r["right_z"],
                "dth_bend":  r["dth_bend_y_mid"],
                "length_global": r["length_global"],
            }
    return None


def _cantilever_setup(rows):
    last = rows[-1]
    e_b   = last["k_n"] * last["bond_len_mid0"] / last["area"]
    return {
        "iben":     last["iben"],
        "area":     last["area"],
        "k_bend":   last["k_bend"],
        "bond_len": last["bond_len_mid0"],
        "r_b":      last["r_b"],
        "length0":  last["length0"],   # = L_c (anchor-to-tip)
        "E":        e_b,
        "vz":       float(os.environ.get("FIBER_BOND_VZ", "-0.5")),
    }


# ── CombinedStress — cantilever bend ───────────────────────────────────────

def validate_combined_stress(rows, out_dir):
    """Anchor bond breaks at `σ_combined = r_b·M_anchor/I = σ_max` (axial part
    is ~0 in pure transverse bending). Predicted tip displacement at break:
    `y_break = σ_max·L_c² / (3·E·r_b)`."""
    s = _cantilever_setup(rows)
    sigma_max = float(os.environ.get("FIBER_BOND_SIGMA_MAX", "5.0e6"))
    y_break = sigma_max * s["length0"]**2 / (3.0 * s["E"] * s["r_b"])
    t_break = abs(y_break / s["vz"])
    ev = _cantilever_break_event(rows)
    print("=== CombinedStress (Guo Eq. 16) — cantilever break prediction ===")
    print(f"  σ_max                       : {sigma_max:.4e} Pa")
    print(f"  L_c, r_b                    : {s['length0']*1e3:.1f} mm, {s['r_b']*1e3:.2f} mm")
    print(f"  E                           : {s['E']:.3e} Pa")
    print(f"  predicted y_break           : {-y_break*1e3:+.4f} mm   (theory −y)")
    print(f"  predicted t_break           : {t_break*1e3:.4f} ms")
    if ev is None:
        print("  measured break event        : NONE")
        print("  status                      : FAIL")
        return False
    err = abs(abs(ev["tip_z"]) - y_break) / y_break
    print(f"  measured tip_z at break     : {ev['tip_z']*1e3:+.4f} mm   at t = {ev['t']*1e3:.4f} ms")
    print(f"  relative error (|tip_z|)    : {err:.3%}")
    # Loose tolerance — small-def quasi-static EB is a rough approximation
    # for the dynamic + discretised chain near the pinned end.
    ok = err < 0.35
    print(f"  status                      : {'PASS' if ok else 'FAIL'}  (tolerance 35%)")
    _plot_cantilever_break(rows, out_dir, -y_break, t_break,
                           title="CombinedStress — cantilever bend",
                           mode_label=f"σ_max = {sigma_max:.2g} Pa",
                           filename="combined_stress.png")
    return ok


# ── CombinedStrain — cantilever bend ───────────────────────────────────────

def validate_combined_strain(rows, out_dir):
    """Anchor bond breaks at `ε_combined = r_b·κ_anchor = ε_max`. Predicted
    tip displacement: `y_break = ε_max·L_c² / (3·r_b)` (identical to the
    stress test when `ε_max = σ_max/E`)."""
    s = _cantilever_setup(rows)
    eps_max = float(os.environ.get("FIBER_BOND_EPS_MAX", "0.005"))
    y_break = eps_max * s["length0"]**2 / (3.0 * s["r_b"])
    t_break = abs(y_break / s["vz"])
    ev = _cantilever_break_event(rows)
    print("=== CombinedStrain (migration doc Eq. 1.7-1) — cantilever break ===")
    print(f"  ε_max                       : {eps_max:.4e}")
    print(f"  L_c, r_b                    : {s['length0']*1e3:.1f} mm, {s['r_b']*1e3:.2f} mm")
    print(f"  predicted y_break           : {-y_break*1e3:+.4f} mm")
    print(f"  predicted t_break           : {t_break*1e3:.4f} ms")
    if ev is None:
        print("  measured break event        : NONE")
        print("  status                      : FAIL")
        return False
    err = abs(abs(ev["tip_z"]) - y_break) / y_break
    print(f"  measured tip_z at break     : {ev['tip_z']*1e3:+.4f} mm   at t = {ev['t']*1e3:.4f} ms")
    print(f"  relative error (|tip_z|)    : {err:.3%}")
    ok = err < 0.35
    print(f"  status                      : {'PASS' if ok else 'FAIL'}  (tolerance 35%)")
    _plot_cantilever_break(rows, out_dir, -y_break, t_break,
                           title="CombinedStrain — cantilever bend",
                           mode_label=f"ε_max = {eps_max:.4f}",
                           filename="combined_strain.png")
    return ok


# ── InteractionLinearStress — cantilever bend ──────────────────────────────

def validate_interaction_linear_stress(rows, out_dir):
    """Bending-only InteractionLinear envelope at the anchor bond.

    Only the bending channel is active in the canonical config (the shear
    channel is omitted because at v_z = -0.5 m/s the dynamic damping term
    `γ_t · v_t` at the tip's bond produces a large transient shear stress —
    not the quasi-static load we want to validate against). So the envelope
    `Σ |X_i|/X_i,c ≥ 1` reduces to a single-channel condition equivalent to
    a pure CombinedStress with `σ_max = σ_bend,c`. Predicted tip
    displacement at break:

        y_break = σ_bend,c · L_c² / (3·E·r_b)

    The same prediction as `combined_stress.toml`; what we're actually
    validating here is that the multi-channel InteractionLinear code path
    (with three of four channels disabled via `None` thresholds) reaches
    the same break event."""
    s = _cantilever_setup(rows)
    sigma_bend = float(os.environ.get("FIBER_BOND_SIGMA_BEND", "5.0e6"))
    y_break = sigma_bend * s["length0"]**2 / (3.0 * s["E"] * s["r_b"])
    t_break = abs(y_break / s["vz"])
    ev = _cantilever_break_event(rows)
    print("=== InteractionLinearStress (bending-only channel) — cantilever break ===")
    print(f"  σ_bend,c                    : {sigma_bend:.4e} Pa")
    print(f"  L_c, r_b                    : {s['length0']*1e3:.1f} mm, {s['r_b']*1e3:.2f} mm")
    print(f"  E                           : {s['E']:.3e} Pa")
    print(f"  predicted y_break           : {-y_break*1e3:+.4f} mm")
    print(f"  predicted t_break           : {t_break*1e3:.4f} ms")
    if ev is None:
        print("  measured break event        : NONE")
        print("  status                      : FAIL")
        return False
    err = abs(abs(ev["tip_z"]) - y_break) / y_break
    print(f"  measured tip_z at break     : {ev['tip_z']*1e3:+.4f} mm   at t = {ev['t']*1e3:.4f} ms")
    print(f"  relative error (|tip_z|)    : {err:.3%}")
    # Loose tolerance — small-def quasi-static EB is a rough approximation
    # for the dynamic + discretised chain near the pinned end.
    ok = err < 0.35
    print(f"  status                      : {'PASS' if ok else 'FAIL'}  (tolerance 35%)")
    _plot_cantilever_break(
        rows, out_dir, -y_break, t_break,
        title="InteractionLinearStress (bending-only) — cantilever bend",
        mode_label=f"σ_bend,c = {sigma_bend:.2g} Pa",
        filename="interaction_linear_stress.png",
    )
    return ok


def _plot_cantilever_break(rows, out_dir, y_break, t_break, title, mode_label, filename):
    plt = try_import_matplotlib()
    if plt is None:
        return
    initial_n = int(rows[0]["bond_count"])
    t_ms   = [r["t"] * 1e3 for r in rows]
    tip_mm = [r["right_z"] * 1e3 for r in rows]
    nbonds = [r["bond_count"] for r in rows]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    ax1.plot(t_ms, tip_mm, "-", lw=1.2, c="C0", label="DEM tip z(t)")
    ax1.axhline(y_break * 1e3, ls="--", lw=1.6, c="C1",
                label=f"Theory  y_break = {y_break*1e3:+.3f} mm")
    ax1.axvline(t_break * 1e3, ls=":", lw=1.0, c="C1",
                label=f"Theory  t_break = {t_break*1e3:.2f} ms")
    ax1.set_ylabel("tip deflection  z  (mm)")
    ax1.set_title(f"{title}\n{mode_label}")
    ax1.legend(loc="upper right", framealpha=0.95)
    ax2.plot(t_ms, nbonds, "-", lw=1.4, c="C2", label="DEM bond_count(t)")
    ax2.axvline(t_break * 1e3, ls=":", lw=1.0, c="C1")
    ax2.set_xlabel("time  t  (ms)")
    ax2.set_ylabel("bonds remaining")
    ax2.set_ylim(-0.5, initial_n + 0.5)
    ax2.legend(loc="lower left", framealpha=0.95)
    fig.tight_layout()
    out = Path(out_dir) / filename
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  plot saved                  : {out}")


# ── Dispatch ─────────────────────────────────────────────────────────────────

VALIDATORS = {
    "axial_stress_constant":      validate_axial_stress_constant,
    "axial_strain_constant":      validate_axial_strain_constant,
    "axial_stress_weibull":       validate_axial_stress_weibull,
    "axial_plastic_stress_constant": validate_axial_plastic_stress_constant,
    "combined_stress":            validate_combined_stress,
    "combined_strain":            validate_combined_strain,
    "interaction_linear_stress":  validate_interaction_linear_stress,
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
    ok = VALIDATORS[mode](rows, out_dir)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
