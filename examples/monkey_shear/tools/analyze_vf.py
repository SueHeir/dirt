#!/usr/bin/env python3
"""analyze_vf — reduce the monkey-barrel LEBC volume-fraction shear campaign.

Reads the per-run `lebc_shear_results.csv` files produced by the `monkey_shear`
binary (one per (type, Phi) cell, under <data>/<type>/phi_<val>/) and extracts
the steady-shear rheology for each cell:

  * shear strain  gamma  = cumulative  gdot * dt  over the shear stage
                  (settle/compress rows carry gdot = 0 and are skipped)
  * steady window = strain in [GAMMA_SS_LO, gamma_max]  (default gamma > 1.0,
                    or the last STEADY_FRAC of the shear if gamma_max < 1.0)
  * p, tau        = window-mean pressure and shear stress  [Pa]
  * mu = tau / p  (macroscopic friction)
  * T             = window-mean granular temperature       [m^2/s^2]
  * I  = gdot * D_eq / sqrt(p / rho_grain)                  (inertial number)

It classifies every cell honestly:
  COMPLETE  — reached gamma >= GAMMA_DONE (default 1.8, target was ~2)
  PARTIAL   — reached shear but aborted before GAMMA_DONE (still has a usable
              steady estimate if gamma_max > GAMMA_SS_LO)
  NO-SHEAR  — aborted during settle/compress; gdot never went positive
              (no rheology point; typically a compression-stage overlap abort)

No point is fabricated. NO-SHEAR cells yield no (mu, I) datum.

Usage:
  python3 analyze_vf.py [--data <dir>] [--csv out.csv] [--plot fig.png]
Default --data is the campaign output tree; override to point at a snapshot.
"""
from __future__ import annotations
import argparse
import csv
import math
import os
import sys

# ── Campaign constants (must match tools/gen_series.py) ──────────────────────
D_EQ = 0.1            # common equivalent-volume diameter [m]
RHO_GRAIN = 2500.0    # grain (glass) density [kg/m^3]
TYPES = ["sphere", "rigid", "bpm"]
PHIS = [0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55]

# ── Steady-state / completeness thresholds ───────────────────────────────────
GAMMA_SS_LO = 1.0     # start of the steady averaging window (strain units)
GAMMA_DONE = 1.8      # >= this counts as COMPLETE (shear target was ~2)
STEADY_FRAC = 0.3     # fallback: last 30% of shear if gamma_max < GAMMA_SS_LO

DEFAULT_DATA = os.path.join(
    os.path.expanduser("~"),
    "projects/worktrees/dirt-monkey-shear-lebc-campaign/examples/monkey_shear/data",
)


def load_rows(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    out = []
    for r in rows:
        try:
            out.append({k: float(v) for k, v in r.items()})
        except (TypeError, ValueError):
            continue  # skip a torn final line from a live run
    return out


def analyze_cell(path):
    """Return a dict of reduced quantities for one run CSV (or a status dict)."""
    rows = load_rows(path)
    if len(rows) < 2:
        return {"status": "NO-SHEAR", "reason": "empty/short csv", "n_rows": len(rows)}

    # shear rows = gdot > 0; accumulate strain gamma = sum gdot*dt across them.
    gamma = 0.0
    shear = []  # (gamma, row)
    prev_t = None
    prev_shear = False
    for r in rows:
        gd = r.get("gdot", 0.0)
        t = r.get("time", 0.0)
        if gd > 0.0:
            if prev_shear and prev_t is not None:
                gamma += gd * (t - prev_t)
            shear.append((gamma, r))
            prev_shear = True
        else:
            prev_shear = False
        prev_t = t

    if not shear:
        return {"status": "NO-SHEAR", "reason": "gdot never positive (aborted in settle/compress)",
                "n_rows": len(rows), "phi_reached": rows[-1].get("phi", float("nan"))}

    gamma_max = shear[-1][0]
    gdot = shear[-1][1].get("gdot", 1.0)

    # steady window
    lo = GAMMA_SS_LO if gamma_max >= GAMMA_SS_LO else gamma_max * (1.0 - STEADY_FRAC)
    win = [r for (g, r) in shear if g >= lo]
    if len(win) < 2:
        win = [r for (g, r) in shear[-max(2, len(shear)//5):]]

    def mean(key):
        vals = [w[key] for w in win if key in w and math.isfinite(w[key])]
        return sum(vals) / len(vals) if vals else float("nan")

    p = mean("p")
    tau = mean("tau")
    T = mean("T")
    phi = mean("phi")
    N1 = mean("N1")
    N2 = mean("N2")
    mu = tau / p if (p and math.isfinite(p) and p != 0) else float("nan")
    I = (gdot * D_EQ / math.sqrt(p / RHO_GRAIN)) if (p and p > 0) else float("nan")

    status = "COMPLETE" if gamma_max >= GAMMA_DONE else "PARTIAL"
    return {
        "status": status, "gamma_max": gamma_max, "gdot": gdot,
        "p": p, "tau": tau, "mu": mu, "T": T, "phi": phi,
        "N1": N1, "N2": N2, "I": I,
        "n_shear_rows": len(shear), "n_win": len(win), "gamma_ss_lo": lo,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=DEFAULT_DATA)
    ap.add_argument("--csv", default=None, help="write reduced table here")
    ap.add_argument("--plot", default=None, help="write mu(I)/phi(I) figure here")
    args = ap.parse_args()

    table = []  # (type, phi_nominal, result)
    for t in TYPES:
        for phi in PHIS:
            path = os.path.join(args.data, t, f"phi_{phi}", "lebc_shear_results.csv")
            if not os.path.exists(path):
                table.append((t, phi, {"status": "MISSING"}))
                continue
            table.append((t, phi, analyze_cell(path)))

    # ── console report ───────────────────────────────────────────────────────
    hdr = f"{'type':7} {'phi_nom':7} {'status':9} {'gamma':6} {'phi':7} {'p[Pa]':11} {'tau[Pa]':11} {'mu':7} {'I':9} {'T':11}"
    print(hdr)
    print("-" * len(hdr))
    for t, phi, r in table:
        s = r["status"]
        if s in ("MISSING", "NO-SHEAR"):
            extra = r.get("reason", "")
            pr = r.get("phi_reached", float("nan"))
            print(f"{t:7} {phi:<7} {s:9} {'-':6} {pr:7.3f} {'-':>11} {'-':>11} {'-':7} {'-':9}   {extra}")
            continue
        print(f"{t:7} {phi:<7} {s:9} {r['gamma_max']:6.2f} {r['phi']:7.4f} "
              f"{r['p']:11.4e} {r['tau']:11.4e} {r['mu']:7.4f} {r['I']:9.4e} {r['T']:11.4e}")

    # ── optional csv dump ────────────────────────────────────────────────────
    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["type", "phi_nom", "status", "gamma_max", "phi_meas",
                        "p", "tau", "mu", "I", "T", "N1", "N2", "gdot"])
            for t, phi, r in table:
                if r["status"] in ("MISSING", "NO-SHEAR"):
                    w.writerow([t, phi, r["status"], "", "", "", "", "", "", "", "", "", ""])
                else:
                    w.writerow([t, phi, r["status"], f"{r['gamma_max']:.4f}",
                                f"{r['phi']:.5f}", f"{r['p']:.6e}", f"{r['tau']:.6e}",
                                f"{r['mu']:.5f}", f"{r['I']:.6e}", f"{r['T']:.6e}",
                                f"{r['N1']:.6e}", f"{r['N2']:.6e}", f"{r['gdot']:.4f}"])
        print(f"\nwrote {args.csv}")

    # ── optional plot ────────────────────────────────────────────────────────
    if args.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except Exception as e:
            print(f"[plot skipped: {e}]")
            return
        colors = {"sphere": "#1f77b4", "rigid": "#d62728", "bpm": "#2ca02c"}
        fig, axes = plt.subplots(1, 4, figsize=(19, 4.6))
        # highest nominal phi that produced a shear datum, per type -> ceiling band
        ceil = {}
        for t in TYPES:
            # collect points with a finite steady-ish datum; use NOMINAL phi as the
            # true solid fraction (recorder 'phi' = sub-sphere sum, ~1.31x for monkeys)
            pts = [(phin, r) for (tt, phin, r) in table
                   if tt == t and r["status"] in ("COMPLETE", "PARTIAL")
                   and math.isfinite(r.get("I", float("nan")))]
            noshear = [phin for (tt, phin, r) in table
                       if tt == t and r["status"] == "NO-SHEAR"]
            if pts:
                ceil[t] = max(p for p, _ in pts)
            if not pts:
                continue
            I = [r["I"] for _, r in pts]
            mu = [r["mu"] for _, r in pts]
            phin = [p for p, _ in pts]
            pres = [r["p"] for _, r in pts]
            part = [r["status"] == "PARTIAL" for _, r in pts]
            axes[0].plot(I, mu, "o-", color=colors[t], label=t, mfc="white", ms=6)
            axes[1].plot(I, phin, "o-", color=colors[t], label=t, mfc="white", ms=6)
            axes[2].plot(phin, pres, "o-", color=colors[t], label=t, mfc="white", ms=6)
            axes[3].plot(phin, mu, "o-", color=colors[t], label=t, mfc="white", ms=6)
            for (ii, mm, pn, pr, isp) in zip(I, mu, phin, pres, part):
                if isp:  # mark partial/aborted-early with an ×
                    axes[0].plot([ii], [mm], "x", color=colors[t], ms=9)
                    axes[1].plot([ii], [pn], "x", color=colors[t], ms=9)
                    axes[2].plot([pn], [pr], "x", color=colors[t], ms=9)
                    axes[3].plot([pn], [mm], "x", color=colors[t], ms=9)
            # mark the first NO-SHEAR (compaction/blow-up ceiling) for monkeys
            if noshear and t != "sphere":
                axes[2].axvline(min(noshear), color=colors[t], ls="--", alpha=0.5)
        axes[0].set_xscale("log"); axes[0].set_yscale("log")
        axes[0].set_xlabel("I"); axes[0].set_ylabel(r"$\mu=\tau/p$"); axes[0].set_title(r"$\mu(I)$")
        axes[1].set_xscale("log"); axes[1].set_xlabel("I"); axes[1].set_ylabel(r"$\phi$ (nominal, true solid fraction)"); axes[1].set_title(r"$\phi(I)$")
        axes[2].set_yscale("log"); axes[2].set_xlabel(r"$\phi$ (nominal)"); axes[2].set_ylabel("p [Pa]")
        axes[2].set_title(r"$p(\phi)$ — interlocking / jamming ceiling (-- = 1st no-shear $\phi$)")
        axes[3].set_yscale("log"); axes[3].set_xlabel(r"$\phi$ (nominal)"); axes[3].set_ylabel(r"$\mu=\tau/p$"); axes[3].set_title(r"$\mu(\phi)$")
        for ax in axes:
            ax.grid(True, alpha=0.3, which="both"); ax.legend(fontsize=8)
        fig.suptitle("Monkey-barrel LEBC VF campaign — rheology  (○ COMPLETE, × PARTIAL/aborted-early; "
                     "monkey $\\phi$≥0.3 never sheared)")
        fig.tight_layout()
        fig.savefig(args.plot, dpi=130)
        print(f"wrote {args.plot}")


if __name__ == "__main__":
    main()
