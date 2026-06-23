#!/usr/bin/env python3
"""
Isotropic-compression bulk modulus K(Φ) of glass beads.

A triperiodic box of glass beads (gravity off) is compressed quasi-statically and
ISOTROPICALLY — the [deform] vel driver pushes all three box faces inward at equal
speed — while the recorder streams pressure p = −trace(VirialStress)/(3·V) and
solid fraction Φ = Σ(π/6)d³/V each thermo interval. The pack passes through random
close packing (Φ_c ≈ 0.58–0.64), where enduring contacts percolate and the
pressure climbs steeply: the equation of state P(Φ).

The bulk compressibility closure the MUD SPH solver consumes is the modulus

    K(Φ) ≈ ΔP / (ΔΦ/Φ) = dP/d(lnΦ),

extracted here by fitting P(Φ) on the dense, jammed branch (P above a small
threshold) and differentiating. We fit an empirical EOS P(Φ) = A·(Φ − Φ_c)^α to
the jammed branch (the contact-mechanics expectation is α ≈ 1.5 for Hertzian
grains near jamming) and report K(Φ) = Φ · dP/dΦ at representative solid fractions.

Commands:
    python3 examples/SPH_glass_sphere_calibration/02_compressibility/sweep.py generate
    python3 examples/SPH_glass_sphere_calibration/02_compressibility/sweep.py start
    python3 examples/SPH_glass_sphere_calibration/02_compressibility/sweep.py graph
    python3 examples/SPH_glass_sphere_calibration/02_compressibility/sweep.py        # all three

Outputs:
    data/compressibility_results.csv     p, Φ time series          (gitignored)
    plots/eos_p_of_phi.png, bulk_modulus.png                       (tracked)

LAMMPS overlay: omitted. LAMMPS's granular package has no native isotropic
box-deform-while-measuring-virial rig comparable to DIRT's [deform] vel driver
(fix deform exists but the pressure/EOS extraction would not be an apples-to-apples
contact-model match), so the closure is validated against the contact-mechanics
EOS form rather than a cross-code overlay.
"""

import os
import sys
import csv
import math
import shutil
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", ".."))
EXAMPLE = "sphcal_compressibility"
SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
CONFIG = os.path.join(SCRIPT_DIR, "config.toml")
DATA = os.path.join(DATA_DIR, "compressibility_results.csv")

# Pressure threshold (Pa) above which the pack is treated as the jammed/dense
# branch the EOS describes; below it the response is the loose-pack transient.
P_JAM = 1.0e3
# Φ values at which to report the bulk modulus K = dP/d(lnΦ).
PHI_REPORT = [0.60, 0.62, 0.63]


# ── data ─────────────────────────────────────────────────────────────────────
def load():
    """Return (phi, p) sorted by phi, with the loose-pack settle transient dropped.

    The settle stage holds a fixed box (constant volume) while the high-pressure
    insertion overlap relaxes — those samples carry large p at low Φ and must NOT
    be mistaken for the jammed branch. We keep only the compression stage, detected
    as samples taken once the box volume has started shrinking below its initial
    value, so Φ increases monotonically thereafter."""
    rows = list(csv.DictReader(open(DATA)))
    vols = [float(r["volume"]) for r in rows]
    v0 = max(vols) if vols else 0.0
    pts = []
    for r, v in zip(rows, vols):
        if v < v0 * (1.0 - 1.0e-4):  # box has begun to shrink → compression stage
            pts.append((float(r["phi"]), float(r["p"])))
    pts.sort()
    phi = [a for a, _ in pts]
    p = [b for _, b in pts]
    return phi, p


def jammed_branch(phi, p):
    """Points on the dense branch: p ≥ P_JAM (enduring contacts percolated)."""
    return zip(*[(a, b) for a, b in zip(phi, p) if b >= P_JAM]) if any(
        b >= P_JAM for b in p) else ([], [])


# ── EOS fit:  P(Φ) = A·(Φ − Φ_c)^α  on the jammed branch ─────────────────────
def fit_eos(phi_j, p_j):
    """Least-squares fit of ln P = ln A + α ln(Φ − Φ_c), scanning Φ_c for best R²."""
    best = None
    phi_min = min(phi_j)
    for k in range(1, 60):
        phi_c = phi_min - 0.001 - 0.004 * k  # Φ_c just below the lowest jammed Φ
        xs, ys = [], []
        for a, b in zip(phi_j, p_j):
            d = a - phi_c
            if d <= 0 or b <= 0:
                continue
            xs.append(math.log(d))
            ys.append(math.log(b))
        if len(xs) < 3:
            continue
        n = len(xs)
        mx = sum(xs) / n
        my = sum(ys) / n
        sxx = sum((x - mx) ** 2 for x in xs)
        sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
        if sxx <= 0:
            continue
        alpha = sxy / sxx
        lnA = my - alpha * mx
        ss_res = sum((y - (lnA + alpha * x)) ** 2 for x, y in zip(xs, ys))
        ss_tot = sum((y - my) ** 2 for y in ys)
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        if best is None or r2 > best[0]:
            best = (r2, phi_c, alpha, math.exp(lnA))
    return best  # (r2, phi_c, alpha, A)


def eos_p(phi, phi_c, alpha, A):
    d = phi - phi_c
    return A * d ** alpha if d > 0 else 0.0


def bulk_modulus(phi, phi_c, alpha, A):
    """K = dP/d(lnΦ) = Φ·dP/dΦ = Φ·A·α·(Φ − Φ_c)^(α−1)."""
    d = phi - phi_c
    return phi * A * alpha * d ** (alpha - 1.0) if d > 0 else 0.0


# ── stages ───────────────────────────────────────────────────────────────────
def generate():
    """Single representative case: copy config.toml into sweep/isotropic/."""
    os.makedirs(SWEEP_DIR, exist_ok=True)
    case = os.path.join(SWEEP_DIR, "isotropic")
    os.makedirs(case, exist_ok=True)
    shutil.copy(CONFIG, os.path.join(case, "config.toml"))
    print(f"wrote {case}/config.toml")


def start():
    subprocess.run(
        ["cargo", "build", "--release", "--no-default-features", "--example", EXAMPLE],
        cwd=REPO_ROOT, check=True)
    exe = os.path.join(REPO_ROOT, "target", "release", "examples", EXAMPLE)
    # Wipe stale results so an old run can never be re-plotted.
    if os.path.isdir(DATA_DIR):
        shutil.rmtree(DATA_DIR)
    os.makedirs(DATA_DIR, exist_ok=True)
    subprocess.run([exe, CONFIG], cwd=REPO_ROOT, check=True)


def validate(phi, p):
    """PASS criteria: pressure rises monotonically on the jammed branch and the
    EOS fit is good (R² ≥ 0.95) with a physically sensible exponent (1 ≤ α ≤ 3)."""
    phi_j, p_j = jammed_branch(phi, p)
    phi_j, p_j = list(phi_j), list(p_j)
    print("\n  Bulk compressibility closure — validation")
    print("  " + "-" * 56)
    ok = True
    if len(phi_j) < 4:
        print(f"  jammed branch (p ≥ {P_JAM:g} Pa): too few points ({len(phi_j)})  FAIL")
        return False, None
    # Monotonic rise of p with Φ on the jammed branch (allow small noise).
    rises = sum(1 for i in range(1, len(p_j)) if p_j[i] >= p_j[i - 1])
    frac_mono = rises / (len(p_j) - 1)
    mono_ok = frac_mono >= 0.8
    ok &= mono_ok
    print(f"  jammed-branch points         : {len(phi_j)}  (Φ {min(phi_j):.3f}–{max(phi_j):.3f})")
    print(f"  monotonic p(Φ) fraction      : {frac_mono:.2f}   {'PASS' if mono_ok else 'FAIL'} (≥0.80)")
    fit = fit_eos(phi_j, p_j)
    if fit is None:
        print("  EOS fit                      : failed  FAIL")
        return False, None
    r2, phi_c, alpha, A = fit
    fit_ok = (r2 >= 0.95) and (1.0 <= alpha <= 3.0)
    ok &= fit_ok
    print(f"  EOS P=A(Φ−Φc)^α  fit         : Φc={phi_c:.3f}, α={alpha:.2f}, A={A:.3e}, R²={r2:.3f}")
    print(f"                                 {'PASS' if fit_ok else 'FAIL'} (R²≥0.95, 1≤α≤3)")
    print("  bulk modulus K = dP/d(lnΦ):")
    for pr in PHI_REPORT:
        if pr > phi_c:
            print(f"    Φ={pr:.2f}  P={eos_p(pr, phi_c, alpha, A):.3e} Pa   K={bulk_modulus(pr, phi_c, alpha, A):.3e} Pa")
    print("  " + "-" * 56)
    print(f"  RESULT: {'PASS' if ok else 'FAIL'}\n")
    return ok, fit


def graph():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(PLOT_DIR, exist_ok=True)
    if not os.path.isfile(DATA):
        print("no data — run `start` first")
        return False
    phi, p = load()
    ok, fit = validate(phi, p)

    phi_j, p_j = jammed_branch(phi, p)
    phi_j, p_j = list(phi_j), list(p_j)

    # --- EOS: P vs Φ ---
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ax.scatter(phi, p, s=10, c="lightgray", label="all samples")
    if phi_j:
        ax.scatter(phi_j, p_j, s=16, c="tab:blue", zorder=3, label="jammed branch (EOS)")
    if fit is not None:
        _, phi_c, alpha, A = fit
        xs = [phi_c + 1e-4 + 0.002 * k for k in range(0, 250)]
        xs = [x for x in xs if x <= (max(phi) + 0.01)]
        ax.plot(xs, [eos_p(x, phi_c, alpha, A) for x in xs], "k-", lw=2,
                label=fr"$P=A(\Phi-\Phi_c)^\alpha$, $\alpha$={alpha:.2f}")
        ax.axvline(phi_c, c="tab:red", ls="--", lw=1, label=fr"$\Phi_c$={phi_c:.3f}")
    ax.set_yscale("log")
    ax.set_xlabel(r"solid fraction $\Phi$")
    ax.set_ylabel(r"pressure $p$ [Pa]")
    ax.set_title("Isotropic-compression equation of state")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "eos_p_of_phi.png"), dpi=130)
    plt.close(fig)

    # --- Bulk modulus K(Φ) ---
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    if fit is not None:
        _, phi_c, alpha, A = fit
        xs = [phi_c + 1e-3 + 0.002 * k for k in range(0, 250)]
        xs = [x for x in xs if x <= (max(phi) + 0.01)]
        ax.plot(xs, [bulk_modulus(x, phi_c, alpha, A) for x in xs], "k-", lw=2,
                label=r"$K = \Phi\,dP/d\Phi$")
        for pr in PHI_REPORT:
            if pr > phi_c:
                ax.scatter([pr], [bulk_modulus(pr, phi_c, alpha, A)], c="tab:red", zorder=3)
    ax.set_xlabel(r"solid fraction $\Phi$")
    ax.set_ylabel(r"bulk modulus $K=dP/d(\ln\Phi)$ [Pa]")
    ax.set_title("Bulk compressibility closure")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "bulk_modulus.png"), dpi=130)
    plt.close(fig)

    print(f"wrote {PLOT_DIR}/eos_p_of_phi.png and bulk_modulus.png")
    return ok


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    if cmd == "generate":
        generate()
    elif cmd == "start":
        start()
    elif cmd == "graph":
        ok = graph()
        sys.exit(0 if ok else 1)
    else:
        generate()
        start()
        ok = graph()
        sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
