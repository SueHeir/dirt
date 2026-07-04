#!/usr/bin/env python3
"""Independent LAMMPS cross-validation of the Hertz normal contact using the
CLASSIC `gran/hertz/history` pair/wall style — a *different* LAMMPS code path
from the `pair_style granular hertz/material` overlay that `sweep.py` already
runs. Running both is deliberate: `gran/hertz/history` (src/GRANULAR/
pair_gran_hertz_history.cpp) and the modern `granular` model (granular_model.cpp)
are separate implementations, so agreement is genuine independent-code provenance
for DIRT's Hertz contact rather than a self-consistency check.

What it validates against DIRT's `bench_hertz_rebound` (data/sweep_results.csv):
  * peak overlap   δ_max
  * contact time   t_c
  * COR            v_rebound / v_impact

Method
------
`gran/hertz/history` uses (LAMMPS doc/src/pair_gran.rst):

    F_hz = sqrt(δ)·sqrt(Ri·Rj/(Ri+Rj))·[ Kn·δ·n_ij − m_eff·γ_n·v_n − (tangential) ]

with the documented material mapping Kn = 4G/(3(1−ν)) = 2E/(3(1−ν²)). For a
sphere on a flat wall (Rj→∞) the geometric prefactor is sqrt(R), so the elastic
term is

    F = Kn·sqrt(R)·δ^{3/2} = (4/3)·E*·sqrt(R)·δ^{3/2},   E* = E/(2(1−ν²)),

i.e. it is IDENTICAL to true Hertz and to DIRT's spring (DIRT uses the same
E* = E/(2(1−ν²)) convention). So Kn = (4/3)·E* and contact-time / peak-overlap
are geometry+stiffness predictions, NOT fitted.

Two comparisons are made:

  1. NEAR-ELASTIC ANCHOR, all four impact speeds. The `gran/hertz/history`
     classic parser sets the tangential coeff to `gammat/gamman`
     (granular_model.cpp), so an EXACTLY-zero normal damping is a 0/0 NaN and is
     not representable. We therefore use the smallest γ_n that realizes COR ≥
     0.9999 (a negligible dashpot): δ_max and t_c then sit within <0.5 % of the
     undamped Hertz closed form, so `gran/hertz/history` must reproduce that form
     AND DIRT's COR≈1.0 row for δ_max, t_c and COR. This is the strongest,
     tolerance-bounded, nearly-damping-free cross-check of the contact stiffness.

  2. DAMPED COR (e = 0.7, 0.9). `gran/hertz/history` damps with a CONSTANT γ_n
     (an older viscoelastic dashpot), NOT the Tsuji polynomial DIRT/`granular`
     use. The two damping *laws* differ: a constant-γ_n Hertz dashpot realizes a
     weakly velocity-dependent COR (∝ v0^{1/5}), whereas Tsuji is velocity-
     independent — a documented, expected model difference, not a code error.
     We therefore calibrate the ONE free damping knob γ_n (via an independent
     RK4 integration of the contact ODE) so `gran/hertz/history` realizes DIRT's
     measured COR at the reference speed v0 = 1.0 m/s, confirm the LAMMPS run
     reproduces that COR, and report that the geometry-determined δ_max and t_c
     STILL match DIRT to tolerance (they are not tuned).

Nothing here weakens the DIRT benchmark: all DIRT reference numbers and its own
tolerances are untouched; this is an additive, independent overlay.

Run:  $BENCH_PYTHON examples/bench_hertz_rebound/xval_gran_hertz_history.py
Reads data/sweep_results.csv (DIRT); writes data/xval_gran_hertz_history.csv.
Exit 0 = all cross-checks within tolerance.
"""

import os
import sys
import csv
import math
import shutil
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
import sweep  # reuse the exact constants + trace parser DIRT's bench uses

REPO_ROOT = sweep.REPO_ROOT
DATA_DIR = sweep.DATA_DIR
XVAL_DIR = os.path.join(SCRIPT_DIR, "xval")
DIRT_CSV = sweep.SWEEP_CSV
OUT_CSV = os.path.join(DATA_DIR, "xval_gran_hertz_history.csv")

# Physical constants — identical to the DIRT bench.
E = sweep.YOUNGS_MOD
NU = sweep.POISSON_RATIO
R = sweep.RADIUS
DENSITY = sweep.DENSITY
E_STAR = sweep.E_STAR          # E/(2(1-nu^2)) — DIRT / sphere-on-same-material-flat
MASS = sweep.MASS
M_EFF = sweep.M_EFF

# gran/hertz/history normal elastic constant (LAMMPS material mapping).
# Kn = 4G/(3(1-nu)) = 2E/(3(1-nu^2)) = (4/3) E*  ->  identical Hertz spring.
G = E / (2.0 * (1.0 + NU))
KN = 4.0 * G / (3.0 * (1.0 - NU))
KN_CHECK = (4.0 / 3.0) * E_STAR    # must equal KN
KT = 4.0 * G / (2.0 - NU)          # documented tangential constant (unused: xmu=0)

VELOCITIES = sweep.VELOCITIES      # [0.1, 0.5, 1.0, 2.0]
REF_V0 = 1.0                       # calibration reference speed for damped cases
DAMPED_CORS = [0.7, 0.9]           # nominal restitutions for the damped comparison
ANCHOR_COR = 0.9999                # near-elastic target (exact 0 damping -> NaN, see below)

# ── tolerances (stated, not back-fitted) ────────────────────────────────────
OVERLAP_TOL = 0.01     # |δ_max_ghh - δ_max_dirt| / δ_max_dirt
CONTACT_TOL = 0.03     # |t_c_ghh - t_c_dirt| / t_c_dirt  (~1 shared timestep)
THEORY_TOL = 0.02      # δ_max / t_c vs analytic Hertz closed form
COR_ELASTIC_TOL = 0.005    # |COR_ghh - 1| at the near-elastic anchor
COR_XCODE_TOL = 0.01       # |COR_ghh - COR_dirt| at the reference speed (damped)


# ── gran/hertz/history 1-DOF normal contact ODE (independent of LAMMPS) ──────
# m δ̈ = -sqrt(R)·[ Kn·δ^{3/2} + m·γ_n·δ^{1/2}·δ̇ ]   during contact (δ>0).
# Used only to CALIBRATE γ_n to a target COR; LAMMPS then verifies it.
def ode_cor(gamma_n, v0):
    sr = math.sqrt(R)
    kterm = sr * KN / MASS          # δ̈ elastic coeff
    cterm = sr * gamma_n            # δ̈ damping coeff (m cancels)

    def acc(d, v):
        d = max(d, 0.0)
        return -(kterm * d ** 1.5 + cterm * math.sqrt(d) * v)

    # step ~ contact_time / 4000, robustly small
    tc = sweep.hertz_contact_duration(v0)
    dt = tc / 4000.0
    d, v, t = 0.0, v0, 0.0
    # advance one step off zero so sqrt(δ) damping engages
    while t < 10.0 * tc:
        a1 = acc(d, v)
        a2 = acc(d + 0.5 * dt * v, v + 0.5 * dt * a1)
        a3 = acc(d + 0.5 * dt * (v + 0.5 * dt * a1), v + 0.5 * dt * a2)
        a4 = acc(d + dt * (v + 0.5 * dt * a2), v + dt * a3)
        d_new = d + dt * (v + dt / 6.0 * (a1 + a2 + a3))
        v_new = v + dt / 6.0 * (a1 + 2 * a2 + 2 * a3 + a4)
        if d > 0.0 and d_new <= 0.0:
            return abs(v_new) / v0
        d, v, t = d_new, v_new, t + dt
    return 0.0


def calibrate_gamma_n(target_cor, v0):
    """Bisection: COR decreases monotonically with γ_n. Returns γ_n for target."""
    lo, hi = 0.0, 1.0
    while ode_cor(hi, v0) > target_cor:   # grow bracket until COR falls below target
        hi *= 2.0
        if hi > 1e12:
            return hi
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if ode_cor(mid, v0) > target_cor:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ── LAMMPS input for gran/hertz/history (classic pair + wall) ────────────────
LMP_TEMPLATE = """\
# Independent cross-validation input: gran/hertz/history (classic GRANULAR path)
# v0 = {v0} m/s, gamma_n = {gamma_n:.6e}, label = {label}
units           si
atom_style      sphere
atom_modify     map array
dimension       3
boundary        f f f
newton          off
comm_modify     vel yes

region          simbox block -0.01 0.01 -0.01 0.01 0.0 0.1 units box
create_box      1 simbox

create_atoms    1 single 0.0 0.0 0.0075 units box
set             group all diameter {diam}
set             group all density {density}

# Kn = 4G/(3(1-nu)) = (4/3)E*  (material mapping, doc/src/pair_gran.rst)
# Kt=NULL -> 2/7 Kn, gamma_t=NULL -> gamma_n/2, xmu=0 (frictionless), dampflag=0
pair_style      gran/hertz/history {Kn:.10e} NULL {gamma_n:.10e} NULL 0.0 0
pair_coeff      * *

fix             wall all wall/gran hertz/history {Kn:.10e} NULL {gamma_n:.10e} NULL 0.0 0 zplane 0.0 NULL
fix             integrate all nve/sphere

velocity        all set 0.0 0.0 -{v0} units box
timestep        {dt}

variable        tnow equal time
variable        zpos equal z[1]
variable        zvel equal vz[1]
fix             rec all print 1 "${{tnow}} ${{zpos}} ${{zvel}}" file {trace} screen no title "t z vz"

thermo          5000
run             {steps}
"""


def find_lammps():
    for b in sweep.LAMMPS_BINS:
        p = shutil.which(b)
        if p:
            return p
    return None


def run_case(lammps, label, v0, gamma_n):
    cdir = os.path.join(XVAL_DIR, label)
    os.makedirs(cdir, exist_ok=True)
    in_path = os.path.join(cdir, "in.lammps")
    log_path = os.path.join(cdir, "lammps.log")
    trace = os.path.join(cdir, "lammps_trace.txt")
    dt = sweep.dt_for(v0)
    with open(in_path, "w") as f:
        f.write(LMP_TEMPLATE.format(
            v0=v0, gamma_n=gamma_n, label=label,
            diam=2.0 * R, density=DENSITY, Kn=KN,
            dt=f"{dt:.10e}", steps=sweep.steps_for(v0), trace=trace,
        ))
    proc = subprocess.run([lammps, "-in", in_path, "-log", log_path],
                          cwd=REPO_ROOT, stdout=subprocess.DEVNULL,
                          stderr=subprocess.STDOUT)
    if proc.returncode != 0 or not os.path.isfile(trace):
        return None
    return sweep.parse_rebound_trace(trace, R, dt)


# ── DIRT reference rows ──────────────────────────────────────────────────────
def load_dirt():
    rows = {}
    with open(DIRT_CSV) as f:
        for r in csv.DictReader(f):
            rows[(float(r["input_v0"]), float(r["input_cor"]))] = r
    return rows


def rel(a, b):
    return abs(a - b) / abs(b) if b else float("inf")


def main():
    if abs(KN - KN_CHECK) / KN_CHECK > 1e-12:
        print(f"FAIL: Kn mapping inconsistent: {KN:.6e} vs (4/3)E*={KN_CHECK:.6e}")
        return 1
    lammps = find_lammps()
    if not lammps:
        print("FAIL: no LAMMPS binary on PATH (lmp/lmp_serial/...).")
        return 1
    if not os.path.isfile(DIRT_CSV):
        print(f"FAIL: DIRT results not found: {DIRT_CSV}\n"
              f"  run: $BENCH_PYTHON examples/bench_hertz_rebound/sweep.py start")
        return 1

    print(f"LAMMPS: {lammps}")
    print(f"Kn = 4G/(3(1-nu)) = {KN:.6e} Pa  ( = (4/3)E*, E*={E_STAR:.6e} )")
    os.makedirs(XVAL_DIR, exist_ok=True)
    dirt = load_dirt()
    out_rows = []
    failures = []

    # ---- 1. NEAR-ELASTIC ANCHOR (smallest gamma_n; exact 0 -> 0/0 NaN) --------
    # granular_model.cpp sets the classic-hertz tangential coeff = gammat/gamman,
    # so gamma_n=0 is a 0/0 NaN. Use the smallest gamma_n giving COR>=0.9999
    # (calibrated at the fastest speed, where const-gamma COR is lowest, so every
    # speed lands >=0.9999). Damping is then negligible: dmax/tc ~ undamped Hertz.
    anchor_gamma = calibrate_gamma_n(ANCHOR_COR, max(VELOCITIES))
    print(f"\n=== Near-elastic anchor: gran/hertz/history (gamma_n={anchor_gamma:.3e}, COR>=~{ANCHOR_COR}) "
          f"vs DIRT COR=1.0 & Hertz theory ===")
    print(f"{'v0':>5} {'COR_ghh':>9} {'dmax_ghh':>11} {'dmax_dirt':>11} {'dmax_th':>11} "
          f"{'tc_ghh':>11} {'tc_dirt':>11} {'tc_th':>11}")
    for v0 in VELOCITIES:
        res = run_case(lammps, f"elastic_v{v0}", v0, anchor_gamma)
        if not res:
            failures.append(f"elastic v0={v0}: LAMMPS run/parse failed")
            continue
        d = dirt.get((v0, 1.0))
        dmax_dirt = float(d["max_overlap"]); tc_dirt = float(d["contact_time"])
        dmax_th = sweep.hertz_max_overlap(res["v_impact"])
        tc_th = sweep.hertz_contact_duration(res["v_impact"])
        cor = res["cor_measured"]; dmax = res["max_overlap"]; tc = res["contact_time"]
        print(f"{v0:>5} {cor:>9.5f} {dmax:>11.4e} {dmax_dirt:>11.4e} {dmax_th:>11.4e} "
              f"{tc:>11.4e} {tc_dirt:>11.4e} {tc_th:>11.4e}")
        out_rows.append(dict(case="elastic", input_v0=v0, input_cor=1.0, gamma_n=anchor_gamma,
                             cor_ghh=cor, cor_dirt=float(d["cor_measured"]),
                             dmax_ghh=dmax, dmax_dirt=dmax_dirt, dmax_theory=dmax_th,
                             tc_ghh=tc, tc_dirt=tc_dirt, tc_theory=tc_th))
        if abs(cor - 1.0) > COR_ELASTIC_TOL:
            failures.append(f"elastic v0={v0}: COR {cor:.5f} not ~1 (tol {COR_ELASTIC_TOL})")
        if rel(dmax, dmax_dirt) > OVERLAP_TOL:
            failures.append(f"elastic v0={v0}: dmax vs DIRT {rel(dmax,dmax_dirt)*100:.2f}% > {OVERLAP_TOL*100:.0f}%")
        if rel(dmax, dmax_th) > THEORY_TOL:
            failures.append(f"elastic v0={v0}: dmax vs theory {rel(dmax,dmax_th)*100:.2f}% > {THEORY_TOL*100:.0f}%")
        if rel(tc, tc_dirt) > CONTACT_TOL:
            failures.append(f"elastic v0={v0}: tc vs DIRT {rel(tc,tc_dirt)*100:.2f}% > {CONTACT_TOL*100:.0f}%")
        if rel(tc, tc_th) > THEORY_TOL:
            failures.append(f"elastic v0={v0}: tc vs theory {rel(tc,tc_th)*100:.2f}% > {THEORY_TOL*100:.0f}%")

    # ---- 2. DAMPED COR (calibrated gamma_n; documented model difference) -------
    print("\n=== Damped COR: gamma_n calibrated to DIRT's COR at v0=1.0 (ODE), verified in LAMMPS ===")
    for cor_nom in DAMPED_CORS:
        d_ref = dirt.get((REF_V0, cor_nom))
        target = float(d_ref["cor_measured"])         # DIRT's realized COR at v0=1
        gamma_n = calibrate_gamma_n(target, REF_V0)
        print(f"\n nominal e={cor_nom}: DIRT realized COR@v0=1 = {target:.5f}"
              f"  ->  calibrated gamma_n = {gamma_n:.4e}  (ODE COR={ode_cor(gamma_n,REF_V0):.5f})")
        print(f"{'v0':>5} {'COR_ghh':>9} {'COR_dirt':>9} {'dmax_ghh':>11} {'dmax_dirt':>11} "
              f"{'tc_ghh':>11} {'tc_dirt':>11}")
        cors = []
        for v0 in VELOCITIES:
            res = run_case(lammps, f"damped_e{cor_nom}_v{v0}", v0, gamma_n)
            if not res:
                failures.append(f"damped e={cor_nom} v0={v0}: LAMMPS run/parse failed")
                continue
            dd = dirt.get((v0, cor_nom))
            cor = res["cor_measured"]; dmax = res["max_overlap"]; tc = res["contact_time"]
            cor_dirt = float(dd["cor_measured"]); dmax_dirt = float(dd["max_overlap"]); tc_dirt = float(dd["contact_time"])
            cors.append((v0, cor))
            print(f"{v0:>5} {cor:>9.5f} {cor_dirt:>9.5f} {dmax:>11.4e} {dmax_dirt:>11.4e} "
                  f"{tc:>11.4e} {tc_dirt:>11.4e}")
            out_rows.append(dict(case="damped", input_v0=v0, input_cor=cor_nom, gamma_n=gamma_n,
                                 cor_ghh=cor, cor_dirt=cor_dirt,
                                 dmax_ghh=dmax, dmax_dirt=dmax_dirt, dmax_theory="",
                                 tc_ghh=tc, tc_dirt=tc_dirt, tc_theory=""))
            # Assert agreement ONLY at the reference speed v0=1.0 (the calibration
            # point). At other speeds the constant-gamma_n COR diverges from Tsuji
            # by the documented ~v0^{1/5} model difference, and dmax/tc diverge
            # consistently with it — that divergence is reported, not asserted.
            if v0 == REF_V0:
                if abs(cor - cor_dirt) > COR_XCODE_TOL:
                    failures.append(f"damped e={cor_nom} v0={v0}: COR {cor:.5f} vs DIRT {cor_dirt:.5f} > {COR_XCODE_TOL}")
                if rel(dmax, dmax_dirt) > OVERLAP_TOL:
                    failures.append(f"damped e={cor_nom} v0={v0}: dmax vs DIRT {rel(dmax,dmax_dirt)*100:.2f}% > {OVERLAP_TOL*100:.0f}%")
                if rel(tc, tc_dirt) > CONTACT_TOL:
                    failures.append(f"damped e={cor_nom} v0={v0}: tc vs DIRT {rel(tc,tc_dirt)*100:.2f}% > {CONTACT_TOL*100:.0f}%")
        if cors:
            spread = max(c for _, c in cors) - min(c for _, c in cors)
            print(f"   COR velocity spread (gran/hertz/history, const gamma_n): {spread:.4f}"
                  f"   [expected nonzero ~v0^1/5; DIRT/Tsuji ~0 — documented model difference]")

    # ---- write results + verdict ----------------------------------------------
    os.makedirs(DATA_DIR, exist_ok=True)
    fields = ["case", "input_v0", "input_cor", "gamma_n", "cor_ghh", "cor_dirt",
              "dmax_ghh", "dmax_dirt", "dmax_theory", "tc_ghh", "tc_dirt", "tc_theory"]
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in out_rows:
            w.writerow(r)
    print(f"\nWrote {len(out_rows)} rows -> {OUT_CSV}")

    if failures:
        print(f"\nFAIL ({len(failures)}):")
        for m in failures:
            print(f"  - {m}")
        return 1
    print("\nPASS: gran/hertz/history matches DIRT bench_hertz_rebound within stated tolerances.")
    print(f"  near-elastic anchor (COR>=~{ANCHOR_COR}): dmax within {OVERLAP_TOL*100:.0f}% (vs DIRT) / {THEORY_TOL*100:.0f}% (vs theory),")
    print(f"  tc within {CONTACT_TOL*100:.0f}% (vs DIRT) / {THEORY_TOL*100:.0f}% (vs theory);")
    print(f"  damped COR/dmax/tc at v0={REF_V0} within tol of DIRT for e in {DAMPED_CORS}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
