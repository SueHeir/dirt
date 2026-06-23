#!/usr/bin/env python3
"""
sphcal_cooling_dissipation — granular-temperature dissipation closure.

A periodic box of glass beads is given a Maxwellian velocity field and left to
cool through inelastic collisions. DIRT's contact has a velocity-independent
restitution (constant e), so the granular temperature follows Haff's law:

    T(t) = T0 / (1 + t/tc)^2          (late-time log-log slope = -2)

The cooling time tc is the dissipation closure for the de-fluidization energy
balance:  dT/dt = -omega0 sqrt(T) * T,  with the inelastic cooling rate
omega0 = 2 / (tc sqrt(T0)).  This script fits the Haff law to extract tc and
reports omega0, comparing it against the kinetic-theory dissipation coefficient
omega0_theory = (4/3) n d^2 g0 sqrt(pi) (1 - e^2)  (Carnahan-Starling g0).

The same gas is also run in LAMMPS (matched normal + tangential contact) as an
independent cross-check. An 800-particle gas is chaotic, so the two codes are
NOT compared trajectory-by-trajectory — only the cooling law (the -2 slope, the
fitted tc, and the dissipation rate omega0).

Commands (from anywhere):
    python3 .../05_cooling_dissipation/sweep.py generate   # write LAMMPS input
    python3 .../05_cooling_dissipation/sweep.py start      # build + run DIRT (+ LAMMPS)
    python3 .../05_cooling_dissipation/sweep.py graph      # validate vs Haff + plot
    python3 .../05_cooling_dissipation/sweep.py            # all three, in order

Outputs:
    sweep/in.lammps            generated LAMMPS input   (gitignored)
    data/cooling.csv           DIRT temperatures        (gitignored, written by the example)
    data/haff_trace.txt        raw LAMMPS reductions    (gitignored)
    data/lammps_cooling.csv    LAMMPS temperatures      (gitignored)
    plots/haff_cooling.png     final figure             (tracked)

Reference: P.K. Haff, "Grain flow as a fluid-mechanical phenomenon", JFM 134 (1983).
"""

import os
import sys
import csv
import math
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", ".."))
EXAMPLE = "sphcal_cooling_dissipation"

SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
CONFIG = os.path.join(SCRIPT_DIR, "config.toml")
DIRT_CSV = os.path.join(DATA_DIR, "cooling.csv")
LMP_INPUT = os.path.join(SWEEP_DIR, "in.lammps")
LMP_TRACE = os.path.join(DATA_DIR, "haff_trace.txt")
LMP_CSV = os.path.join(DATA_DIR, "lammps_cooling.csv")
PLOT = os.path.join(PLOT_DIR, "haff_cooling.png")

LAMMPS_BINS = ["lmp_serial", "lmp", "lmp_mpi", "lammps"]

# Physics — must match config.toml (canonical glass calibration material).
N = 800
L = 0.04             # box side [m]
RADIUS = 0.0011      # m
DENSITY = 2500.0     # kg/m^3
YOUNGS_MOD = 7.0e7   # Pa
POISSON = 0.245
RESTITUTION = 0.926
FRICTION = 0.16
SIGMA = 0.5          # initial Gaussian velocity sigma per component [m/s]
STEPS = 700_000
OUTPUT_EVERY = 2000  # matches the example's record() cadence
KB = 1.380649e-23    # Boltzmann constant [J/K] (for LAMMPS velocity create)

MASS = DENSITY * (4.0 / 3.0) * math.pi * RADIUS**3


def dt_rayleigh_fraction():
    """DIRT's timestep = 0.15 x Rayleigh critical timestep (solver default)."""
    g = YOUNGS_MOD / (2.0 * (1.0 + POISSON))
    alpha = 0.1631 * POISSON + 0.876605
    return 0.15 * math.pi * RADIUS / alpha * (DENSITY / g) ** 0.5


def find_lammps():
    from shutil import which
    for b in LAMMPS_BINS:
        p = which(b)
        if p:
            return p
    return None


# ── LAMMPS input ─────────────────────────────────────────────────────────────
LMP_TEMPLATE = """\
# Auto-generated LAMMPS Haff-cooling counterpart to {example}
# Matched to DIRT: Hertz/material normal + Mindlin tangential, tsuji damping (e),
# no rolling friction. Granular temperature reductions are written to a trace.
units           si
atom_style      sphere
atom_modify     map array
boundary        p p p
newton          off
comm_modify     vel yes

region          box block 0 {L} 0 {L} 0 {L} units box
create_box      1 box
create_atoms    1 random {n} 12345 box overlap {overlap} maxtry 200000
set             type 1 diameter {diam}
set             type 1 density {density}

pair_style      granular
pair_coeff      1 1 hertz/material {E} {e} {nu} tangential mindlin NULL 1.0 {mu} damping tsuji rolling none twisting none

fix             integrate all nve/sphere
velocity        all create {t_lammps:.6e} 98765 dist gaussian mom yes rot no
timestep        {dt:.10e}

compute         om all property/atom omegax omegay omegaz
variable        v2 atom vx*vx+vy*vy+vz*vz
variable        w2 atom c_om[1]*c_om[1]+c_om[2]*c_om[2]+c_om[3]*c_om[3]
compute         svx all reduce sum vx
compute         svy all reduce sum vy
compute         svz all reduce sum vz
compute         sv2 all reduce sum v_v2
compute         sw2 all reduce sum v_w2
fix             out all ave/time {every} 1 {every} c_svx c_svy c_svz c_sv2 c_sw2 file {trace} mode scalar

thermo          {thermo}
run             {steps}
"""


def write_lammps_input(dt):
    os.makedirs(SWEEP_DIR, exist_ok=True)
    os.makedirs(DATA_DIR, exist_ok=True)
    # LAMMPS `velocity create` wants a thermal temperature; pick the one that
    # yields per-component sigma: <v^2>=3 sigma^2 → T = m sigma^2 / kB.
    t_lammps = MASS * SIGMA**2 / KB
    with open(LMP_INPUT, "w") as f:
        f.write(LMP_TEMPLATE.format(
            example=EXAMPLE, L=L, n=N, overlap=2.0 * RADIUS * 1.05,
            diam=2.0 * RADIUS, density=DENSITY,
            E=f"{YOUNGS_MOD:.6e}", e=RESTITUTION, nu=POISSON, mu=FRICTION,
            t_lammps=t_lammps, dt=dt, every=OUTPUT_EVERY, thermo=10 * OUTPUT_EVERY,
            trace=LMP_TRACE, steps=STEPS,
        ))


def parse_lammps_trace(dt):
    """Convert the LAMMPS reduction trace into granular temperatures."""
    rows = []
    with open(LMP_TRACE) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 6:
                continue
            step, svx, svy, svz, sv2, sw2 = (float(x) for x in parts)
            # T_trans = <(v-<v>)^2>/3 ;  T_rot = (2/5) r^2 <w^2> / 3  (monodisperse)
            t_trans = (sv2 / N - (svx**2 + svy**2 + svz**2) / N**2) / 3.0
            t_rot = (0.4 / 3.0) * RADIUS**2 * sw2 / N
            rows.append({
                "step": int(step), "time": step * dt,
                "T_trans": t_trans, "T_rot": t_rot,
                "T_total": t_trans + t_rot,
            })
    return rows


# ── generate ─────────────────────────────────────────────────────────────────
def generate():
    write_lammps_input(dt_rayleigh_fraction())
    print(f"Generated LAMMPS input -> {LMP_INPUT}")


# ── start ────────────────────────────────────────────────────────────────────
def start():
    os.makedirs(DATA_DIR, exist_ok=True)
    print(f"Building {EXAMPLE} (release)...", flush=True)
    subprocess.run(
        ["cargo", "build", "--release", "--example", EXAMPLE, "--no-default-features"],
        cwd=REPO_ROOT, check=True,
    )

    print(f"Running DIRT ({N} spheres, {STEPS} steps)...", flush=True)
    if os.path.exists(DIRT_CSV):
        os.remove(DIRT_CSV)
    log = os.path.join(DATA_DIR, "dirt_run.log")
    with open(log, "w") as lf:
        proc = subprocess.run(
            ["cargo", "run", "--release", "--example", EXAMPLE,
             "--no-default-features", "--", CONFIG],
            cwd=REPO_ROOT, stdout=lf, stderr=subprocess.STDOUT,
        )
    if proc.returncode != 0 or not os.path.isfile(DIRT_CSV):
        print(f"ERROR: DIRT run failed (see {log}).")
        sys.exit(1)
    # Use DIRT's actual timestep for LAMMPS so the time axes coincide.
    with open(DIRT_CSV) as f:
        rows = list(csv.DictReader(f))
    dt = next((float(r["time"]) / float(r["step"]) for r in rows if int(r["step"]) > 0),
              dt_rayleigh_fraction())
    print(f"  DIRT done. dt = {dt:.6e} s, {len(rows)} samples.")

    lammps = find_lammps()
    if not lammps:
        print("LAMMPS not found on PATH — skipping the cross-check (DIRT only).")
        return
    print(f"Running LAMMPS ({lammps})...", flush=True)
    write_lammps_input(dt)
    if os.path.exists(LMP_TRACE):
        os.remove(LMP_TRACE)
    lmp_log = os.path.join(DATA_DIR, "lammps.log")
    proc = subprocess.run([lammps, "-in", LMP_INPUT, "-log", lmp_log],
                          cwd=REPO_ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    if proc.returncode != 0 or not os.path.isfile(LMP_TRACE):
        print(f"WARNING: LAMMPS run failed (see {lmp_log}); continuing DIRT-only.")
        return
    lrows = parse_lammps_trace(dt)
    with open(LMP_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["step", "time", "T_trans", "T_rot", "T_total"])
        w.writeheader()
        w.writerows(lrows)
    print(f"  LAMMPS done. {len(lrows)} samples -> {LMP_CSV}")


# ── graph (validate + plot) ──────────────────────────────────────────────────
def load(path, cols):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return {c: [float(r[c]) for r in rows] for c in cols}


def haff_fit(t, T):
    """Fit Haff's law via its linearization: 1/sqrt(T) = 1/sqrt(T0) + t/(sqrt(T0) tc).
    Returns (T0, tc, r2, slope_loglog) using the clean cooling window (T/T0 > 1e-3)."""
    import numpy as np
    t = np.asarray(t)
    T = np.asarray(T)
    T0_raw = T[0] if T[0] > 0 else T[T > 0][0]
    win = (t >= 0) & (T > 1e-3 * T0_raw) & np.isfinite(T)
    tw, Tw = t[win], T[win]
    y = 1.0 / np.sqrt(Tw)
    b, a = np.polyfit(tw, y, 1)          # y = b*t + a
    T0 = 1.0 / a**2
    tc = a / b
    resid = y - (b * tw + a)
    r2 = 1.0 - np.sum(resid**2) / np.sum((y - y.mean())**2)
    # late-time log-log slope (second half of the window, t>0)
    pos = tw > 0
    tp, Tp = tw[pos], Tw[pos]
    half = slice(len(tp) // 2, None)
    slope = np.polyfit(np.log(tp[half]), np.log(Tp[half]), 1)[0]
    return T0, tc, r2, slope


def dissipation_theory():
    """Kinetic-theory inelastic cooling rate omega0 (Carnahan-Starling g0).
    Energy balance: dT/dt = -omega0 sqrt(T) T, so tc_theory = 2/(omega0 sqrt(T0))."""
    n = N / L**3
    phi = n * (4.0 / 3.0) * math.pi * RADIUS**3
    g0 = (1.0 - phi / 2.0) / (1.0 - phi) ** 3
    omega0 = (4.0 / 3.0) * n * (2 * RADIUS)**2 * g0 * math.sqrt(math.pi) * (1 - RESTITUTION**2)
    return phi, g0, omega0


def validate(dirt):
    import numpy as np
    t = np.array(dirt["time"])
    Tt = np.array(dirt["T_total"])
    print("=" * 64)
    print("Granular-Temperature Dissipation Closure — Glass Beads (DIRT)")
    print("=" * 64)

    total = passed = 0

    def check(name, ok, detail=""):
        nonlocal total, passed
        total += 1
        passed += bool(ok)
        print(f"  {name:<28}{'PASS' if ok else 'FAIL'}   {detail}")

    check("finite temperatures", bool(np.all(np.isfinite(Tt))))
    check("non-negative T", bool(np.all(Tt >= 0)))
    check("cooling (Tf < Ti)", Tt[-1] < Tt[0],
          f"Ti={Tt[0]:.3e}  Tf={Tt[-1]:.3e}")
    check("no energy growth", float(np.max(Tt)) < 1.5 * Tt[0],
          f"max={np.max(Tt):.3e}")

    T0, tc, r2, slope = haff_fit(t, Tt)
    check("Haff law (1/sqrt(T) linear)", r2 > 0.99, f"R^2={r2:.4f}")
    check("late-time slope ~ -2", -2.3 < slope < -1.7, f"slope={slope:.3f}")

    # Extracted dissipation closure: omega0 = 2 / (tc sqrt(T0)).
    omega0_fit = 2.0 / (tc * math.sqrt(T0)) if (tc > 0 and T0 > 0) else float("nan")
    phi, g0, omega0_th = dissipation_theory()
    tc_theory = 2.0 / (omega0_th * math.sqrt(Tt[0]))
    print(f"\n  Haff fit:        T0={T0:.3e}  tc={tc:.3e} s")
    print(f"  dissipation rate omega0_fit   = {omega0_fit:.3e} (m^2/s^2)^-1/2 s^-1")
    print(f"  kinetic theory   omega0_theory= {omega0_th:.3e}  "
          f"(phi={phi:.3f}, g0={g0:.3f}, 1-e^2={1-RESTITUTION**2:.4f})")
    print(f"  tc_theory={tc_theory:.3e} s   (fit/theory: tc={tc/tc_theory:.2f}, "
          f"omega0={omega0_fit/omega0_th:.2f})")

    print(f"\nResult: {passed}/{total} checks passed")
    print("ALL CHECKS PASSED" if passed == total
          else f"WARNING: {total - passed} check(s) failed")
    return passed == total


def graph():
    if not os.path.isfile(DIRT_CSV):
        print(f"ERROR: {DIRT_CSV} not found. Run 'start' first.")
        sys.exit(1)
    dirt = load(DIRT_CSV, ["time", "T_trans", "T_rot", "T_total"])
    ok = validate(dirt)
    lammps = load(LMP_CSV, ["time", "T_trans", "T_rot", "T_total"]) \
        if os.path.isfile(LMP_CSV) else None

    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"\n(matplotlib/numpy unavailable, skipped plot: {e})")
        return ok

    plt.rcParams.update({"font.size": 11, "axes.labelsize": 12,
                         "figure.dpi": 150, "savefig.dpi": 150})
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # ── Panel 1: cooling-law comparison (normalized, log-log) ──────────────
    ax = axes[0]
    td = np.array(dirt["time"]); Td = np.array(dirt["T_total"])
    md = td > 0
    T0d, tcd, r2d, sd = haff_fit(td, Td)
    omega0_fit = 2.0 / (tcd * math.sqrt(T0d))
    ax.loglog(td[md], Td[md] / Td[0], "o", ms=3, color="#1f77b4", alpha=0.6,
              label=f"DIRT  (slope {sd:.2f})")
    if lammps:
        tl = np.array(lammps["time"]); Tl = np.array(lammps["T_total"])
        ml = tl > 0
        _, _, _, sl = haff_fit(tl, Tl)
        ax.loglog(tl[ml], Tl[ml] / Tl[0], "s", ms=3, color="#ff7f0e", alpha=0.6,
                  label=f"LAMMPS  (slope {sl:.2f})")
    # Haff fit (DIRT) and -2 reference
    tf = np.linspace(td[md][0], td[md][-1], 300)
    ax.loglog(tf, (T0d / Td[0]) / (1 + tf / tcd) ** 2, "-", color="black", lw=1.5,
              label=f"Haff fit  tc={tcd:.2e}s\nomega0={omega0_fit:.2e}")
    tref = np.array([td[md][len(td[md]) // 3], td[md][-1]])
    yref = (Td[md][len(td[md]) // 3] / Td[0]) * (tref / tref[0]) ** -2
    ax.loglog(tref, yref, "--", color="gray", lw=1.2, label="slope -2 (Haff)")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$T_\mathrm{total}/T_0$")
    ax.set_title("Cooling law / dissipation closure")
    ax.legend(fontsize=8)

    # ── Panel 2: DIRT energy partition (rough spheres) ─────────────────────
    ax = axes[1]
    ax.semilogy(dirt["time"], dirt["T_trans"], color="#1f77b4", label=r"$T_\mathrm{trans}$")
    ax.semilogy(dirt["time"], dirt["T_rot"], color="#d62728", label=r"$T_\mathrm{rot}$")
    ax.semilogy(dirt["time"], dirt["T_total"], color="black", label=r"$T_\mathrm{total}$")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"Granular temperature [m$^2$/s$^2$]")
    ax.set_title("DIRT energy partition (rough glass beads)")
    ax.legend(fontsize=8)

    fig.suptitle("Granular-Temperature Dissipation — Glass Beads", y=1.02)
    fig.tight_layout()
    os.makedirs(PLOT_DIR, exist_ok=True)
    fig.savefig(PLOT, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved: {PLOT}")
    return ok


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    if cmd == "generate":
        generate()
    elif cmd == "start":
        start()
    elif cmd == "graph":
        sys.exit(0 if graph() else 1)
    elif cmd == "all":
        generate()
        start()
        print()
        sys.exit(0 if graph() else 1)
    else:
        print(f"Unknown command: {cmd!r}")
        print("Usage: sweep.py [generate|start|graph]   (no arg = all three)")
        sys.exit(2)


if __name__ == "__main__":
    main()
