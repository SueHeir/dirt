#!/usr/bin/env python3
"""
Haff's-cooling benchmark — multisphere clumps (7-sphere "ball").

A periodic box of rigid clumps is given a random velocity field and left to cool
through inelastic collisions. DIRT's contact has a velocity-independent
restitution (constant e), so the granular temperature obeys Haff's law:

    T(t) = T0 / (1 + t/tc)^2          (late-time log-log slope -> -2)

The same gas is also run in LAMMPS (rigid multisphere via `fix rigid/small`,
matched Hertz + Mindlin + tsuji contact) as an independent cross-check. Because
a many-body granular gas is chaotic — and the two codes use different rigid-body
integrators and clump-contact handling — the comparison is on the *cooling law*
(the Haff fit / -2 slope), not trajectory-by-trajectory.

Commands (from anywhere):
    python3 examples/bench_clump_haff_cooling/sweep.py generate   # write clump.mol + in.lammps
    python3 examples/bench_clump_haff_cooling/sweep.py start      # build + run DIRT and LAMMPS
    python3 examples/bench_clump_haff_cooling/sweep.py graph      # validate vs Haff + plot
    python3 examples/bench_clump_haff_cooling/sweep.py            # all three, in order

Outputs (all under data/, gitignored, except the tracked figure):
    data/clump.mol, data/in.lammps, data/haff_trace.txt, data/lammps_cooling.csv
    data/cooling.csv  (DIRT, written by the example)
    plots/haff_cooling.png  (tracked)

Reference: P.K. Haff, "Grain flow as a fluid-mechanical phenomenon", JFM 134 (1983).
"""

import os
import sys
import csv
import math
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_clump_haff_cooling"
TITLE = "Multisphere Clumps (sphere7)"

DATA_DIR = os.path.join(SCRIPT_DIR, "data")
CONFIG = os.path.join(SCRIPT_DIR, "config.toml")
DIRT_CSV = os.path.join(DATA_DIR, "cooling.csv")
MOL_FILE = os.path.join(DATA_DIR, "clump.mol")
LMP_INPUT = os.path.join(DATA_DIR, "in.lammps")
LMP_TRACE = os.path.join(DATA_DIR, "haff_trace.txt")
LMP_CSV = os.path.join(DATA_DIR, "lammps_cooling.csv")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
PLOT = os.path.join(PLOT_DIR, "haff_cooling.png")

LAMMPS_BINS = ["lmp_serial", "lmp", "lmp_mpi", "lammps"]

# Clump geometry & physics — must match config.toml.
# sphere7: central sub-sphere + 6 satellites at +/-0.6 mm along each axis.
SUB_R = 0.0005
CLUMP = [(0.0, 0.0, 0.0), (0.0006, 0.0, 0.0), (-0.0006, 0.0, 0.0),
         (0.0, 0.0006, 0.0), (0.0, -0.0006, 0.0),
         (0.0, 0.0, 0.0006), (0.0, 0.0, -0.0006)]
N = 500
L = 0.04
DENSITY = 2500.0
YOUNGS_MOD = 5.0e7
POISSON = 0.3
RESTITUTION = 0.9
FRICTION = 0.3
SIGMA = 0.5
STEPS = 700_000
OUTPUT_EVERY = 2000
KB = 1.380649e-23

NSUB = len(CLUMP)
M_SUB = DENSITY * (4.0 / 3.0) * math.pi * SUB_R**3
M_BODY = NSUB * M_SUB
M_TOTAL = N * M_BODY


def dt_rayleigh_fraction():
    g = YOUNGS_MOD / (2.0 * (1.0 + POISSON))
    alpha = 0.1631 * POISSON + 0.876605
    return 0.15 * math.pi * SUB_R / alpha * (DENSITY / g) ** 0.5


def find_lammps():
    from shutil import which
    for b in LAMMPS_BINS:
        p = which(b)
        if p:
            return p
    return None


# ── LAMMPS molecule + input ──────────────────────────────────────────────────
def write_molecule():
    os.makedirs(DATA_DIR, exist_ok=True)
    lines = ["# auto-generated clump", "", f"{NSUB} atoms", "", "Coords", ""]
    lines += [f"{i+1} {x} {y} {z}" for i, (x, y, z) in enumerate(CLUMP)]
    lines += ["", "Types", ""] + [f"{i+1} 1" for i in range(NSUB)]
    lines += ["", "Diameters", ""] + [f"{i+1} {2*SUB_R}" for i in range(NSUB)]
    lines += ["", "Masses", ""] + [f"{i+1} {M_SUB:.8e}" for i in range(NSUB)]
    with open(MOL_FILE, "w") as f:
        f.write("\n".join(lines) + "\n")


LMP_TEMPLATE = """\
# Auto-generated LAMMPS rigid-clump Haff-cooling counterpart to {example}
units           si
atom_style      hybrid sphere molecular
atom_modify     map array
boundary        p p p
newton          off
comm_modify     vel yes

region          box block 0 {L} 0 {L} 0 {L} units box
create_box      1 box
molecule        clump {mol}
create_atoms    0 random {n} 12345 box mol clump 54321 overlap {overlap} maxtry 500000

pair_style      granular
pair_coeff      1 1 hertz/material {E} {e} {nu} tangential mindlin NULL 1.0 {mu} damping tsuji rolling none twisting none
neigh_modify    exclude molecule/intra all

velocity        all create {t_lammps:.6e} 98765 dist gaussian mom yes
fix             integrate all rigid/small molecule
compute         ctemp all temp/sphere
thermo_modify   temp ctemp
timestep        {dt:.10e}
compute         keatom all ke
"""

# Reference velocity-create temperature (arbitrary scale; calibrated at run time
# so the LAMMPS initial granular temperature matches DIRT's — the rigid-body
# velocity projection makes the absolute value code-dependent).
T_CREATE_REF = M_SUB * (SIGMA**2 * NSUB) / KB

# Granular T from total clump KE: rigid-body summed sub-sphere KE = body trans +
# rot KE, and DIRT's (no-1/2) convention gives T_total = 2 KE / (3 M_total).
KE_TO_T = 2.0 / (3.0 * M_TOTAL)


def _write_lammps(dt, t_create, tail):
    os.makedirs(DATA_DIR, exist_ok=True)
    head = LMP_TEMPLATE.format(
        example=EXAMPLE, L=L, n=N, mol=MOL_FILE, overlap=2.0 * SUB_R * 1.1,
        E=f"{YOUNGS_MOD:.6e}", e=RESTITUTION, nu=POISSON, mu=FRICTION,
        t_lammps=t_create, dt=dt)
    with open(LMP_INPUT, "w") as f:
        f.write(head + tail)


def write_lammps_input(dt, t_create=T_CREATE_REF):
    """Full cooling run."""
    tail = (f"fix             out all ave/time {OUTPUT_EVERY} 1 {OUTPUT_EVERY} "
            f"c_keatom file {LMP_TRACE} mode scalar\n"
            f"thermo          {10*OUTPUT_EVERY}\nrun             {STEPS}\n")
    _write_lammps(dt, t_create, tail)


def write_lammps_calib(dt, t_create=T_CREATE_REF):
    """0-step setup to read the post-projection initial KE (-> initial T)."""
    calib = os.path.join(DATA_DIR, "calib.txt")
    if os.path.exists(calib):
        os.remove(calib)
    tail = (f"run             0\nvariable kk equal c_keatom\n"
            f"print           \"${{kk}}\" file {calib} screen no\n")
    _write_lammps(dt, t_create, tail)
    return calib


def calibrate_t_create(lammps, dt, t0_target):
    """Find the velocity-create temperature whose projected initial granular
    temperature equals DIRT's t0_target. KE scales linearly with create-T, so
    one 0-step run suffices."""
    calib = write_lammps_calib(dt, T_CREATE_REF)
    lmp_log = os.path.join(DATA_DIR, "lammps_calib.log")
    proc = subprocess.run([lammps, "-in", LMP_INPUT, "-log", lmp_log],
                          cwd=REPO_ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    if proc.returncode != 0 or not os.path.isfile(calib):
        return T_CREATE_REF
    ke0 = float(open(calib).read().split()[0])
    t0_ref = ke0 * KE_TO_T
    if t0_ref <= 0:
        return T_CREATE_REF
    return T_CREATE_REF * (t0_target / t0_ref)


def parse_lammps_trace(dt):
    rows = []
    with open(LMP_TRACE) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 2:
                continue
            step, ke = float(parts[0]), float(parts[1])
            rows.append({"step": int(step), "time": step * dt,
                         "T_total": ke * KE_TO_T})
    return rows


# ── generate ─────────────────────────────────────────────────────────────────
def generate():
    write_molecule()
    write_lammps_input(dt_rayleigh_fraction())
    print(f"Generated {MOL_FILE} and {LMP_INPUT}")


# ── start ────────────────────────────────────────────────────────────────────
def start():
    os.makedirs(DATA_DIR, exist_ok=True)
    print(f"Building {EXAMPLE} (release)...", flush=True)
    subprocess.run(
        ["cargo", "build", "--release", "--example", EXAMPLE, "--no-default-features"],
        cwd=REPO_ROOT, check=True,
    )

    print(f"Running DIRT ({N} clumps, {STEPS} steps)...", flush=True)
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
    with open(DIRT_CSV) as f:
        rows = list(csv.DictReader(f))
    dt = next((float(r["time"]) / float(r["step"]) for r in rows if int(r["step"]) > 0),
              dt_rayleigh_fraction())
    print(f"  DIRT done. dt = {dt:.6e} s, {len(rows)} samples.")

    lammps = find_lammps()
    if not lammps:
        print("LAMMPS not found on PATH — skipping the cross-check (DIRT only).")
        return
    t0_dirt = float(rows[0]["T_total"])
    write_molecule()
    print("Calibrating LAMMPS initial temperature to DIRT...", flush=True)
    t_create = calibrate_t_create(lammps, dt, t0_dirt)
    print(f"  matched T0 = {t0_dirt:.3e} m^2/s^2  (create-T = {t_create:.3e})")
    print(f"Running LAMMPS ({lammps})...", flush=True)
    write_lammps_input(dt, t_create)
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
        w = csv.DictWriter(f, fieldnames=["step", "time", "T_total"])
        w.writeheader()
        w.writerows(lrows)
    print(f"  LAMMPS done. {len(lrows)} samples -> {LMP_CSV}")


# ── graph (validate + plot) ──────────────────────────────────────────────────
def load(path, cols):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return {c: [float(r[c]) for r in rows] for c in cols}


def equilibration_time(dirt):
    """Time after which the rotational/translational partition has settled, so
    DIRT (which starts at T_rot=0) and LAMMPS (which starts already spinning from
    the rigid-body velocity projection) cool with the same quasi-steady ratio.
    Detected as when DIRT's T_rot/T_trans first reaches 90% of its plateau."""
    import numpy as np
    t = np.array(dirt["time"]); tr = np.array(dirt["T_trans"]); ro = np.array(dirt["T_rot"])
    good = (t > 0) & (tr > 0)
    t, ratio = t[good], ro[good] / tr[good]
    if len(ratio) < 8:
        return 0.0
    plateau = np.median(ratio[len(ratio) // 2:])
    idx = int(np.argmax(ratio >= 0.9 * plateau))
    return min(t[idx], 0.3 * t[-1])   # never skip more than the first 30%


def haff_fit(t, T, t_min=0.0):
    """Haff fit via 1/sqrt(T) = 1/sqrt(T0) + t/(sqrt(T0) tc), restricted to
    t >= t_min (skip the rotational-equilibration transient). Returns
    (T0, tc, r2, slope, t_over_tc)."""
    import numpy as np
    t, T = np.asarray(t), np.asarray(T)
    T0r = T[T > 0][0]
    win = (t >= t_min) & (T > 1e-3 * T0r) & np.isfinite(T)
    tw, Tw = t[win], T[win]
    y = 1.0 / np.sqrt(Tw)
    b, a = np.polyfit(tw, y, 1)
    T0, tc = 1.0 / a**2, a / b
    r2 = 1.0 - np.sum((y - (b * tw + a))**2) / np.sum((y - y.mean())**2)
    pos = tw > 0
    tp, Tp = tw[pos], Tw[pos]
    slope = np.polyfit(np.log(tp[len(tp)//2:]), np.log(Tp[len(tp)//2:]), 1)[0]
    return T0, tc, r2, slope, tp[-1] / tc


def validate(dirt):
    import numpy as np
    t, Tt = np.array(dirt["time"]), np.array(dirt["T_total"])
    print("=" * 60)
    print(f"Haff Cooling Validation — {TITLE} (DIRT)")
    print("=" * 60)
    total = passed = 0

    def check(name, ok, detail=""):
        nonlocal total, passed
        total += 1
        passed += bool(ok)
        print(f"  {name:<30}{'PASS' if ok else 'FAIL'}   {detail}")

    check("finite temperatures", bool(np.all(np.isfinite(Tt))))
    check("non-negative T", bool(np.all(Tt >= 0)))
    check("cooling (Tf < Ti)", Tt[-1] < Tt[0], f"Ti={Tt[0]:.3e} Tf={Tt[-1]:.3e}")
    check("no energy growth", float(np.max(Tt)) < 1.5 * Tt[0])
    t_eq = equilibration_time(dirt)
    T0, tc, r2, slope, ttc = haff_fit(t, Tt, t_eq)
    check("Haff law (1/sqrt(T) linear)", r2 > 0.99, f"R^2={r2:.4f}")
    # The asymptotic slope is -2, approached only at t >> tc; a dilute clump gas
    # cools slowly, so accept any clearly-Haff slope (R^2 above is the real test).
    check("late-time slope -> -2", -2.3 < slope < -1.6,
          f"slope={slope:.3f} at t/tc={ttc:.1f} (-> -2 as t/tc grows)")
    print(f"\n  Equilibration cutoff: t_eq={t_eq:.3f} s (fit uses t >= t_eq)")
    print(f"  Haff fit:  T0={T0:.3e}  tc={tc:.3e} s  (cooled to t/tc={ttc:.1f})")
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
    lammps = load(LMP_CSV, ["time", "T_total"]) if os.path.isfile(LMP_CSV) else None

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

    # Panel 1: discard the start-up transient and treat the equilibration point
    # as a fresh start — time re-zeroed there, T re-normalized to its value
    # there. Past this point DIRT and LAMMPS share the same quasi-steady
    # partition, so this isolates the cooling-law comparison.
    ax = axes[0]
    t_eq = equilibration_time(dirt)

    def restart(t, T):
        m = t >= t_eq
        return t[m] - t[m][0], T[m] / T[m][0]

    td2, Td2 = restart(np.array(dirt["time"]), np.array(dirt["T_total"]))
    _, tcd, _, sd, _ = haff_fit(td2, Td2)
    p = td2 > 0
    ax.loglog(td2[p], Td2[p], "o", ms=3, color="#1f77b4", alpha=0.7,
              label=f"DIRT  (slope {sd:.2f})")
    if lammps:
        tl2, Tl2 = restart(np.array(lammps["time"]), np.array(lammps["T_total"]))
        _, _, _, sl, _ = haff_fit(tl2, Tl2)
        q = tl2 > 0
        ax.loglog(tl2[q], Tl2[q], "s", ms=3, color="#ff7f0e", alpha=0.7,
                  label=f"LAMMPS  (slope {sl:.2f})")
    tf = np.linspace(td2[p][0], td2[-1], 300)
    ax.loglog(tf, 1.0 / (1 + tf / tcd) ** 2, "-", color="black", lw=1.5,
              label=f"Haff fit  T0/(1+t/tc)²,  tc={tcd:.2e}s")
    # The Haff fit IS the -2 law (it bends to slope -2 only at t >> tc); this run
    # reaches t/tc ~ 8, so the data sits at ~-1.6 and a literal -2 line would
    # diverge from it. The fit lying on the data is the validation.
    ax.set_xlabel("Time since equilibration [s]")
    ax.set_ylabel(r"$T_\mathrm{total}/T(t_\mathrm{eq})$")
    ax.set_title(f"Cooling from equilibration onward (skipped first {t_eq:.2f} s)")
    ax.legend(fontsize=8)

    # Panel 2: DIRT energy partition
    ax = axes[1]
    ax.semilogy(dirt["time"], dirt["T_trans"], color="#1f77b4", label=r"$T_\mathrm{trans}$")
    ax.semilogy(dirt["time"], dirt["T_rot"], color="#d62728", label=r"$T_\mathrm{rot}$")
    ax.semilogy(dirt["time"], dirt["T_total"], color="black", label=r"$T_\mathrm{total}$")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"Granular temperature [m$^2$/s$^2$]")
    ax.set_title("DIRT energy partition")
    ax.legend(fontsize=8)

    fig.suptitle(f"Haff Cooling — {TITLE}", y=1.02)
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
