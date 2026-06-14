#!/usr/bin/env python3
"""
Rolling-resistance decay benchmark driver.

A single sphere is launched in pure rolling (v = omega*R) on a flat frictional
floor [[wall]] (z = 0, normal +z) under gravity. With rolling-friction
coefficient mu_r the wall contact applies a decelerating rolling-resistance
couple tau_r = mu_r*F_n*r_eff that opposes the spin while Mindlin static
(sliding) friction enforces the rolling constraint, so the sphere decelerates at
a CONSTANT rate. For DIRT's constant-torque rolling model (couple + I = 2/5 m
R^2), and because a flat wall has r_eff = R (no curvature correction), the rate
is EXACT:

    a = (5/7) * mu_r * g

The 5/7 comes from the rolling constraint: friction supplies both the
translational deceleration and the spin-down torque, so the inertia enters as
I/R + mR = (7/5) mR. (dirt_wall now carries the full friction trio — normal +
Mindlin sliding + rolling resistance — on every wall type, so the floor is a
real wall plane; no giant frozen sphere is needed.)

This script sweeps a few mu_r, fits the v(t) deceleration slope, and validates
it against a_pred; PASS requires the fit to match within tolerance AND the
sphere to stay in pure rolling (vx ~= omega*R) throughout.

Commands (from anywhere):
    python3 examples/bench_rolling_decay/sweep.py generate   # write per-case configs
    python3 examples/bench_rolling_decay/sweep.py start      # build + run all sims -> CSV
    python3 examples/bench_rolling_decay/sweep.py graph      # validate + plot
    python3 examples/bench_rolling_decay/sweep.py            # all three, in order

If a LAMMPS binary (lmp_serial / lmp / lmp_mpi / lammps) is on PATH, each case is
also run in LAMMPS's granular model on a flat 'wall/gran' floor with a stiff
'rolling sds' contact (a close approximation of DIRT's constant-torque couple in
the saturated regime) and overlaid on the plots. LAMMPS is optional — without
it, only DIRT runs and the benchmark still validates against theory.

Outputs:
    sweep/<case>/config.toml    DIRT configs                          (gitignored)
    sweep/<case>/in.lammps      LAMMPS inputs                         (gitignored)
    data/decay_<mu_r>.csv       per-case DIRT time series             (gitignored)
    data/sweep.csv              fitted-slope summary (DIRT)           (gitignored)
    data/sweep_lammps.csv       fitted-slope summary (LAMMPS)         (gitignored)
    plots/*.png                 final figures                         (tracked)

Reference: the analytical pure-rolling deceleration of a sphere under a
rolling-resistance couple — e.g. J. Ai, J.-F. Chen, J.M. Rotter, J.Y. Ooi,
"Assessment of rolling resistance models in discrete element simulations",
Powder Technology 206 (2011) 269-282 (model A, "constant directional torque").
"""

import os
import sys
import csv
import math
import shutil
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_rolling_decay"

SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
SWEEP_CSV = os.path.join(DATA_DIR, "sweep.csv")            # DIRT fitted slopes
LAMMPS_CSV = os.path.join(DATA_DIR, "sweep_lammps.csv")    # LAMMPS fitted slopes

# LAMMPS binary candidates, in preference order. LAMMPS is optional.
LAMMPS_BINS = ["lmp_serial", "lmp", "lmp_mpi", "lammps"]

# ── Material / geometry ───────────────────────────────────────────────────────
NU = 0.3              # Poisson ratio
MU = 0.5              # sliding friction — healthy, so static friction holds rolling
E_N = 0.5             # normal restitution (only damps the tiny settling transient)
YOUNGS_MOD = 1.0e8    # Pa — soft (keeps dt large; the decay rate is E-independent)
DENSITY = 2500.0      # kg/m^3
RADIUS = 0.005        # m — sphere radius R
GRAVITY = 9.81        # m/s^2
V0 = 0.03             # m/s — initial rolling speed
DT = 1.0e-5           # s
STEPS = 40000         # plenty for the slowest (smallest mu_r) case to stop

# Swept rolling-friction coefficients (modest, keeps every case < few seconds).
MU_R_LIST = [0.02, 0.05, 0.10]

# ── Theory ────────────────────────────────────────────────────────────────────
def a_pred(mu_r):
    """Exact constant pure-rolling deceleration for DIRT's constant-torque model
    on a flat wall (r_eff = R, no curvature correction)."""
    return (5.0 / 7.0) * mu_r * GRAVITY


# ── DIRT config template ──────────────────────────────────────────────────────
TOML_TEMPLATE = """[comm]
processors_x = 1
processors_y = 1
processors_z = 1
[domain]
x_low = -0.05
x_high = 0.10
y_low = -0.02
y_high = 0.02
z_low = 0.0
z_high = 0.05
boundary_x = "fixed"
boundary_y = "fixed"
boundary_z = "fixed"
[neighbor]
skin_fraction = 1.1
bin_size = 0.02
every = 1
[gravity]
gz = -{g}
[dem]
contact_model = "hertz"
rolling_model = "constant"
[[dem.materials]]
name = "grain"
youngs_mod = {youngs:.6e}
poisson_ratio = {nu}
restitution = {e_n}
friction = {mu}
rolling_friction = {mu_r}
[[particles.insert]]
material = "grain"
count = 1
radius = {radius}
density = {density}
velocity_x = {v0}
region = {{ type = "block", min = [-0.005, -0.005, 0.0055], max = [0.005, 0.005, 0.0065] }}
[[wall]]
point_x = 0.0
point_y = 0.0
point_z = 0.0
normal_x = 0.0
normal_y = 0.0
normal_z = 1.0
material = "grain"
[output]
dir = "{outdir}"
[run]
steps = {steps}
thermo = {steps}
dt = {dt:.6e}
"""

# ── LAMMPS template ───────────────────────────────────────────────────────────
# Mirrors the DIRT setup: a single sphere rolling in pure rolling on a flat
# granular floor wall (z = 0) under gravity. LAMMPS's 'rolling sds'
# (spring-dashpot-slider) with a stiff spring saturates at tau = mu_r*F_n*r_eff,
# the same cap DIRT's constant couple applies — and for a flat wall r_eff = R in
# both codes — so in steady rolling the two decelerations agree. We give the
# sphere its velocity AND the matching spin so it starts in pure rolling.
LMP_TEMPLATE = """units           si
atom_style      sphere
boundary        f f f
comm_modify     vel yes
region          box block -0.1 0.2 -0.05 0.05 0.0 0.1 units box
create_box      1 box
create_atoms    1 single 0.0 0.0 {ztop:.9f} units box          # rolling sphere
group           mover id 1
set             group mover diameter {dmover} density {density}
pair_style      granular
pair_coeff      1 1 hertz/material {youngs:.6e} {e_n} {nu} &
                tangential mindlin NULL 1.0 {mu} &
                rolling sds {kroll:.6e} {droll:.6e} {mu_r} &
                twisting none
# Flat granular floor wall at z = 0, same contact law as the pair (rolling on).
fix             floor all wall/gran granular &
                hertz/material {youngs:.6e} {e_n} {nu} &
                tangential mindlin NULL 1.0 {mu} &
                rolling sds {kroll:.6e} {droll:.6e} {mu_r} &
                twisting none &
                zplane 0.0 NULL
# Pure rolling start: translational v0 plus matching spin omega = v0/R about +y.
velocity        mover set {v0} 0.0 0.0 units box
set             group mover omega 0.0 {omega:.9f} 0.0
fix             grav all gravity {g} vector 0 0 -1
fix             integ mover nve/sphere
timestep        {dt:.6e}
thermo          {steps}
compute         vxm mover reduce sum vx
fix             rec mover ave/time 1 1 {every} c_vxm file {out} mode scalar
run             {steps}
"""


# ── helpers ───────────────────────────────────────────────────────────────────
def case_tag(mu_r):
    return f"mu_r_{mu_r:g}"


def case_dir(mu_r):
    return os.path.join(SWEEP_DIR, case_tag(mu_r))


def decay_csv(mu_r):
    return os.path.join(DATA_DIR, f"decay_{mu_r:g}.csv")


def find_lammps():
    for b in LAMMPS_BINS:
        path = shutil.which(b)
        if path:
            return path
    return None


def _dirt_config(mu_r, outdir):
    return TOML_TEMPLATE.format(
        g=GRAVITY, youngs=YOUNGS_MOD, nu=NU, e_n=E_N, mu=MU, mu_r=mu_r,
        radius=RADIUS, density=DENSITY, v0=V0,
        outdir=outdir, steps=STEPS, dt=DT,
    )


# ── generate ──────────────────────────────────────────────────────────────────
def generate():
    n = 0
    for mu_r in MU_R_LIST:
        cdir = case_dir(mu_r)
        os.makedirs(cdir, exist_ok=True)
        with open(os.path.join(cdir, "config.toml"), "w") as f:
            f.write(_dirt_config(mu_r, cdir))
        n += 1
    print(f"Generated {n} DIRT sweep configs under {SWEEP_DIR}")


# ── start ─────────────────────────────────────────────────────────────────────
SWEEP_FIELDS = ["mu_r", "a_fit", "a_pred", "rel_err", "max_slip", "npts"]


def _read_timeseries(path):
    """Read t,x,vx,omega from a DIRT rolling_decay_results.csv (skip '#' lines)."""
    rows = []
    with open(path) as f:
        rdr = csv.reader(f)
        header_seen = False
        for parts in rdr:
            if not parts or parts[0].startswith("#"):
                continue
            if not header_seen and parts[0] == "t":
                header_seen = True
                continue
            if len(parts) >= 4:
                rows.append(tuple(float(p) for p in parts[:4]))
    return rows  # list of (t, x, vx, omega)


def _fit_decay(rows):
    """Linear-fit vx(t) over the moving portion; return (a_fit, max_slip, npts).

    a_fit = -slope (deceleration, positive). max_slip = max|vx - omega*R| / V0
    over the fitted window (rolling-purity check)."""
    # Use only the initial monotone decay, before the sphere first reaches near
    # zero. (LAMMPS's SDS rolling spring can store energy and make vx oscillate
    # slightly about zero after stopping; that tail must not enter the fit.)
    lo, hi = 0.05 * V0, 0.95 * V0
    prefix = []
    for r in rows:
        if r[2] < lo:
            break
        prefix.append(r)
    win = [r for r in prefix if lo <= r[2] <= hi]
    if len(win) < 10:
        win = [r for r in prefix if r[2] > 1e-4] or prefix
    n = len(win)
    if n < 2:
        return None
    ts = [r[0] for r in win]
    vs = [r[2] for r in win]
    tbar = sum(ts) / n
    vbar = sum(vs) / n
    sxx = sum((t - tbar) ** 2 for t in ts)
    sxy = sum((t - tbar) * (v - vbar) for t, v in zip(ts, vs))
    slope = sxy / sxx if sxx > 0 else 0.0
    a_fit = -slope
    max_slip = max(abs(r[2] - r[3] * RADIUS) for r in win) / V0
    return a_fit, max_slip, n


def _run_dirt(mu_r):
    cdir = case_dir(mu_r)
    config = os.path.join(cdir, "config.toml")
    res = os.path.join(cdir, "data", "rolling_decay_results.csv")
    if os.path.exists(res):
        os.remove(res)
    log = os.path.join(cdir, "run.log")
    with open(log, "w") as lf:
        proc = subprocess.run(
            ["cargo", "run", "--release", "--example", EXAMPLE,
             "--no-default-features", "--", config],
            cwd=REPO_ROOT, stdout=lf, stderr=subprocess.STDOUT,
        )
    if proc.returncode != 0 or not os.path.isfile(res):
        return None
    rows = _read_timeseries(res)
    # Persist the time series for plotting.
    with open(decay_csv(mu_r), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t", "x", "vx", "omega"])
        w.writerows(rows)
    return rows


def _lammps_config(mu_r, cdir):
    # Stiff rolling spring: pick k_roll so the spring saturates almost
    # immediately (overshoots the cap within a step or two), then the slider
    # holds tau = mu_r*F_n*r_eff. Critically damped is unnecessary; use light
    # damping. Units: k_roll in N*m/rad ~ a few * (mu_r*m*g*r_eff)/(small angle).
    kroll = 1.0e2
    droll = 0.0
    out = os.path.join(cdir, "vx.txt")
    in_path = os.path.join(cdir, "in.lammps")
    with open(in_path, "w") as f:
        f.write(LMP_TEMPLATE.format(
            ztop=RADIUS - 2.0e-6,
            dmover=2.0 * RADIUS, density=DENSITY,
            youngs=YOUNGS_MOD, e_n=E_N, nu=NU, mu=MU,
            kroll=kroll, droll=droll, mu_r=mu_r,
            v0=V0, omega=V0 / RADIUS, g=GRAVITY, dt=DT,
            steps=STEPS, every=max(1, STEPS // 2000), out=out,
        ))
    return in_path, out


def _run_lammps(lammps, mu_r):
    cdir = case_dir(mu_r)
    in_path, out = _lammps_config(mu_r, cdir)
    if os.path.exists(out):
        os.remove(out)
    proc = subprocess.run([lammps, "-in", in_path], cwd=cdir,
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if proc.returncode != 0 or not os.path.isfile(out):
        return None
    rows = []
    with open(out) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split()
            if len(p) >= 2:
                # ave/time scalar file: col0 = timestep, col1 = vx
                rows.append((float(p[0]) * DT, 0.0, float(p[1]), 0.0))
    fit = _fit_decay(rows)
    return fit


def _write_csv(path, fields, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})


def start():
    os.makedirs(DATA_DIR, exist_ok=True)
    print(f"Building {EXAMPLE} (release)...", flush=True)
    subprocess.run(["cargo", "build", "--release", "--example", EXAMPLE,
                    "--no-default-features"], cwd=REPO_ROOT, check=True)

    lammps = find_lammps()
    print(f"LAMMPS: {lammps}" if lammps else
          "LAMMPS: not found on PATH — running DIRT only.")

    # Wipe stale summaries so old results can never be re-plotted.
    for stale in (SWEEP_CSV, LAMMPS_CSV):
        if os.path.exists(stale):
            os.remove(stale)

    dirt_rows, lmp_rows = [], []
    n = len(MU_R_LIST)
    for i, mu_r in enumerate(MU_R_LIST, 1):
        cdir = case_dir(mu_r)
        if not os.path.isfile(os.path.join(cdir, "config.toml")):
            print(f"  [{i}/{n}] missing config for mu_r={mu_r} — run 'generate' first.")
            continue
        print(f"  [{i}/{n}] mu_r={mu_r:<5}", end="  ", flush=True)
        rows = _run_dirt(mu_r)
        if rows is None:
            print("DIRT FAILED")
            continue
        fit = _fit_decay(rows)
        if fit is None:
            print("DIRT (no fit)")
            continue
        a_fit, max_slip, npts = fit
        ap = a_pred(mu_r)
        rel = abs(a_fit - ap) / ap
        dirt_rows.append({"mu_r": mu_r, "a_fit": a_fit, "a_pred": ap,
                          "rel_err": rel, "max_slip": max_slip, "npts": npts})
        print(f"DIRT a_fit={a_fit:.4f} a_pred={ap:.4f} rel={rel*100:.1f}% "
              f"slip={max_slip:.2e}", end="")
        if lammps:
            lf = _run_lammps(lammps, mu_r)
            if lf:
                la, lslip, lnp = lf
                lmp_rows.append({"mu_r": mu_r, "a_fit": la, "a_pred": ap,
                                 "rel_err": abs(la - ap) / ap,
                                 "max_slip": lslip, "npts": lnp})
                print(f"   LAMMPS a_fit={la:.4f}", end="")
        print()

    if not dirt_rows:
        print("\nERROR: no DIRT results collected.")
        sys.exit(1)
    _write_csv(SWEEP_CSV, SWEEP_FIELDS, dirt_rows)
    print(f"\nDIRT:   {len(dirt_rows)}/{n} cases -> {SWEEP_CSV}")
    if lmp_rows:
        _write_csv(LAMMPS_CSV, SWEEP_FIELDS, lmp_rows)
        print(f"LAMMPS: {len(lmp_rows)}/{n} cases -> {LAMMPS_CSV}")


# ── graph (validate + plot) ───────────────────────────────────────────────────
SLOPE_TOL = 0.02    # 2% relative error on the fitted deceleration (theory is exact)
SLIP_TOL = 1.0e-2   # |vx - omega*R| must stay < 1% of V0 (pure rolling)


def _load(path):
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]


def validate(rows):
    print("\n=== Rolling-resistance decay validation ===")
    print(f"  R={RADIUS} m  flat wall (r_eff = R)  mu={MU}  g={GRAVITY}")
    print(f"  model: a = (5/7) mu_r g   (exact)")
    print(f"  {'mu_r':>6}{'a_fit':>10}{'a_pred':>10}{'rel_err':>9}"
          f"{'max_slip':>11}  note")
    ok = True
    for r in sorted(rows, key=lambda x: x["mu_r"]):
        note = ""
        if r["rel_err"] > SLOPE_TOL:
            note = "SLOPE MISMATCH"; ok = False
        if r["max_slip"] > SLIP_TOL:
            note = (note + " NOT-PURE-ROLLING").strip(); ok = False
        print(f"  {r['mu_r']:>6.3f}{r['a_fit']:>10.4f}{r['a_pred']:>10.4f}"
              f"{r['rel_err']*100:>8.1f}%{r['max_slip']:>11.2e}  {note}")
    print(f"\n  tolerances: slope <= {SLOPE_TOL*100:.0f}% rel, "
          f"slip <= {SLIP_TOL*100:.0f}% of V0")
    print("RESULT:", "PASS" if ok else "FAIL")
    return ok


def compare_codes(dirt, lammps):
    lmp = {round(r["mu_r"], 4): r for r in lammps}
    print("\n=== DIRT vs LAMMPS (fitted deceleration) ===")
    print(f"  {'mu_r':>6}{'DIRT':>10}{'LAMMPS':>10}{'d_a':>10}")
    for r in sorted(dirt, key=lambda x: x["mu_r"]):
        l = lmp.get(round(r["mu_r"], 4))
        if not l:
            continue
        d = r["a_fit"] - l["a_fit"]
        print(f"  {r['mu_r']:>6.3f}{r['a_fit']:>10.4f}{l['a_fit']:>10.4f}{d:>+10.4f}")


def plot(dirt, lammps):
    os.makedirs(PLOT_DIR, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150, "font.size": 11})

    # ── v(t) decay curves: measured vs predicted, per mu_r ──
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    colors = plt.cm.viridis([0.15, 0.5, 0.82])
    for c, mu_r in zip(colors, MU_R_LIST):
        path = decay_csv(mu_r)
        if not os.path.isfile(path):
            continue
        with open(path) as f:
            ts, vs = [], []
            for row in csv.DictReader(f):
                ts.append(float(row["t"]))
                vs.append(float(row["vx"]))
        ax.plot(ts, vs, "-", color=c, lw=1.6,
                label=fr"DIRT $\mu_r$={mu_r}")
        # Predicted line: v(t) = V0 - a_pred*t, clamped at 0.
        ap = a_pred(mu_r)
        tstop = V0 / ap
        ax.plot([0, tstop], [V0, 0.0], "k--", lw=1.0,
                label="theory" if mu_r == MU_R_LIST[0] else None)
    ax.set_xlabel("time t (s)")
    ax.set_ylabel("translational velocity $v_x$ (m/s)")
    ax.set_title(r"Pure-rolling decay: measured vs $a=(5/7)\mu_r g$")
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "velocity_decay.png"))
    plt.close(fig)

    # ── fitted deceleration vs theory (and LAMMPS overlay) ──
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    d = sorted(dirt, key=lambda x: x["mu_r"])
    mu = [r["mu_r"] for r in d]
    ax.plot(mu, [r["a_pred"] for r in d], "k-", label="theory (5/7)")
    ax.plot(mu, [r["a_fit"] for r in d], "o", ms=8, color="tab:blue",
            label="DIRT (fit)")
    if lammps:
        l = sorted(lammps, key=lambda x: x["mu_r"])
        ax.plot([r["mu_r"] for r in l], [r["a_fit"] for r in l], "s", ms=8,
                mfc="none", color="tab:red", label="LAMMPS (fit)")
    ax.set_xlabel(r"rolling friction $\mu_r$")
    ax.set_ylabel(r"deceleration $a$ (m/s$^2$)")
    ax.set_title("Fitted pure-rolling deceleration vs theory")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "deceleration_vs_mu_r.png"))
    plt.close(fig)
    print(f"\nFigures -> {PLOT_DIR}/velocity_decay.png, deceleration_vs_mu_r.png")


def graph():
    dirt = _load(SWEEP_CSV)
    if not dirt:
        print(f"No {SWEEP_CSV} — run 'start' first.")
        return False
    lammps = _load(LAMMPS_CSV)
    ok = validate(dirt)
    if lammps:
        compare_codes(dirt, lammps)
    else:
        print("\n(no LAMMPS sweep — plotting DIRT only)")
    plot(dirt, lammps)
    return ok


# ── dispatch ──────────────────────────────────────────────────────────────────
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
