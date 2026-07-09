#!/usr/bin/env python3
"""
Cundall non-viscous global damping benchmark driver.

Validates the [[cundall]] fix (dirt_fixes) — DIRT's implementation of LAMMPS
`fix damping/cundall` — against its EXACT analytical effect on a single free
particle, for several damping coefficients.

A lone sphere (no walls, no contacts) is launched straight up (+z) under gravity
and given an initial spin about +z with a constant applied torque. The Cundall
fix damps force and torque component-by-component keyed to the sign of the
mechanical power (LAMMPS `fix damping/cundall`, Cundall 1987):

    F_k <- F_k * (1 - gamma_l * sign(F_k * v_k))
    T_k <- T_k * (1 - gamma_a * sign(T_k * omega_k))

Because the sign flips as the sphere passes its apex (rising -> falling) and as
the spin passes zero (spin-down -> spin-up), each motion is piecewise
constant-acceleration with the EXACT, mass-independent rates

    a_up       = -g (1 + gamma_l)          (v_z > 0, gravity opposes motion)
    a_down     = -g (1 - gamma_l)          (v_z < 0, gravity along motion)
    alpha_down =  T_z (1 + gamma_a) / I    (omega_z > 0, torque opposes spin)
    alpha_up   =  T_z (1 - gamma_a) / I    (omega_z < 0, torque along spin)

The script sweeps a few gamma, linear-fits each of the four phases, and checks
them against theory; PASS requires ALL four rates (for every gamma) within
tolerance.

If a LAMMPS binary (lmp_serial / lmp / lmp_mpi / lammps) is on PATH, the LINEAR
part of each case is ALSO run in LAMMPS with the native `fix damping/cundall`
(single sphere, `fix gravity` before the damping fix so gravity is damped) and
its fitted a_up/a_down are overlaid — a true independent cross-code check. The
angular part is validated against theory only (LAMMPS has no built-in constant
body torque). LAMMPS is optional; without it the theory check still gates.

Commands (from anywhere):
    python3 examples/bench_cundall_damping/sweep.py generate
    python3 examples/bench_cundall_damping/sweep.py start
    python3 examples/bench_cundall_damping/sweep.py graph
    python3 examples/bench_cundall_damping/sweep.py            # all three

Outputs:
    sweep/<case>/config.toml     DIRT configs                        (gitignored)
    sweep/<case>/in.lammps       LAMMPS inputs                       (gitignored)
    data/case_<gamma>.csv        per-case DIRT time series           (gitignored)
    data/sweep.csv               fitted-rate summary (DIRT)          (gitignored)
    data/sweep_lammps.csv        fitted-rate summary (LAMMPS linear) (gitignored)
    plots/*.png                  figures                             (tracked)

Reference: LAMMPS `fix damping/cundall` (doc/src/fix_damping_cundall.rst;
src/GRANULAR/fix_damping_cundall.cpp) after P.A. Cundall, "Distinct element
models of rock and soil structure", 1987; as in Yade-DEM and PFC.
"""

import os
import sys
import csv
import shutil
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_cundall_damping"

SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
SWEEP_CSV = os.path.join(DATA_DIR, "sweep.csv")
LAMMPS_CSV = os.path.join(DATA_DIR, "sweep_lammps.csv")

LAMMPS_BINS = ["lmp_serial", "lmp", "lmp_mpi", "lammps"]

# ── Fixed scenario parameters (must match main.rs constants) ──────────────────
GRAVITY = 9.81      # m/s^2 (config [gravity] gz = -GRAVITY)
V0 = 3.0            # m/s   initial +z launch speed (config velocity_z)
OMEGA0 = 5.0        # rad/s initial +z spin (main.rs OMEGA0_Z)
TORQUE_Z = -1.0e-6  # N*m   constant applied torque about z (main.rs TORQUE_Z)
RADIUS = 0.005      # m
DENSITY = 2500.0    # kg/m^3
DT = 1.0e-4         # s
STEPS = 6000        # 0.6 s — both v_z and omega_z cross zero well within this

# Swept damping coefficients (gamma_l = gamma_a = gamma for each case).
GAMMA_LIST = [0.2, 0.5, 0.8]

# ── DIRT config template ──────────────────────────────────────────────────────
TOML_TEMPLATE = """[comm]
processors_x = 1
processors_y = 1
processors_z = 1
[domain]
x_low = -0.05
x_high = 0.05
y_low = -0.05
y_high = 0.05
z_low = -0.6
z_high = 0.6
boundary_x = "fixed"
boundary_y = "fixed"
boundary_z = "fixed"
[neighbor]
skin_fraction = 1.1
bin_size = 0.05
every = 1
[gravity]
gz = -{g}
[dem]
contact_model = "hertz"
[[dem.materials]]
name = "grain"
youngs_mod = 1.0e8
poisson_ratio = 0.3
restitution = 0.5
friction = 0.5
[[particles.insert]]
material = "grain"
count = 1
radius = {radius}
density = {density}
velocity_z = {v0}
region = {{ type = "block", min = [-0.002, -0.002, -0.002], max = [0.002, 0.002, 0.002] }}
[[cundall]]
group = "all"
gamma_l = {gamma}
gamma_a = {gamma}
[output]
dir = "{outdir}"
[run]
steps = {steps}
thermo = {steps}
dt = {dt:.6e}
"""

# ── LAMMPS template (LINEAR cross-check only) ─────────────────────────────────
# Single sphere launched +z under gravity, with the native fix damping/cundall.
# fix gravity is declared BEFORE the damping fix so its force IS damped (per the
# LAMMPS doc note). No pair interactions (lone particle). No spin/torque here —
# LAMMPS validates the linear branch; the angular branch is checked vs theory.
LMP_TEMPLATE = """units           si
atom_style      sphere
boundary        f f f
comm_modify     vel yes
region          box block -0.05 0.05 -0.05 0.05 -0.6 0.6 units box
create_box      1 box
create_atoms    1 single 0.0 0.0 0.0 units box
set             group all diameter {d} density {density}
pair_style      none
velocity        all set 0.0 0.0 {v0} units box
fix             grav all gravity {g} vector 0 0 -1
fix             damp all damping/cundall {gl} {ga}
fix             integ all nve/sphere
timestep        {dt:.6e}
thermo          {steps}
compute         vzc all reduce sum vz
fix             rec all ave/time 1 1 {every} c_vzc file {out} mode scalar
run             {steps}
"""


def case_tag(gamma):
    return f"gamma_{gamma:g}"


def case_dir(gamma):
    return os.path.join(SWEEP_DIR, case_tag(gamma))


def series_csv(gamma):
    return os.path.join(DATA_DIR, f"case_{gamma:g}.csv")


def find_lammps():
    for b in LAMMPS_BINS:
        p = shutil.which(b)
        if p:
            return p
    return None


# ── generate ──────────────────────────────────────────────────────────────────
def generate():
    n = 0
    for gamma in GAMMA_LIST:
        cdir = case_dir(gamma)
        os.makedirs(cdir, exist_ok=True)
        with open(os.path.join(cdir, "config.toml"), "w") as f:
            f.write(TOML_TEMPLATE.format(
                g=GRAVITY, radius=RADIUS, density=DENSITY, v0=V0,
                gamma=gamma, outdir=cdir, steps=STEPS, dt=DT,
            ))
        n += 1
    print(f"Generated {n} DIRT sweep configs under {SWEEP_DIR}")


# ── fitting ───────────────────────────────────────────────────────────────────
def _read_header_and_rows(path):
    """Return (params dict, rows[(t,vz,omega_z)]) from a cundall_results.csv."""
    params, rows = {}, []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                for kv in s.lstrip("#").split():
                    if "=" in kv:
                        k, v = kv.split("=", 1)
                        params[k] = float(v)
                continue
            if s.startswith("t,"):
                continue
            p = s.split(",")
            if len(p) >= 3:
                rows.append((float(p[0]), float(p[1]), float(p[2])))
    return params, rows


def _fit_slope(ts, ys):
    """Least-squares slope dy/dt over the given points (None if < 3)."""
    n = len(ts)
    if n < 3:
        return None
    tb = sum(ts) / n
    yb = sum(ys) / n
    sxx = sum((t - tb) ** 2 for t in ts)
    sxy = sum((t - tb) * (y - yb) for t, y in zip(ts, ys))
    return sxy / sxx if sxx > 0 else None


def _fit_phase(rows, col, lo, hi):
    """Fit slope of column `col` (1=vz, 2=omega_z) over samples whose value lies
    in [lo, hi]. Excludes the sign-flip transition and the launch/settle ends."""
    ts = [r[0] for r in rows if lo <= r[col] <= hi]
    ys = [r[col] for r in rows if lo <= r[col] <= hi]
    return _fit_slope(ts, ys)


def _fit_all(rows):
    """Fit the four phase rates from a DIRT time series."""
    # Linear: rising = vz in [+0.1 V0, +0.9 V0]; falling = vz in [-0.9 V0, -0.1 V0].
    a_up = _fit_phase(rows, 1, 0.1 * V0, 0.9 * V0)
    a_down = _fit_phase(rows, 1, -0.9 * V0, -0.1 * V0)
    # Angular: spin-down = omega in [+0.1 W0, +0.9 W0]; the spin keeps growing
    # negative (torque never reverses), so spin-up window is [-N, -0.2 W0].
    al_down = _fit_phase(rows, 2, 0.1 * OMEGA0, 0.9 * OMEGA0)
    al_up = _fit_phase(rows, 2, -1.0e30, -0.2 * OMEGA0)
    return a_up, a_down, al_down, al_up


def _predict(params):
    g = params["g"]
    gl = params["gamma_l"]
    ga = params["gamma_a"]
    tz = params["torque_z"]
    invI = params["inv_inertia"]
    return {
        "a_up": -g * (1.0 + gl),
        "a_down": -g * (1.0 - gl),
        "alpha_down": tz * (1.0 + ga) * invI,
        "alpha_up": tz * (1.0 - ga) * invI,
    }


# ── start ─────────────────────────────────────────────────────────────────────
SWEEP_FIELDS = [
    "gamma", "a_up", "a_up_pred", "a_down", "a_down_pred",
    "alpha_down", "alpha_down_pred", "alpha_up", "alpha_up_pred",
]
LMP_FIELDS = ["gamma", "a_up", "a_up_pred", "a_down", "a_down_pred"]


def _run_dirt(gamma):
    cdir = case_dir(gamma)
    config = os.path.join(cdir, "config.toml")
    res = os.path.join(cdir, "data", "cundall_results.csv")
    if os.path.exists(res):
        os.remove(res)
    log = os.path.join(cdir, "run.log")
    with open(log, "w") as lf:
        proc = subprocess.run(
            ["cargo", "run", "--release", "--example", EXAMPLE,
             "--no-default-features", "--features", "precision-double", "--", config],
            cwd=REPO_ROOT, stdout=lf, stderr=subprocess.STDOUT,
        )
    if proc.returncode != 0 or not os.path.isfile(res):
        return None
    params, rows = _read_header_and_rows(res)
    with open(series_csv(gamma), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t", "vz", "omega_z"])
        w.writerows(rows)
    return params, rows


def _run_lammps(lammps, gamma):
    cdir = case_dir(gamma)
    out = os.path.join(cdir, "vz.txt")
    in_path = os.path.join(cdir, "in.lammps")
    if os.path.exists(out):
        os.remove(out)
    with open(in_path, "w") as f:
        f.write(LMP_TEMPLATE.format(
            d=2.0 * RADIUS, density=DENSITY, v0=V0, g=GRAVITY,
            gl=gamma, ga=gamma, dt=DT, steps=STEPS,
            every=max(1, STEPS // 2000), out=out,
        ))
    proc = subprocess.run([lammps, "-in", in_path], cwd=cdir,
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if proc.returncode != 0 or not os.path.isfile(out):
        return None
    rows = []
    with open(out) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            p = s.split()
            if len(p) >= 2:
                rows.append((float(p[0]) * DT, float(p[1]), 0.0))
    a_up = _fit_phase(rows, 1, 0.1 * V0, 0.9 * V0)
    a_down = _fit_phase(rows, 1, -0.9 * V0, -0.1 * V0)
    return a_up, a_down


def _write_csv(path, fields, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def start():
    os.makedirs(DATA_DIR, exist_ok=True)
    print(f"Building {EXAMPLE} (release)...", flush=True)
    subprocess.run(["cargo", "build", "--release", "--example", EXAMPLE,
                    "--no-default-features", "--features", "precision-double"],
                   cwd=REPO_ROOT, check=True)

    lammps = find_lammps()
    print(f"LAMMPS: {lammps}" if lammps else
          "LAMMPS: not found on PATH — running DIRT only (linear cross-check skipped).")

    for stale in (SWEEP_CSV, LAMMPS_CSV):
        if os.path.exists(stale):
            os.remove(stale)

    dirt_rows, lmp_rows = [], []
    n = len(GAMMA_LIST)
    for i, gamma in enumerate(GAMMA_LIST, 1):
        cdir = case_dir(gamma)
        if not os.path.isfile(os.path.join(cdir, "config.toml")):
            print(f"  [{i}/{n}] missing config for gamma={gamma} — run 'generate' first.")
            continue
        print(f"  [{i}/{n}] gamma={gamma:<4}", end="  ", flush=True)
        got = _run_dirt(gamma)
        if got is None:
            print("DIRT FAILED")
            continue
        params, rows = got
        a_up, a_down, al_down, al_up = _fit_all(rows)
        if None in (a_up, a_down, al_down, al_up):
            print("DIRT (incomplete phases)")
            continue
        pred = _predict(params)
        dirt_rows.append({
            "gamma": gamma,
            "a_up": a_up, "a_up_pred": pred["a_up"],
            "a_down": a_down, "a_down_pred": pred["a_down"],
            "alpha_down": al_down, "alpha_down_pred": pred["alpha_down"],
            "alpha_up": al_up, "alpha_up_pred": pred["alpha_up"],
        })
        print(f"DIRT a_up={a_up:.3f}/{pred['a_up']:.3f} "
              f"a_dn={a_down:.3f}/{pred['a_down']:.3f} "
              f"al_dn={al_down:.1f}/{pred['alpha_down']:.1f} "
              f"al_up={al_up:.1f}/{pred['alpha_up']:.1f}", end="")
        if lammps:
            lf = _run_lammps(lammps, gamma)
            if lf and None not in lf:
                lmp_rows.append({
                    "gamma": gamma,
                    "a_up": lf[0], "a_up_pred": pred["a_up"],
                    "a_down": lf[1], "a_down_pred": pred["a_down"],
                })
                print(f"   LAMMPS a_up={lf[0]:.3f} a_dn={lf[1]:.3f}", end="")
        print()

    if not dirt_rows:
        print("\nERROR: no DIRT results collected.")
        sys.exit(1)
    _write_csv(SWEEP_CSV, SWEEP_FIELDS, dirt_rows)
    print(f"\nDIRT:   {len(dirt_rows)}/{n} cases -> {SWEEP_CSV}")
    if lmp_rows:
        _write_csv(LAMMPS_CSV, LMP_FIELDS, lmp_rows)
        print(f"LAMMPS: {len(lmp_rows)}/{n} cases -> {LAMMPS_CSV}")


# ── graph (validate + plot) ───────────────────────────────────────────────────
REL_TOL = 0.01   # 1% relative error on every fitted rate (theory is exact)


def _load(path):
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        return [{k: (float(v) if v != "" else None) for k, v in r.items()}
                for r in csv.DictReader(f)]


def _rel(fit, pred):
    return abs(fit - pred) / abs(pred) if pred != 0 else abs(fit)


def validate(rows):
    print("\n=== Cundall non-viscous damping validation ===")
    print(f"  g={GRAVITY}  V0={V0}  Omega0={OMEGA0}  Tz={TORQUE_Z}")
    print("  exact: a_up=-g(1+gl)  a_down=-g(1-gl)  "
          "alpha_down=Tz(1+ga)/I  alpha_up=Tz(1-ga)/I")
    hdr = f"  {'gamma':>6}{'a_up':>18}{'a_down':>18}{'alpha_down':>20}{'alpha_up':>20}  note"
    print(hdr)
    ok = True
    keys = [("a_up", "a_up_pred"), ("a_down", "a_down_pred"),
            ("alpha_down", "alpha_down_pred"), ("alpha_up", "alpha_up_pred")]
    for r in sorted(rows, key=lambda x: x["gamma"]):
        bad = []
        cells = []
        for fit_k, pred_k in keys:
            rel = _rel(r[fit_k], r[pred_k])
            cells.append(f"{r[fit_k]:>9.3f}/{r[pred_k]:<7.3f}")
            if rel > REL_TOL:
                bad.append(f"{fit_k}({rel*100:.1f}%)")
        if bad:
            ok = False
        note = "" if not bad else "MISMATCH:" + ",".join(bad)
        print(f"  {r['gamma']:>6.2f}" + "".join(f"{c:>18}" for c in cells) + f"  {note}")
    print(f"\n  tolerance: <= {REL_TOL*100:.0f}% relative error on all four rates, all gamma")
    print("RESULT:", "PASS" if ok else "FAIL")
    return ok


def compare_codes(dirt, lammps):
    lmp = {round(r["gamma"], 4): r for r in lammps}
    print("\n=== DIRT vs LAMMPS fix damping/cundall (linear a_up / a_down) ===")
    print(f"  {'gamma':>6}{'DIRT a_up':>12}{'LMP a_up':>12}"
          f"{'DIRT a_dn':>12}{'LMP a_dn':>12}")
    for r in sorted(dirt, key=lambda x: x["gamma"]):
        l = lmp.get(round(r["gamma"], 4))
        if not l:
            continue
        print(f"  {r['gamma']:>6.2f}{r['a_up']:>12.3f}{l['a_up']:>12.3f}"
              f"{r['a_down']:>12.3f}{l['a_down']:>12.3f}")


def plot(dirt, lammps):
    os.makedirs(PLOT_DIR, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150, "font.size": 11})

    # ── v_z(t) and omega_z(t) traces vs piecewise-linear theory ──
    fig, (axv, axw) = plt.subplots(1, 2, figsize=(11.0, 4.4))
    colors = plt.cm.viridis([0.15, 0.5, 0.82])
    for c, gamma in zip(colors, GAMMA_LIST):
        path = series_csv(gamma)
        if not os.path.isfile(path):
            continue
        ts, vs, ws = [], [], []
        with open(path) as f:
            for row in csv.DictReader(f):
                ts.append(float(row["t"]))
                vs.append(float(row["vz"]))
                ws.append(float(row["omega_z"]))
        axv.plot(ts, vs, "-", color=c, lw=1.6, label=fr"$\gamma$={gamma}")
        axw.plot(ts, ws, "-", color=c, lw=1.6, label=fr"$\gamma$={gamma}")
    axv.axhline(0, color="k", lw=0.6)
    axw.axhline(0, color="k", lw=0.6)
    axv.set_xlabel("time t (s)"); axv.set_ylabel(r"$v_z$ (m/s)")
    axv.set_title(r"Linear: throw-up under damped gravity")
    axw.set_xlabel("time t (s)"); axw.set_ylabel(r"$\omega_z$ (rad/s)")
    axw.set_title(r"Angular: spin-down under damped torque")
    axv.legend(fontsize=9); axw.legend(fontsize=9)
    fig.suptitle("Cundall non-viscous damping: sign-flip at apex / zero-spin")
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "cundall_traces.png"))
    plt.close(fig)

    # ── fitted rates vs theory (and LAMMPS overlay for linear) ──
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(11.0, 4.4))
    d = sorted(dirt, key=lambda x: x["gamma"])
    gam = [r["gamma"] for r in d]

    a_up_pred = [r["a_up_pred"] for r in d]
    a_down_pred = [r["a_down_pred"] for r in d]
    axa.plot(gam, a_up_pred, "k-", label="theory")
    axa.plot(gam, a_down_pred, "k-")
    axa.plot(gam, [r["a_up"] for r in d], "o", ms=8, color="tab:blue", label="DIRT a_up")
    axa.plot(gam, [r["a_down"] for r in d], "^", ms=8, color="tab:green", label="DIRT a_down")
    if lammps:
        l = sorted(lammps, key=lambda x: x["gamma"])
        axa.plot([r["gamma"] for r in l], [r["a_up"] for r in l], "s", ms=9,
                 mfc="none", color="tab:red", label="LAMMPS a_up")
        axa.plot([r["gamma"] for r in l], [r["a_down"] for r in l], "D", ms=8,
                 mfc="none", color="tab:orange", label="LAMMPS a_down")
    axa.set_xlabel(r"$\gamma_l$"); axa.set_ylabel(r"linear rate (m/s$^2$)")
    axa.set_title("Linear acceleration vs theory / LAMMPS")
    axa.legend(fontsize=8)

    alpha_down_pred = [r["alpha_down_pred"] for r in d]
    alpha_up_pred = [r["alpha_up_pred"] for r in d]
    axb.plot(gam, alpha_down_pred, "k-", label="theory")
    axb.plot(gam, alpha_up_pred, "k-")
    axb.plot(gam, [r["alpha_down"] for r in d], "o", ms=8, color="tab:blue", label=r"DIRT $\alpha_{down}$")
    axb.plot(gam, [r["alpha_up"] for r in d], "^", ms=8, color="tab:green", label=r"DIRT $\alpha_{up}$")
    axb.set_xlabel(r"$\gamma_a$"); axb.set_ylabel(r"angular rate (rad/s$^2$)")
    axb.set_title("Angular acceleration vs theory")
    axb.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "cundall_rates.png"))
    plt.close(fig)
    print(f"\nFigures -> {PLOT_DIR}/cundall_traces.png, cundall_rates.png")


def graph():
    dirt = _load(SWEEP_CSV)
    if not dirt:
        print(f"No {SWEEP_CSV} — run 'start' first.")
        return False
    lammps = _load(LAMMPS_CSV)
    ok = validate(dirt)
    if lammps:
        compare_codes(dirt, lammps)
        # LAMMPS is an independent cross-check: it too must match theory.
        for r in lammps:
            for fit_k, pred_k in (("a_up", "a_up_pred"), ("a_down", "a_down_pred")):
                if _rel(r[fit_k], r[pred_k]) > REL_TOL:
                    print(f"  LAMMPS {fit_k} mismatch at gamma={r['gamma']}: "
                          f"{r[fit_k]:.3f} vs {r[pred_k]:.3f}")
                    ok = False
    else:
        print("\n(no LAMMPS sweep — plotting DIRT vs theory only)")
    plot(dirt, lammps)
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
