#!/usr/bin/env python3
"""
Oblique-impact tangential-contact benchmark driver.

Fires a sphere obliquely at a frozen target across a range of incidence angles
and compares the tangential coefficient of restitution beta(psi1) to Maw, Barber
& Fawcett (1976) theory and (if LAMMPS is available) to LAMMPS's Hertz-Mindlin
granular model. Also traces the per-step contact force through a single impact.

Commands (from anywhere):
    python3 examples/bench_oblique_impact/sweep.py generate   # write per-case configs
    python3 examples/bench_oblique_impact/sweep.py start      # build + run all sims -> CSV
    python3 examples/bench_oblique_impact/sweep.py graph      # validate + plot
    python3 examples/bench_oblique_impact/sweep.py            # all three, in order

If a LAMMPS binary (lmp_serial / lmp / lmp_mpi / lammps) is on PATH, each case is
also run in LAMMPS with the equivalent granular Hertz-Mindlin model and overlaid
on the plots. Both codes aim the projectile dead-center (impact normal = +z) and
decompose in that frame, so the comparison isolates the shared contact physics.
LAMMPS is optional — without it, only DIRT runs.

Outputs:
    sweep/<case>/config.toml    DIRT configs                          (gitignored)
    sweep/<case>/in.lammps      LAMMPS inputs                         (gitignored)
    data/sweep.csv              DIRT sweep results                    (gitignored)
    data/sweep_lammps.csv       LAMMPS sweep results                  (gitignored)
    data/trace_dirt.csv         DIRT per-step contact trace           (gitignored)
    data/trace_lammps.csv       LAMMPS per-step contact trace         (gitignored)
    plots/*.png                 final figures                         (tracked)

Non-dim incidence angle:  psi1   = [2(1-nu)/(mu(2-nu))] * (v_t / v_n)
Tangential restitution:   beta   = -v_s'/v_s,  v_s' = v_t' - R*omega'  (contact point)
Gross-sliding branch:     beta_gs= -1 + 7(1+e_n)(1-nu) / [(2-nu) psi1]

Reference: Maw, Barber & Fawcett, "The oblique impact of elastic spheres",
           Wear 38 (1976) 101-114; Kharaz, Gorham & Salman, Powder Tech. 120 (2001) 281-291.
"""

import os
import sys
import csv
import math
import shutil
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_oblique_impact"

SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
SWEEP_CSV = os.path.join(DATA_DIR, "sweep.csv")              # DIRT sweep
LAMMPS_CSV = os.path.join(DATA_DIR, "sweep_lammps.csv")      # LAMMPS sweep
TRACE_DIRT_CSV = os.path.join(DATA_DIR, "trace_dirt.csv")    # DIRT per-step trace
TRACE_LAMMPS_CSV = os.path.join(DATA_DIR, "trace_lammps.csv")  # LAMMPS per-step trace

# LAMMPS binary candidates, in preference order. LAMMPS is optional.
LAMMPS_BINS = ["lmp_serial", "lmp", "lmp_mpi", "lammps"]

# ── Sweep parameters — Kharaz, Gorham & Salman (2001): alumina, near-elastic ──
NU = 0.23           # Poisson ratio
MU = 0.092          # sliding friction
E_N = 0.98          # nominal normal restitution
YOUNGS_MOD = 380.0e9  # Pa
DENSITY = 4000.0    # kg/m^3
RADIUS = 0.005      # m
V_Z = 1.0           # fixed normal (impact) velocity, m/s
# Tangential velocities span psi1 ~ 0.3 (sticking) -> ~9.5 (gross slip).
V_X_LIST = [0.03, 0.05, 0.08, 0.12, 0.18, 0.25, 0.35, 0.5, 0.7, 1.0]
DT = 1.0e-7
STEPS = 120000      # DIRT: launched from z0=0.065, ~50k steps to contact
TRACE_VX = 0.18     # psi1 ~ 1.7 — representative microslip case for the trace figure

Z0 = 0.065          # DIRT projectile launch height (target at z=0.05)
DESCENT = Z0 - 0.06  # vertical drop to first contact (center reaches z=0.06)
PSI_PREF = 2.0 * (1.0 - NU) / (MU * (2.0 - NU))

# LAMMPS launches close (z=0.0601, gap 1e-4) so it needs far fewer steps.
LMP_STEPS = 8000

# ── DIRT config template (one tangential velocity) ───────────────────────────
TOML_TEMPLATE = """[comm]
processors_x = 1
processors_y = 1
processors_z = 1
[domain]
x_low = -0.3
x_high = 0.3
y_low = -0.02
y_high = 0.02
z_low = 0.0
z_high = 0.2
boundary_x = "fixed"
boundary_y = "fixed"
boundary_z = "fixed"
[neighbor]
skin_fraction = 1.1
bin_size = 0.02
every = 1
[dem]
contact_model = "hertz"
[[dem.materials]]
name = "alumina"
youngs_mod = {youngs:.6e}
poisson_ratio = {nu}
restitution = {e_n}
friction = {mu}
[[particles.insert]]
material = "alumina"
count = 1
radius = {radius}
density = {density}
region = {{ type = "block", min = [-1.0e-6, -1.0e-6, 0.049999], max = [1.0e-6, 1.0e-6, 0.050001] }}
[[particles.insert]]
material = "alumina"
count = 1
radius = {radius}
density = {density}
velocity_x = {vx}
velocity_z = -{v_z}
region = {{ type = "block", min = [{xlo}, -1.0e-6, 0.064999], max = [{xhi}, 1.0e-6, 0.065001] }}
[[group]]
name = "target"
region = {{ type = "block", min = [-0.002, -0.002, 0.045], max = [0.002, 0.002, 0.055] }}
[[freeze]]
group = "target"
[output]
dir = "{outdir}"
[run]
steps = {steps}
thermo = {steps}
dt = {dt:.6e}
"""

# ── LAMMPS templates (aimed dead-center: x0 = -vx*1e-4 so impact normal = +z) ─
LMP_HEAD = """units           si
atom_style      sphere
boundary        f f f
comm_modify     vel yes
region          box block -0.1 0.1 -0.02 0.02 0.0 0.15 units box
create_box      1 box
variable        xstart equal -{vx}*1.0e-4
create_atoms    1 single 0.0 0.0 0.05      units box   # id 1 = target
create_atoms    1 single ${{xstart}} 0.0 0.0601 units box   # id 2 = projectile
set             group all diameter {diam}
set             group all density {density}
pair_style      granular
pair_coeff      1 1 hertz/material {youngs:.6e} {e_n} {nu} tangential mindlin NULL 1.0 {mu} &
                damping coeff_restitution rolling none twisting none
group           target id 1
group           proj   id 2
velocity        proj set {vx} 0.0 -{v_z} units box
fix             hold  target setforce 0.0 0.0 0.0
fix             integ proj   nve/sphere
timestep        {dt:.6e}
"""

LMP_SWEEP_TEMPLATE = LMP_HEAD + """run             {steps}
# Final (post-rebound) projectile translational + angular velocity.
write_dump      proj custom {out} id vx vz omegay modify sort id
"""

LMP_TRACE_TEMPLATE = LMP_HEAD + """# Per-step contact force + state of the projectile.
dump            trc proj custom 1 {out} id x z vx vz fx fz omegay
run             {steps}
"""


# ── helpers ──────────────────────────────────────────────────────────────────
def case_tag(vx):
    return f"vx_{vx:g}"


def case_dir(vx):
    return os.path.join(SWEEP_DIR, case_tag(vx))


def aim_x0(vx):
    """DIRT launch x so the projectile drifts to x=0 at contact (normal = +z)."""
    return -vx * DESCENT / V_Z


def find_lammps():
    for b in LAMMPS_BINS:
        path = shutil.which(b)
        if path:
            return path
    return None


def _dirt_config(vx, outdir, steps):
    x0 = aim_x0(vx)
    return TOML_TEMPLATE.format(
        youngs=YOUNGS_MOD, nu=NU, e_n=E_N, mu=MU, radius=RADIUS, density=DENSITY,
        vx=vx, v_z=V_Z, xlo=x0 - 1e-6, xhi=x0 + 1e-6, outdir=outdir, steps=steps, dt=DT,
    )


# ── generate ─────────────────────────────────────────────────────────────────
def generate():
    n = 0
    for vx in V_X_LIST:
        cdir = case_dir(vx)
        os.makedirs(cdir, exist_ok=True)
        with open(os.path.join(cdir, "config.toml"), "w") as f:
            f.write(_dirt_config(vx, cdir, STEPS))
        n += 1
    # The per-step trace case (single representative velocity).
    tdir = os.path.join(SWEEP_DIR, "trace")
    os.makedirs(tdir, exist_ok=True)
    with open(os.path.join(tdir, "config.toml"), "w") as f:
        f.write(_dirt_config(TRACE_VX, tdir, STEPS))
    print(f"Generated {n} DIRT sweep configs + 1 trace config under {SWEEP_DIR}")


# ── start ────────────────────────────────────────────────────────────────────
SWEEP_FIELDS = ["psi1", "beta", "beta_gross_slip", "e_n", "ke_out_over_in", "v_t", "v_n"]
TRACE_FIELDS = ["overlap", "fn", "ft", "omega"]


def _beta_row(r):
    """Turn one DIRT oblique_results.csv row into a sweep-result dict."""
    v_n = float(r["vn_impact"]); v_t = float(r["vt_impact"])
    vn2 = float(r["vn_rebound"]); vt2 = float(r["vt_rebound"])
    wy = float(r["omega_y_rebound"]); R = float(r["radius"])
    psi1 = PSI_PREF * v_t / v_n
    v_s2 = vt2 - R * wy                        # rebound contact-point tangential
    beta = -v_s2 / v_t
    beta_gs = -1.0 + 7.0 * (1.0 + E_N) * (1.0 - NU) / ((2.0 - NU) * psi1)
    ke_in = 0.5 * (v_t**2 + v_n**2)
    ke_out = 0.5 * (vt2**2 + vn2**2) + 0.5 * (2.0 / 5.0) * R**2 * wy**2
    return {"psi1": psi1, "beta": beta, "beta_gross_slip": beta_gs,
            "e_n": vn2 / v_n, "ke_out_over_in": ke_out / ke_in, "v_t": v_t, "v_n": v_n}


def _run_dirt(cdir, trace=False):
    """Run the DIRT example for a prepared case dir; return its results-row dict."""
    config = os.path.join(cdir, "config.toml")
    res = os.path.join(cdir, "data", "oblique_results.csv")
    for stale in (res, os.path.join(cdir, "contact_trace.csv")):
        if os.path.exists(stale):
            os.remove(stale)
    env = dict(os.environ, DIRT_TRACE="1") if trace else os.environ
    log = os.path.join(cdir, "run.log")
    with open(log, "w") as lf:
        proc = subprocess.run(
            ["cargo", "run", "--release", "--example", EXAMPLE, "--no-default-features", "--", config],
            cwd=REPO_ROOT, stdout=lf, stderr=subprocess.STDOUT, env=env,
        )
    if proc.returncode != 0 or not os.path.isfile(res):
        return None
    with open(res) as f:
        return next(csv.DictReader(f))


def _run_lammps_sweep(lammps, vx, cdir):
    """Run one LAMMPS sweep case; return (vxf, vzf, wyf) or None."""
    in_path = os.path.join(cdir, "in.lammps")
    dump = os.path.join(cdir, "lammps.dump")
    with open(in_path, "w") as f:
        f.write(LMP_SWEEP_TEMPLATE.format(
            vx=vx, v_z=V_Z, diam=2.0 * RADIUS, density=DENSITY,
            youngs=YOUNGS_MOD, e_n=E_N, nu=NU, mu=MU, dt=DT, steps=LMP_STEPS, out=dump))
    proc = subprocess.run([lammps, "-in", in_path], cwd=cdir,
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if proc.returncode != 0 or not os.path.isfile(dump):
        return None
    return _parse_lammps_final(dump)


def _parse_lammps_final(dump):
    with open(dump) as f:
        lines = f.readlines()
    for k, ln in enumerate(lines):
        if ln.startswith("ITEM: ATOMS"):
            p = lines[k + 1].split()
            return float(p[1]), float(p[2]), float(p[3])  # vx, vz, omegay
    return None


def _lammps_sweep_row(vx, vxf, vzf, wyf):
    v_n, v_t = V_Z, vx
    psi1 = PSI_PREF * v_t / v_n
    v_s2 = vxf - RADIUS * wyf
    beta = -v_s2 / v_t
    beta_gs = -1.0 + 7.0 * (1.0 + E_N) * (1.0 - NU) / ((2.0 - NU) * psi1)
    ke_in = 0.5 * (v_t**2 + v_n**2)
    ke_out = 0.5 * (vxf**2 + vzf**2) + 0.5 * (2.0 / 5.0) * RADIUS**2 * wyf**2
    return {"psi1": psi1, "beta": beta, "beta_gross_slip": beta_gs,
            "e_n": vzf / v_n, "ke_out_over_in": ke_out / ke_in, "v_t": v_t, "v_n": v_n}


def _parse_dirt_trace(path):
    """DIRT contact_trace.csv: step,overlap,fn,ft_mag,ft_signed,vt_x,omega_y."""
    rows = []
    with open(path) as f:
        for line in f:
            p = line.split(",")
            if len(p) >= 7:
                rows.append({"overlap": float(p[1]), "fn": float(p[2]),
                             "ft": float(p[4]), "omega": float(p[6])})
    return rows


def _parse_lammps_trace(dump):
    """LAMMPS per-step dump: id x z vx vz fx fz omegay -> overlap/fn/ft/omega rows."""
    rows = []
    lines = open(dump).read().split("\n")
    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: ATOMS"):
            i += 1
            if i < len(lines) and lines[i].strip():
                p = lines[i].split()
                x, z, vx_, vz, fx, fz, wy = (float(v) for v in p[1:8])
                ov = (2.0 * RADIUS) - (z - 0.05)
                if ov > 0.0:
                    rows.append({"overlap": ov, "fn": fz, "ft": fx, "omega": wy})
        i += 1
    return rows


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
    subprocess.run(["cargo", "build", "--release", "--example", EXAMPLE, "--no-default-features"],
                   cwd=REPO_ROOT, check=True)

    lammps = find_lammps()
    print(f"LAMMPS: {lammps}" if lammps else "LAMMPS: not found on PATH — running DIRT only.")

    # DIRT sweep.
    dirt_rows, lmp_rows = [], []
    n = len(V_X_LIST)
    for i, vx in enumerate(V_X_LIST, 1):
        cdir = case_dir(vx)
        if not os.path.isfile(os.path.join(cdir, "config.toml")):
            print(f"  [{i:2d}/{n}] missing config for vx={vx} — run 'generate' first.")
            continue
        print(f"  [{i:2d}/{n}] vx={vx:<5}", end="  ", flush=True)
        r = _run_dirt(cdir)
        if r is None:
            print("DIRT FAILED")
            continue
        row = _beta_row(r)
        dirt_rows.append(row)
        print(f"DIRT psi1={row['psi1']:.2f} beta={row['beta']:+.3f}", end="")
        if lammps:
            lf = _run_lammps_sweep(lammps, vx, cdir)
            if lf:
                lrow = _lammps_sweep_row(vx, *lf)
                lmp_rows.append(lrow)
                print(f"   LAMMPS beta={lrow['beta']:+.3f}", end="")
        print()

    if not dirt_rows:
        print("\nERROR: no DIRT results collected.")
        sys.exit(1)
    _write_csv(SWEEP_CSV, SWEEP_FIELDS, dirt_rows)
    print(f"\nDIRT:   {len(dirt_rows)}/{n} cases -> {SWEEP_CSV}")
    if lmp_rows:
        _write_csv(LAMMPS_CSV, SWEEP_FIELDS, lmp_rows)
        print(f"LAMMPS: {len(lmp_rows)}/{n} cases -> {LAMMPS_CSV}")

    # Per-step trace (single representative velocity).
    print(f"\nTrace case (vx={TRACE_VX}, psi1~{PSI_PREF*TRACE_VX:.1f}):")
    tdir = os.path.join(SWEEP_DIR, "trace")
    if _run_dirt(tdir, trace=True) is not None:
        tr = _parse_dirt_trace(os.path.join(tdir, "contact_trace.csv"))
        _write_csv(TRACE_DIRT_CSV, TRACE_FIELDS, tr)
        print(f"  DIRT   {len(tr)} steps -> {TRACE_DIRT_CSV}")
    if lammps:
        in_path = os.path.join(tdir, "in.lammps_trace")
        dump = os.path.join(tdir, "lammps_trace.dump")
        with open(in_path, "w") as f:
            f.write(LMP_TRACE_TEMPLATE.format(
                vx=TRACE_VX, v_z=V_Z, diam=2.0 * RADIUS, density=DENSITY,
                youngs=YOUNGS_MOD, e_n=E_N, nu=NU, mu=MU, dt=DT, steps=LMP_STEPS, out=dump))
        proc = subprocess.run([lammps, "-in", in_path], cwd=tdir,
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if proc.returncode == 0 and os.path.isfile(dump):
            tr = _parse_lammps_trace(dump)
            _write_csv(TRACE_LAMMPS_CSV, TRACE_FIELDS, tr)
            print(f"  LAMMPS {len(tr)} steps -> {TRACE_LAMMPS_CSV}")


# ── graph (validate + plot) ──────────────────────────────────────────────────
def _load(path):
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]


def validate(rows):
    """Check the gross-slip asymptote, normal-restitution constancy, and energy."""
    print("\n=== Oblique-impact validation ===")
    print(f"  nu={NU}  mu={MU}  e_n(nominal)={E_N}")
    print(f"  {'psi1':>7}{'beta':>9}{'beta_gs':>9}{'e_n':>8}{'KE_out/in':>11}  note")
    ok = True
    e_ns = [r["e_n"] for r in rows]
    for r in sorted(rows, key=lambda x: x["psi1"]):
        note = ""
        # In the gross-slip regime (large psi1) beta must track the rigid branch.
        if r["psi1"] > 5.0 and abs(r["beta"] - r["beta_gross_slip"]) > 0.05:
            note = "GROSS-SLIP MISMATCH"; ok = False
        if r["ke_out_over_in"] > 1.001:
            note = (note + " ENERGY-GAIN").strip(); ok = False
        print(f"  {r['psi1']:>7.3f}{r['beta']:>9.3f}{r['beta_gross_slip']:>9.3f}"
              f"{r['e_n']:>8.3f}{r['ke_out_over_in']:>11.4f}  {note}")
    # Normal restitution must be (nearly) independent of tangential velocity.
    e_spread = max(e_ns) - min(e_ns)
    if e_spread > 0.01:
        print(f"  e_n spread {e_spread:.4f} > 0.01 — normal restitution not constant"); ok = False
    print(f"\n  e_n spread across sweep: {e_spread:.4f} (constant => normal model decoupled from tangential)")
    print("RESULT:", "PASS" if ok else "FAIL")
    return ok


def compare_codes(dirt, lammps):
    lmp = {round(r["psi1"], 1): r for r in lammps}
    print("\n=== DIRT vs LAMMPS ===")
    print(f"  {'psi1':>7}{'DIRT':>9}{'LAMMPS':>9}{'d_beta':>9}")
    mx = 0.0
    for r in sorted(dirt, key=lambda x: x["psi1"]):
        l = lmp.get(round(r["psi1"], 1))
        if not l:
            continue
        d = r["beta"] - l["beta"]; mx = max(mx, abs(d))
        print(f"  {r['psi1']:>7.2f}{r['beta']:>9.3f}{l['beta']:>9.3f}{d:>+9.3f}")
    print(f"  max |d_beta| = {mx:.4f}")


def plot(dirt, lammps, trace_d, trace_l):
    os.makedirs(PLOT_DIR, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150, "font.size": 11})

    # ── beta(psi1): DIRT vs LAMMPS vs gross-slip theory ──
    d = sorted(dirt, key=lambda x: x["psi1"])
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.plot([r["psi1"] for r in d], [r["beta"] for r in d], "o-", label="DIRT (Hertz-Mindlin)")
    if lammps:
        l = sorted(lammps, key=lambda x: x["psi1"])
        ax.plot([r["psi1"] for r in l], [r["beta"] for r in l], "s--", label="LAMMPS (Hertz-Mindlin)")
    # Gross-slip branch, only where it is physically valid (|beta_gs| <= peak).
    gs = [(r["psi1"], r["beta_gross_slip"]) for r in d if r["beta_gross_slip"] <= 0.5]
    if gs:
        ax.plot([p for p, _ in gs], [b for _, b in gs], "k:", label="gross-slip (rigid)")
    ax.axhline(0, color="gray", lw=0.5)
    ax.set_xlabel(r"non-dim incidence angle $\psi_1$")
    ax.set_ylabel(r"tangential restitution $\beta = -v_s'/v_s$")
    ax.set_title("Oblique impact vs Maw (1976) — Kharaz conditions")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "beta_vs_psi1.png"))
    plt.close(fig)

    # ── per-step contact-force trace (normal + tangential loops) ──
    if trace_d:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.5))
        um = 1e6
        ax1.plot([r["overlap"] * um for r in trace_d], [r["fn"] for r in trace_d], "-", label="DIRT")
        ax2.plot([r["overlap"] * um for r in trace_d], [r["ft"] for r in trace_d], "-", label="DIRT")
        if trace_l:
            ax1.plot([r["overlap"] * um for r in trace_l], [r["fn"] for r in trace_l], "--", label="LAMMPS")
            ax2.plot([r["overlap"] * um for r in trace_l], [r["ft"] for r in trace_l], "--", label="LAMMPS")
        ax1.set_xlabel("overlap (µm)"); ax1.set_ylabel("normal force F_n (N)")
        ax1.set_title("Normal"); ax1.legend()
        ax2.axhline(0, color="gray", lw=0.5)
        ax2.set_xlabel("overlap (µm)"); ax2.set_ylabel("tangential force F_t,x (N)")
        ax2.set_title(f"Tangential (ψ₁≈{PSI_PREF*TRACE_VX:.1f})"); ax2.legend()
        fig.suptitle("Per-step contact-force trace: DIRT vs LAMMPS")
        fig.tight_layout()
        fig.savefig(os.path.join(PLOT_DIR, "contact_trace.png"))
        plt.close(fig)
    print(f"\nFigures -> {PLOT_DIR}/beta_vs_psi1.png" + (", contact_trace.png" if trace_d else ""))


def graph():
    dirt = _load(SWEEP_CSV)
    if not dirt:
        print(f"No {SWEEP_CSV} — run 'start' first.")
        return False
    lammps = _load(LAMMPS_CSV)
    trace_d = _load(TRACE_DIRT_CSV)
    trace_l = _load(TRACE_LAMMPS_CSV)
    ok = validate(dirt)
    if lammps:
        compare_codes(dirt, lammps)
    else:
        print("\n(no LAMMPS sweep — plotting DIRT only)")
    plot(dirt, lammps, trace_d, trace_l)
    return ok


# ── dispatch ─────────────────────────────────────────────────────────────────
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
