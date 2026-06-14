#!/usr/bin/env python3
"""
Fiber-crossover Coulomb-friction benchmark driver.

Two perpendicular bonded-sphere fibers cross at a single contact. The lower
fiber is frozen; the upper fiber is dragged tangentially at constant velocity
with its height held fixed, imposing a known Hertzian normal load N at the
crossover. The tangential contact force rises as the Mindlin spring loads, then
plateaus at the Coulomb limit F_slide = mu * N. This validates inter-fiber
CONTACT + FRICTION — distinct from the intra-fiber bonds the fiber_bond
examples cover.

We sweep the imposed normal load N (by varying the upper fiber's height, hence
the crossover overlap) at fixed mu, and validate that the sliding tangential
force F_slide equals mu * N. Because N is measured directly (the summed +z
contact reaction on the upper fiber), the validation is the ratio
F_slide / N == mu and is insensitive to the exact Hertz overlap.

Commands (from anywhere):
    python3 examples/bench_fiber_crossover/sweep.py generate   # write per-case configs
    python3 examples/bench_fiber_crossover/sweep.py start      # build + run all sims -> CSV
    python3 examples/bench_fiber_crossover/sweep.py graph      # validate + plot
    python3 examples/bench_fiber_crossover/sweep.py            # all three, in order

LAMMPS is not used here (the validation is the exact mu*N Coulomb limit); the
example runs and validates entirely against theory.

Reference: Coulomb friction limit F_t <= mu * F_n; Hertz normal contact
F_n = (4/3) E* sqrt(R*) delta^(3/2). Mindlin tangential spring with a mu|F_n|
slider cap (Mindlin & Deresiewicz, 1953).
"""

import os
import sys
import csv
import shutil
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_fiber_crossover"

SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
SWEEP_CSV = os.path.join(DATA_DIR, "sweep.csv")        # one row per case (N, F_slide)
TRACE_CSV = os.path.join(DATA_DIR, "trace.csv")        # full F_t(displacement) for one case

# LAMMPS binary candidates (unused here, kept for structural consistency).
LAMMPS_BINS = ["lmp_serial", "lmp", "lmp_mpi", "lammps"]

# ── Material / geometry ──────────────────────────────────────────────────────
R = 1.0e-3            # sphere radius (m)
SPACING = 2.0 * R     # centre-to-centre spacing within a fiber
N_SPHERES = 7         # spheres per fiber
YOUNGS = 1.0e7        # Pa (soft polymer fiber)
POISSON = 0.3
DENSITY = 1200.0      # kg/m^3
MU = 0.4              # inter-fiber sliding friction coefficient
DRAG_V = 1.0e-3       # tangential drag speed (m/s)
DT = 2.0e-6           # s
STEPS = 30000         # drag distance = DRAG_V * STEPS * DT = 60 um

# ── Sweep: imposed crossover overlap (um) -> different normal load N ──────────
# delta in metres. Larger overlap -> larger Hertzian N. F_slide should track muN.
OVERLAP_UM_LIST = [8.0, 12.0, 16.0, 20.0, 26.0, 32.0, 40.0]
TRACE_OVERLAP_UM = 20.0   # representative case for the F_t(displacement) figure

# Plateau window (fraction of the run) over which F_slide and N are averaged.
# After the spring saturates (~first few um) and before the contact drifts far.
PLATEAU_LO_FRAC = 0.20
PLATEAU_HI_FRAC = 0.45


# ── helpers ──────────────────────────────────────────────────────────────────
def case_tag(ov_um):
    return f"ov_{ov_um:g}um"


def case_dir(ov_um):
    return os.path.join(SWEEP_DIR, case_tag(ov_um))


def find_lammps():
    for b in LAMMPS_BINS:
        if shutil.which(b):
            return shutil.which(b)
    return None


def _half_offsets(n, spacing):
    """Centred coordinates for n spheres at the given spacing."""
    return [(i - (n - 1) / 2.0) * spacing for i in range(n)]


def _lower_fiber_csv():
    ys = _half_offsets(N_SPHERES, SPACING)
    lines = ["# lower fiber: spheres along y at z=0 (frozen)", "# x, y, z"]
    lines += [f"0.0, {y:.6e}, 0.0" for y in ys]
    return "\n".join(lines) + "\n"


def _upper_fiber_csv(ov_m):
    """Upper fiber along x, at z = 2R - overlap (imposes crossover normal load)."""
    z = SPACING - ov_m
    xs = _half_offsets(N_SPHERES, SPACING)
    lines = [f"# upper fiber: spheres along x at z={z:.6e} (overlap {ov_m*1e6:g} um)", "# x, y, z"]
    lines += [f"{x:.6e}, 0.0, {z:.6e}" for x in xs]
    return "\n".join(lines) + "\n"


def _bonds_data():
    """Explicit intra-fiber bond list (no cross-fiber bonds)."""
    lines = ["# intra-fiber bonds only (tags 0-6 lower, 7-13 upper)", ""]
    nb = 2 * (N_SPHERES - 1)
    lines += [f"{nb} bonds", "1 bond types", "", "Bonds", ""]
    bid = 1
    for base in (0, N_SPHERES):
        for k in range(N_SPHERES - 1):
            a, b = base + k, base + k + 1
            lines.append(f"{bid} 1 {a} {b}")
            bid += 1
    return "\n".join(lines) + "\n"


CONFIG_TEMPLATE = """[comm]
processors_x = 1
processors_y = 1
processors_z = 1
[domain]
x_low  = -0.012
x_high =  0.012
y_low  = -0.010
y_high =  0.010
z_low  = -0.004
z_high =  0.006
boundary_x = "fixed"
boundary_y = "fixed"
boundary_z = "fixed"
[neighbor]
skin_fraction = 1.2
bin_size = 0.004
every = 1
[dem]
contact_model = "hertz"
[[dem.materials]]
name = "polymer"
youngs_mod = {youngs:.6e}
poisson_ratio = {poisson}
restitution = 0.5
friction = {mu}
[[particles.insert]]
source = "file"
file = "{lower_csv}"
format = "csv"
material = "polymer"
radius = {radius:.6e}
density = {density}
columns = {{ x = 0, y = 1, z = 2 }}
[[particles.insert]]
source = "file"
file = "{upper_csv}"
format = "csv"
material = "polymer"
radius = {radius:.6e}
density = {density}
columns = {{ x = 0, y = 1, z = 2 }}
[bonds]
auto_bond = false
file = "{bonds_file}"
format = "lammps_data"
bond_radius_ratio = 1.0
youngs_modulus = {youngs:.6e}
shear_modulus  = {shear:.6e}
beta_normal  = 1.0
beta_shear   = 1.0
beta_twist   = 1.0
beta_bending = 1.0
[[group]]
name = "lower_fiber"
region = {{ type = "block", min = [-0.0015, -0.008, -0.0015], max = [0.0015, 0.008, 0.0015] }}
dynamic = false
[[group]]
name = "upper_fiber"
region = {{ type = "block", min = [-0.008, -0.0015, 0.0005], max = [0.008, 0.0015, 0.0035] }}
dynamic = false
[[freeze]]
group = "lower_fiber"
[[move_linear]]
group = "upper_fiber"
vx = {drag_v:.6e}
vy = 0.0
vz = 0.0
[output]
dir = "{outdir}"
[run]
steps = {steps}
thermo = {thermo}
dt = {dt:.6e}
"""


def _write_case(ov_um, mu=MU):
    cdir = case_dir(ov_um)
    os.makedirs(cdir, exist_ok=True)
    ov_m = ov_um * 1e-6
    lower_csv = os.path.join(cdir, "lower_fiber.csv")
    upper_csv = os.path.join(cdir, "upper_fiber.csv")
    bonds_file = os.path.join(cdir, "bonds.data")
    with open(lower_csv, "w") as f:
        f.write(_lower_fiber_csv())
    with open(upper_csv, "w") as f:
        f.write(_upper_fiber_csv(ov_m))
    with open(bonds_file, "w") as f:
        f.write(_bonds_data())
    shear = YOUNGS / (2.0 * (1.0 + POISSON))
    cfg = CONFIG_TEMPLATE.format(
        youngs=YOUNGS, poisson=POISSON, mu=mu, radius=R, density=DENSITY,
        shear=shear, lower_csv=lower_csv, upper_csv=upper_csv,
        bonds_file=bonds_file, drag_v=DRAG_V, outdir=cdir, steps=STEPS,
        thermo=STEPS, dt=DT,
    )
    with open(os.path.join(cdir, "config.toml"), "w") as f:
        f.write(cfg)
    return cdir


# ── generate ─────────────────────────────────────────────────────────────────
def generate():
    n = 0
    for ov in OVERLAP_UM_LIST:
        _write_case(ov)
        n += 1
    print(f"Generated {n} crossover configs under {SWEEP_DIR}")


# ── start ────────────────────────────────────────────────────────────────────
SWEEP_FIELDS = ["overlap_um", "N", "F_slide", "ratio", "mu"]
TRACE_FIELDS = ["displacement", "f_normal", "f_tangential", "overlap"]


def _run_dirt(cdir):
    config = os.path.join(cdir, "config.toml")
    res = os.path.join(cdir, "data", "fiber_crossover_results.csv")
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
    with open(res) as f:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]


def _plateau_average(rows):
    """Average N (=f_normal) and F_slide (=|f_tangential|) over the plateau window."""
    n = len(rows)
    lo = int(n * PLATEAU_LO_FRAC)
    hi = int(n * PLATEAU_HI_FRAC)
    window = rows[lo:hi] if hi > lo else rows
    fns = [r["f_normal"] for r in window]
    fts = [abs(r["f_tangential"]) for r in window]
    N = sum(fns) / len(fns)
    F_slide = sum(fts) / len(fts)
    return N, F_slide


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

    if find_lammps():
        print("LAMMPS found but not used by this benchmark (theory-only validation).")

    rows = []
    trace_rows = None
    n = len(OVERLAP_UM_LIST)
    for i, ov in enumerate(OVERLAP_UM_LIST, 1):
        cdir = case_dir(ov)
        if not os.path.isfile(os.path.join(cdir, "config.toml")):
            print(f"  [{i}/{n}] missing config for ov={ov} um — run 'generate' first.")
            continue
        print(f"  [{i}/{n}] overlap={ov:<5} um", end="  ", flush=True)
        data = _run_dirt(cdir)
        if data is None:
            print("DIRT FAILED")
            continue
        N, F_slide = _plateau_average(data)
        ratio = F_slide / N if N != 0 else 0.0
        rows.append({"overlap_um": ov, "N": N, "F_slide": F_slide,
                     "ratio": ratio, "mu": MU})
        print(f"N={N:.4e} N  F_slide={F_slide:.4e} N  F_slide/N={ratio:.3f}")
        if abs(ov - TRACE_OVERLAP_UM) < 1e-9:
            trace_rows = data

    if not rows:
        print("\nERROR: no DIRT results collected.")
        sys.exit(1)
    _write_csv(SWEEP_CSV, SWEEP_FIELDS, rows)
    print(f"\nDIRT: {len(rows)}/{n} cases -> {SWEEP_CSV}")
    if trace_rows is not None:
        _write_csv(TRACE_CSV, TRACE_FIELDS, trace_rows)
        print(f"Trace (overlap={TRACE_OVERLAP_UM} um) -> {TRACE_CSV}")


# ── graph (validate + plot) ──────────────────────────────────────────────────
def _load(path):
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]


# Validation tolerances.
RATIO_TOL = 0.05      # |F_slide/N - mu| must be within this of mu
SLOPE_TOL = 0.06      # fitted slope of F_slide vs N must be within this of mu


def _fit_slope_through_origin(xs, ys):
    """Least-squares slope of y = m x (forced through origin)."""
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    return sxy / sxx if sxx != 0 else 0.0


def validate(rows):
    print("\n=== Fiber-crossover Coulomb-friction validation ===")
    print(f"  mu = {MU}")
    print(f"  {'overlap(um)':>11}{'N(N)':>12}{'F_slide(N)':>12}{'F_slide/N':>11}  note")
    ok = True
    for r in sorted(rows, key=lambda x: x["overlap_um"]):
        note = ""
        if abs(r["ratio"] - MU) > RATIO_TOL:
            note = "RATIO MISMATCH"
            ok = False
        print(f"  {r['overlap_um']:>11.1f}{r['N']:>12.4e}{r['F_slide']:>12.4e}"
              f"{r['ratio']:>11.3f}  {note}")

    xs = [r["N"] for r in rows]
    ys = [r["F_slide"] for r in rows]
    slope = _fit_slope_through_origin(xs, ys)
    print(f"\n  Fitted slope F_slide = m*N : m = {slope:.4f}  (expected mu = {MU})")
    if abs(slope - MU) > SLOPE_TOL:
        print(f"  slope deviates from mu by {abs(slope-MU):.4f} > {SLOPE_TOL}")
        ok = False
    print("RESULT:", "PASS" if ok else "FAIL")
    return ok


def plot(rows, trace):
    os.makedirs(PLOT_DIR, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150, "font.size": 11})

    # ── F_slide vs N (slope = mu) ──
    rs = sorted(rows, key=lambda x: x["N"])
    Ns = [r["N"] for r in rs]
    Fs = [r["F_slide"] for r in rs]
    slope = _fit_slope_through_origin(Ns, Fs)
    fig, ax = plt.subplots(figsize=(6.0, 4.5))
    ax.plot(Ns, Fs, "o", label="DIRT (measured)")
    xline = [0.0, max(Ns) * 1.05]
    ax.plot(xline, [MU * x for x in xline], "k--", label=f"Coulomb limit $\\mu N$ ($\\mu={MU}$)")
    ax.plot(xline, [slope * x for x in xline], "r:", label=f"fit slope $m={slope:.3f}$")
    ax.set_xlabel("normal load $N$ (N)")
    ax.set_ylabel("sliding tangential force $F_{slide}$ (N)")
    ax.set_title("Crossover Coulomb friction: $F_{slide} = \\mu N$")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "fslide_vs_N.png"))
    plt.close(fig)

    # ── F_t vs displacement: static rise then muN plateau ──
    if trace:
        d = [r["displacement"] * 1e6 for r in trace]
        ft = [abs(r["f_tangential"]) for r in trace]
        fn = [r["f_normal"] for r in trace]
        muN = [MU * v for v in fn]
        fig, ax = plt.subplots(figsize=(6.5, 4.5))
        ax.plot(d, ft, "-", label=r"$|F_t|$ (DIRT)")
        ax.plot(d, muN, "k--", label=r"$\mu N$ (instantaneous)")
        ax.set_xlabel(r"tangential displacement ($\mu$m)")
        ax.set_ylabel("tangential force (N)")
        ax.set_title(f"Static rise then $\\mu N$ plateau (overlap {TRACE_OVERLAP_UM:g} $\\mu$m)")
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(PLOT_DIR, "ft_vs_displacement.png"))
        plt.close(fig)

    print(f"\nFigures -> {PLOT_DIR}/fslide_vs_N.png" + (", ft_vs_displacement.png" if trace else ""))


def graph():
    rows = _load(SWEEP_CSV)
    if not rows:
        print(f"No {SWEEP_CSV} — run 'start' first.")
        return False
    trace = _load(TRACE_CSV)
    ok = validate(rows)
    plot(rows, trace)
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
