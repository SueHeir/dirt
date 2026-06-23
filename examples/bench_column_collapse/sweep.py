#!/usr/bin/env python3
"""
Granular column-collapse benchmark driver.

Releases a quasi-2D rectangular column of grains (initial width L0, height H) on a
flat floor for a range of aspect ratios a = H/L0, then extracts the final runout
L_f from each settled deposit and checks the dimensionless runout against the
experimental aspect-ratio scaling laws (Lube et al. 2004; Lajeunesse et al. 2004):

    (L_f - L0)/L0 ~ 1.2 * a          (a <~ 2-3, linear regime)
    (L_f - L0)/L0 ~ 1.6 * a^(2/3)    (a >~ 3,   power-law regime)

Commands (from anywhere):
    python3 examples/bench_column_collapse/sweep.py generate   # write per-case configs
    python3 examples/bench_column_collapse/sweep.py start      # build + run all sims -> CSV
    python3 examples/bench_column_collapse/sweep.py graph       # extract L_f, validate + plot
    python3 examples/bench_column_collapse/sweep.py            # all three, in order

The aspect ratio is swept by changing the particle count (settled column height H)
at fixed column width L0. Each DIRT run dumps the rest-state deposit
(x, y, z, radius) to data/<case>/column_collapse_results.csv; this script reads
those, computes L_f as the furthest x where the local deposit height exceeds one
particle diameter, and fits the runout exponent in each regime.

If a LAMMPS binary (lmp_serial / lmp / lmp_mpi / lammps) is on PATH, each aspect
ratio is ALSO run in LAMMPS with the equivalent granular model (pair_style
granular hertz/material ... tangential mindlin ... damping tsuji, same E/nu/e/mu,
gravity, and frictional floor + back + side + removable-gate walls via
fix wall/gran). LAMMPS's final deposit is parsed into the SAME (x, y, z, radius)
form and runout is extracted with the SAME measure_column() the DIRT leg uses, so
the two codes are compared on equal footing and overlaid (open markers) on
plots/runout_scaling.png. LAMMPS is optional: with no binary present, only DIRT
runs and the validation (DIRT-vs-theory) is unchanged.

Outputs:
    sweep/<case>/config.toml            DIRT configs                  (gitignored)
    sweep/<case>/in.lammps              LAMMPS inputs                 (gitignored)
    data/<case>/column_collapse_*.csv   per-case DIRT deposits        (gitignored)
    data/runout.csv                     L0, H, a, L_f per DIRT case   (gitignored)
    data/lammps_results.csv             L0, H, a, L_f per LAMMPS case (gitignored)
    plots/*.png                         final figures                 (tracked)
"""

import os
import sys
import csv
import math
import shutil
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_column_collapse"

SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
RUNOUT_CSV = os.path.join(DATA_DIR, "runout.csv")          # DIRT runout per aspect
LAMMPS_CSV = os.path.join(DATA_DIR, "lammps_results.csv")  # LAMMPS runout per aspect

# LAMMPS binary candidates, in preference order. LAMMPS is optional: if none is
# found, the LAMMPS leg is skipped and only DIRT is run/plotted.
LAMMPS_BINS = ["lmp_serial", "lmp", "lmp_mpi", "lammps"]

# ── Geometry / material (shared by every case) ───────────────────────────────
RADIUS = 0.0015            # m (d = 3 mm; Lajeunesse used ~1–3 mm glass beads)
DENSITY = 2500.0           # kg/m^3 (glass)
L0 = 0.024                 # initial column width [m] (= 8 diameters)
W = 0.009                  # slab width in y [m] (= 3 diameters, quasi-2D)
# Canonical glass-bead (ballotini) material — measured properties, shared across
# all DIRT calibrations (shear/cooling/conduction/collapse). E softened from the
# real ~65 GPa (rigid-grain limit; keeps dt tractable). e and μ_p are measured
# glass–glass values (Wu et al. 2019, Meas. of restitution & friction for glass beads).
YOUNGS_MOD = 7.0e7         # Pa (softened from ~65 GPa real glass)
POISSON = 0.245
RESTITUTION = 0.926        # measured glass–glass COR
FRICTION = 0.16            # measured glass–glass sliding friction
DT = 4.0e-6               # s
SETTLE_STEPS = 80000
COLLAPSE_STEPS = 200000

PACKING = 0.60             # settled solid fraction used to size the particle count

# Aspect ratios to sweep. Spans both regimes (linear a<~2-3, power-law a>~3) so a
# regime change is resolvable.
ASPECTS = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0]

# Validation tolerances on the fitted runout exponent per regime.
EXP_TOL = 0.25             # |fitted exponent - target| pass band
LINEAR_TARGET = 1.0        # (L_f-L0)/L0 ~ a^1   for a <~ 2-3
POWER_TARGET = 2.0 / 3.0   # (L_f-L0)/L0 ~ a^2/3 for a >~ 3
REGIME_SPLIT = 3.0         # aspect ratio dividing the two regimes


def n_particles(aspect):
    """Particle count whose settled column (width L0, slab W, packing PACKING)
    has height H = aspect * L0."""
    h = aspect * L0
    vol_particle = (4.0 / 3.0) * math.pi * RADIUS**3
    return max(1, int(round(PACKING * L0 * W * h / vol_particle)))


def case_tag(aspect):
    return f"a{aspect:g}".replace(".", "p")


def case_dir(aspect):
    return os.path.join(SWEEP_DIR, case_tag(aspect))


def data_case_dir(aspect):
    return os.path.join(DATA_DIR, case_tag(aspect))


# ── DIRT config template ─────────────────────────────────────────────────────
TOML_TEMPLATE = """\
# Auto-generated column-collapse config — aspect a = {aspect}, N = {count}
[comm]
processors_x = 1
processors_y = 1
processors_z = 1

[domain]
x_low = -0.01
x_high = 0.60
y_low = -0.003
y_high = 0.012
z_low = 0.0
z_high = {z_high}
boundary_x = "fixed"
boundary_y = "fixed"
boundary_z = "fixed"

[neighbor]
skin_fraction = 1.1
bin_size = 0.005
every = 1

[gravity]
gx = 0.0
gy = 0.0
gz = -9.81

[dem]
contact_model = "hertz"

[[dem.materials]]
name = "glass"
youngs_mod = {youngs:.6e}
poisson_ratio = {poisson}
restitution = {restitution}
friction = {friction}

[[particles.insert]]
material = "glass"
count = {count}
radius = {radius}
density = {density}
region = {{ type = "block", min = [0.0015, 0.0015, 0.0015], max = [0.0225, 0.0075, {insert_top}] }}

[[wall]]
type = "plane"
point_z = 0.0
normal_z = 1.0
material = "glass"
name = "floor"

[[wall]]
type = "plane"
point_x = 0.0
normal_x = 1.0
material = "glass"
name = "back"

[[wall]]
type = "plane"
point_y = 0.0
normal_y = 1.0
material = "glass"
name = "side_lo"

[[wall]]
type = "plane"
point_y = 0.009
normal_y = -1.0
material = "glass"
name = "side_hi"

[[wall]]
type = "plane"
point_x = 0.024
normal_x = -1.0
material = "glass"
name = "gate"

[output]
dir = "{output_dir}"

[vtp]
interval = 1000000

[[run]]
name = "settle"
steps = {settle_steps}
thermo = 20000
dt = {dt}

[[run]]
name = "collapse"
steps = {collapse_steps}
thermo = 20000
dt = {dt}
"""


def generate():
    os.makedirs(SWEEP_DIR, exist_ok=True)
    for a in ASPECTS:
        cdir = case_dir(a)
        os.makedirs(cdir, exist_ok=True)
        n = n_particles(a)
        h = a * L0
        # Loose insert column ~1.6x the settled height; cap the box height.
        insert_top = min(0.18, max(0.02, 1.6 * h))
        z_high = max(0.2, insert_top + 0.05)
        with open(os.path.join(cdir, "config.toml"), "w") as f:
            f.write(TOML_TEMPLATE.format(
                aspect=a, count=n,
                youngs=YOUNGS_MOD, poisson=POISSON,
                restitution=RESTITUTION, friction=FRICTION,
                radius=RADIUS, density=DENSITY,
                insert_top=f"{insert_top:.4f}", z_high=f"{z_high:.4f}",
                output_dir=cdir, dt=f"{DT:.3e}",
                settle_steps=SETTLE_STEPS, collapse_steps=COLLAPSE_STEPS,
            ))
    print(f"Generated {len(ASPECTS)} configs under {SWEEP_DIR}")


# ── LAMMPS leg (optional cross-code overlay) ─────────────────────────────────
def find_lammps():
    """Return the first available LAMMPS binary on PATH, or None."""
    for b in LAMMPS_BINS:
        path = shutil.which(b)
        if path:
            return path
    return None


# LAMMPS counterpart of the DIRT column collapse. Same material, same geometry,
# same two-stage protocol (settle against a gate, remove the gate, collapse):
#   pair_style granular hertz/material E e nu  -> Young's modulus, restitution,
#       Poisson ratio (E/nu/e identical to DIRT's [dem.materials]).
#   tangential mindlin NULL {damp} {mu}        -> Mindlin tangential spring with
#       k_t = 8 G* sqrt(R* delta) derived from the normal contact (NULL) — exactly
#       DIRT's k_t — Coulomb friction coefficient {mu}.
#   damping tsuji                              -> viscoelastic normal+tangential
#       damping from the restitution e (DIRT uses the same e-> damping mapping).
#   fix grav ... 9.81 vector 0 0 -1            -> gravity g_z = -9.81, matching
#       [gravity] in the DIRT config.
#   fix wall/gran ... zplane/xplane/yplane     -> frictional floor, back wall, and
#       the quasi-2D side walls — same granular model (incl. friction) as the pair
#       style, the LAMMPS analogue of DIRT's frictional dirt_wall planes.
#   fix gate ... xplane NULL {L0}; unfix gate  -> removable gate at x = L0, present
#       during 'settle', unfix-ed at the start of 'collapse' (mirrors
#       Walls::deactivate_by_name on the first collapse step in DIRT).
# Atoms are seeded overlap-free into a tall loose column and settle under gravity,
# the same loose-insert-then-settle that DIRT performs. The final deposit is dumped
# as (id, x, y, z, radius); runout is then extracted with the SAME measure_column().
LMP_TEMPLATE = """\
# Auto-generated LAMMPS input for the column-collapse sweep — aspect a = {aspect}
units           si
atom_style      sphere
dimension       3
boundary        f f f
newton          off
comm_modify     vel yes

region          simbox block {x_low} {x_high} {y_low} {y_high} 0.0 {z_high} units box
create_box      1 simbox

region          colreg block 0.0015 0.0225 0.0015 0.0075 0.0015 {insert_top} units box
create_atoms    1 random {count} {seed} colreg overlap {min_sep} maxtry 500 units box
set             group all diameter {diam}
set             group all density {density}

pair_style      granular
pair_coeff      1 1 hertz/material {E} {e} {nu} tangential mindlin NULL {tdamp} {mu} damping tsuji rolling none twisting none

fix             grav all gravity {g} vector 0 0 -1
fix             floor all wall/gran granular hertz/material {E} {e} {nu} tangential mindlin NULL {tdamp} {mu} damping tsuji rolling none twisting none zplane 0.0 NULL
fix             back all wall/gran granular hertz/material {E} {e} {nu} tangential mindlin NULL {tdamp} {mu} damping tsuji rolling none twisting none xplane 0.0 NULL
fix             sides all wall/gran granular hertz/material {E} {e} {nu} tangential mindlin NULL {tdamp} {mu} damping tsuji rolling none twisting none yplane 0.0 {W}
fix             gate all wall/gran granular hertz/material {E} {e} {nu} tangential mindlin NULL {tdamp} {mu} damping tsuji rolling none twisting none xplane NULL {L0}
fix             integrate all nve/sphere

thermo_modify   lost warn flush yes
timestep        {dt}
thermo          {thermo}

# Stage 1: settle the loose column against the gate.
run             {settle_steps}

# Stage 2: remove the gate; the column collapses and spreads to rest.
unfix           gate
run             {collapse_steps}

write_dump      all custom {dump} id x y z radius modify sort id
"""


def lammps_dump_path(aspect):
    return os.path.join(case_dir(aspect), "lammps_deposit.txt")


def write_lammps_input(path, aspect):
    """Write the LAMMPS input for one aspect ratio (same geometry as DIRT)."""
    n = n_particles(aspect)
    h = aspect * L0
    # Loose insert column. Unlike DIRT's inserter, LAMMPS 'create_atoms random'
    # rejects overlapping placements, so the loose region must be tall enough to
    # hold all N grains — otherwise it silently places fewer than N (skewing the
    # effective column height and the runout).
    footprint = (0.0225 - 0.0015) * (0.0075 - 0.0015)   # x*y of the insert region
    vol_particle = (4.0 / 3.0) * math.pi * RADIUS**3
    # 'create_atoms random' rejects overlaps, so a denser-but-plausible loose pack
    # (~0.45) sets the required loose-column height; the small initial overlaps the
    # min-separation permits are resolved during 'settle'. Size the column to hold
    # all N grains so the count matches DIRT exactly.
    loose_pack = 0.45
    h_needed = n * vol_particle / (loose_pack * footprint)
    insert_top = max(0.02, 1.6 * h, h_needed + 0.0015)
    z_high = insert_top + 0.05
    # min center separation: slightly below d so the loose column packs densely
    # enough to place every grain (overlaps relax in the settle stage).
    min_sep = 0.85 * 2.0 * RADIUS
    with open(path, "w") as f:
        f.write(LMP_TEMPLATE.format(
            aspect=aspect, count=n, seed=12345,
            x_low=-0.01, x_high=0.60, y_low=-0.003, y_high=0.012,
            z_high=f"{z_high:.4f}", insert_top=f"{insert_top:.4f}",
            min_sep=f"{min_sep:.6f}",
            diam=2.0 * RADIUS, density=DENSITY,
            E=f"{YOUNGS_MOD:.6e}", e=RESTITUTION, nu=POISSON,
            tdamp=1.0, mu=FRICTION, g=9.81,
            W=W, L0=L0, dt=f"{DT:.3e}", thermo=40000,
            settle_steps=SETTLE_STEPS, collapse_steps=COLLAPSE_STEPS,
            dump=lammps_dump_path(aspect),
        ))


def lammps_dump_to_csv(dump_path, csv_path):
    """Convert a LAMMPS 'id x y z radius' dump to the same x,y,z,radius CSV that
    the DIRT recorder writes, so measure_column() can read it unchanged."""
    with open(dump_path) as f:
        lines = f.readlines()
    # Find the 'ITEM: ATOMS' header; columns follow it.
    start = None
    cols = []
    for i, line in enumerate(lines):
        if line.startswith("ITEM: ATOMS"):
            cols = line.split()[2:]
            start = i + 1
            break
    if start is None:
        return False
    ix, iy, iz, ir = (cols.index(c) for c in ("x", "y", "z", "radius"))
    with open(csv_path, "w", newline="") as out:
        out.write("x,y,z,radius\n")
        for line in lines[start:]:
            p = line.split()
            if len(p) < len(cols):
                continue
            out.write(f"{p[ix]},{p[iy]},{p[iz]},{p[ir]}\n")
    return True


def run_lammps_sweep(lammps):
    """Run every aspect ratio in LAMMPS, parse each deposit with the SAME
    measure_column() the DIRT leg uses, and return runout rows."""
    rows = []
    for i, a in enumerate(ASPECTS, 1):
        cdir = case_dir(a)
        os.makedirs(cdir, exist_ok=True)
        in_path = os.path.join(cdir, "in.lammps")
        log_path = os.path.join(cdir, "lammps.log")
        dump = lammps_dump_path(a)
        deposit_csv = os.path.join(cdir, "lammps_deposit.csv")
        for stale in (dump, deposit_csv):
            if os.path.isfile(stale):
                os.remove(stale)
        write_lammps_input(in_path, a)
        print(f"  [LAMMPS {i}/{len(ASPECTS)}] a={a:<4} N={n_particles(a)}", flush=True)
        proc = subprocess.run(
            [lammps, "-in", in_path, "-log", log_path],
            cwd=REPO_ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT,
        )
        if proc.returncode != 0 or not os.path.isfile(dump):
            print(f"    a={a}: LAMMPS run failed (see {log_path}).")
            continue
        if not lammps_dump_to_csv(dump, deposit_csv):
            print(f"    a={a}: could not parse LAMMPS dump.")
            continue
        h, lf = measure_column(deposit_csv)
        rows.append({"aspect": a, "L0": L0, "H": h, "L_f": lf,
                     "runout_norm": (lf - L0) / L0})
    return rows


# ── start ────────────────────────────────────────────────────────────────────
def start():
    os.makedirs(DATA_DIR, exist_ok=True)
    print(f"Building {EXAMPLE} (release)...", flush=True)
    env = dict(os.environ)
    # macOS: ensure system libffi is found if the workspace needs it.
    subprocess.run(
        ["cargo", "build", "--release", "--example", EXAMPLE, "--no-default-features"],
        cwd=REPO_ROOT, check=True, env=env,
    )

    for i, a in enumerate(ASPECTS, 1):
        cdir = case_dir(a)
        config = os.path.join(cdir, "config.toml")
        if not os.path.isfile(config):
            print(f"  [{i}/{len(ASPECTS)}] missing {config} — run 'generate' first.")
            continue
        # Wipe stale deposit so old results can't be re-plotted.
        deposit = os.path.join(cdir, "data", "column_collapse_results.csv")
        if os.path.isfile(deposit):
            os.remove(deposit)
        print(f"  [{i}/{len(ASPECTS)}] a={a:<4} N={n_particles(a)}", flush=True)
        log = os.path.join(cdir, "run.log")
        with open(log, "w") as lf:
            subprocess.run(
                ["cargo", "run", "--release", "--example", EXAMPLE,
                 "--no-default-features", "--", config],
                cwd=REPO_ROOT, stdout=lf, stderr=subprocess.STDOUT, env=env,
            )

    rows = []
    for a in ASPECTS:
        deposit = os.path.join(case_dir(a), "data", "column_collapse_results.csv")
        if not os.path.isfile(deposit):
            print(f"  a={a}: no deposit produced.")
            continue
        h, lf = measure_column(deposit)
        rows.append({"aspect": a, "L0": L0, "H": h, "L_f": lf,
                     "runout_norm": (lf - L0) / L0})

    if not rows:
        print("\nERROR: no deposits collected.")
        sys.exit(1)
    os.makedirs(DATA_DIR, exist_ok=True)
    _write_runout(RUNOUT_CSV, rows)
    print(f"\nDIRT:   wrote {len(rows)} runout rows -> {RUNOUT_CSV}")

    # LAMMPS leg — optional cross-code overlay. Skipped entirely with no binary.
    lammps = find_lammps()
    if lammps:
        print(f"LAMMPS: {lammps} — running cross-code overlay.")
        if os.path.isfile(LAMMPS_CSV):
            os.remove(LAMMPS_CSV)
        lrows = run_lammps_sweep(lammps)
        if lrows:
            _write_runout(LAMMPS_CSV, lrows)
            print(f"LAMMPS: wrote {len(lrows)} runout rows -> {LAMMPS_CSV}")
        else:
            print("LAMMPS: no deposits collected — skipping overlay.")
    else:
        print("LAMMPS: not found on PATH — running DIRT only.")


def _write_runout(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["aspect", "L0", "H", "L_f", "runout_norm"])
        w.writeheader()
        w.writerows(rows)


# ── deposit analysis ─────────────────────────────────────────────────────────
def measure_column(deposit_path):
    """Return (H_initial_estimate, L_f) from a settled deposit.

    H is estimated from the particle count and footprint; L_f is the furthest x
    at which the local deposit height (binned in x) exceeds one particle diameter.
    """
    xs, zs, rs = [], [], []
    with open(deposit_path) as f:
        for r in csv.DictReader(f):
            xs.append(float(r["x"]))
            zs.append(float(r["z"]))
            rs.append(float(r["radius"]))
    if not xs:
        return 0.0, L0

    d = 2.0 * (sum(rs) / len(rs))      # mean particle diameter
    # Initial column height from solids volume in the L0 x W footprint.
    n = len(xs)
    vol = n * (4.0 / 3.0) * math.pi * (sum(rs) / len(rs)) ** 3
    h_init = vol / (PACKING * L0 * W)

    # Bin the deposit in x; the local height is the max z (top of column) in the
    # bin. L_f = furthest bin whose height exceeds one diameter.
    bin_w = d
    x_min = min(xs)
    bins = {}
    for x, z in zip(xs, zs):
        b = int((x - x_min) / bin_w)
        bins[b] = max(bins.get(b, 0.0), z)
    lf = L0
    for b, top in bins.items():
        if top >= d:                   # column at least one grain tall here
            x_edge = x_min + (b + 1) * bin_w
            lf = max(lf, x_edge)
    return h_init, lf


# ── graph (validate + plot) ──────────────────────────────────────────────────
def load_runout():
    if not os.path.isfile(RUNOUT_CSV):
        print(f"ERROR: {RUNOUT_CSV} not found.")
        print("Run the sweep first: python3 examples/bench_column_collapse/sweep.py start")
        sys.exit(1)
    with open(RUNOUT_CSV) as f:
        rows = list(csv.DictReader(f))
    if not rows:
        print("ERROR: no runout data.")
        sys.exit(1)
    return rows


def fit_loglog(pairs):
    """Least-squares slope of log(y) vs log(x). Returns (exponent, prefactor)."""
    pts = [(a, y) for a, y in pairs if a > 0 and y > 0]
    if len(pts) < 2:
        return float("nan"), float("nan")
    lx = [math.log(a) for a, _ in pts]
    ly = [math.log(y) for _, y in pts]
    n = len(pts)
    sx, sy = sum(lx), sum(ly)
    sxx = sum(v * v for v in lx)
    sxy = sum(a * b for a, b in zip(lx, ly))
    denom = n * sxx - sx * sx
    if abs(denom) < 1e-30:
        return float("nan"), float("nan")
    slope = (n * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / n
    return slope, math.exp(intercept)


def validate(rows):
    print("=" * 66)
    print("Granular Column-Collapse Runout Validation")
    print("=" * 66)
    print(f"  L0 = {L0*1000:.1f} mm, slab W = {W*1000:.1f} mm, d = {2*RADIUS*1000:.1f} mm")
    print(f"  E = {YOUNGS_MOD:.1e} Pa, e = {RESTITUTION}, mu = {FRICTION}\n")
    print(f"  {'a':>5} {'H[mm]':>8} {'L_f[mm]':>9} {'(Lf-L0)/L0':>12}")

    pairs = []
    for r in rows:
        a = float(r["aspect"])
        h = float(r["H"])
        lf = float(r["L_f"])
        rn = float(r["runout_norm"])
        pairs.append((a, rn))
        print(f"  {a:>5.2f} {h*1000:>8.2f} {lf*1000:>9.2f} {rn:>12.3f}")

    low = [(a, rn) for a, rn in pairs if a <= REGIME_SPLIT]
    high = [(a, rn) for a, rn in pairs if a >= REGIME_SPLIT]
    e_low, _ = fit_loglog(low)
    e_high, _ = fit_loglog(high)

    low_ok = abs(e_low - LINEAR_TARGET) <= EXP_TOL
    high_ok = abs(e_high - POWER_TARGET) <= EXP_TOL

    print()
    print(f"  Linear regime (a <= {REGIME_SPLIT}): fitted exponent = {e_low:.3f} "
          f"(target {LINEAR_TARGET:.2f})  [{'PASS' if low_ok else 'FAIL'}]")
    print(f"  Power regime  (a >= {REGIME_SPLIT}): fitted exponent = {e_high:.3f} "
          f"(target {POWER_TARGET:.2f})  [{'PASS' if high_ok else 'FAIL'}]")

    ok = low_ok and high_ok
    if not ok:
        print()
        print("  NOTE: dirt_wall has no particle-wall sliding friction, so the")
        print("  basal floor cannot arrest the deposit — the column slides into a")
        print("  thin sheet and the runout does not follow the experimental laws.")
        print("  This is a core limitation; see README. Core crates were NOT edited.")
    print("\nALL CHECKS PASSED" if ok else "VALIDATION FAILED (see note above)")
    return ok


def compare_codes(dirt_rows, lammps_rows):
    """Print a per-aspect DIRT-vs-LAMMPS normalized-runout comparison and the
    fitted exponents for both codes."""
    dirt = {float(r["aspect"]): float(r["runout_norm"]) for r in dirt_rows}
    lammps = {float(r["aspect"]): float(r["runout_norm"]) for r in lammps_rows}
    print("\n" + "=" * 58)
    print("Normalized runout (L_f-L0)/L0: DIRT vs LAMMPS")
    print("=" * 58)
    print(f"  {'a':>5} | {'DIRT':>8} {'LAMMPS':>8} | {'diff':>8}")
    for a in sorted(set(dirt) & set(lammps)):
        d, l = dirt[a], lammps[a]
        print(f"  {a:>5.2f} | {d:>8.3f} {l:>8.3f} | {l - d:>+8.3f}")

    def fits(data):
        pairs = [(a, data[a]) for a in sorted(data)]
        low = [(a, v) for a, v in pairs if a <= REGIME_SPLIT]
        high = [(a, v) for a, v in pairs if a >= REGIME_SPLIT]
        return fit_loglog(low)[0], fit_loglog(high)[0]

    dl, dh = fits(dirt)
    ll, lh = fits(lammps)
    print("\n  Fitted exponents:        linear (a<=3)   power (a>=3)")
    print(f"    DIRT   :               {dl:>10.3f}    {dh:>10.3f}")
    print(f"    LAMMPS :               {ll:>10.3f}    {lh:>10.3f}")
    print(f"    targets:               {LINEAR_TARGET:>10.2f}    {POWER_TARGET:>10.3f}")


def plot(rows, lammps_rows=None):
    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"\n(matplotlib/numpy unavailable, skipped plots: {e})")
        return

    os.makedirs(PLOT_DIR, exist_ok=True)
    plt.rcParams.update({"font.size": 12, "figure.dpi": 150, "savefig.dpi": 150})

    a = np.array([float(r["aspect"]) for r in rows])
    rn = np.array([float(r["runout_norm"]) for r in rows])

    # ── Plot 1: normalized runout vs aspect ratio (log-log) with scaling lines.
    fig, ax = plt.subplots(figsize=(7, 5.2))
    ax.plot(a, rn, "o", color="#1f77b4", markersize=7, label="DIRT")
    if lammps_rows:
        la = np.array([float(r["aspect"]) for r in lammps_rows])
        lrn = np.array([float(r["runout_norm"]) for r in lammps_rows])
        ax.plot(la, lrn, "s", color="#d62728", markersize=8,
                markerfacecolor="none", markeredgewidth=1.6, label="LAMMPS")
    aa = np.logspace(math.log10(0.4), math.log10(6.0), 100)
    ax.plot(aa, 1.2 * aa, "k--", linewidth=1.4, label=r"$1.2\,a$ (linear, $a\lesssim3$)")
    ax.plot(aa, 1.6 * aa ** (2.0 / 3.0), "k:", linewidth=1.4,
            label=r"$1.6\,a^{2/3}$ (power, $a\gtrsim3$)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Aspect ratio  a = H / L0")
    ax.set_ylabel(r"Normalized runout  $(L_f - L_0)/L_0$")
    ax.set_title("Column-Collapse Runout vs Aspect Ratio")
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.3)
    fig.savefig(os.path.join(PLOT_DIR, "runout_scaling.png"), bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {PLOT_DIR}/runout_scaling.png")

    # ── Plot 2: deposit-profile snapshot for the representative a = 2 case.
    target = min(rows, key=lambda r: abs(float(r["aspect"]) - 2.0))
    a_t = float(target["aspect"])
    deposit = os.path.join(case_dir(a_t), "data", "column_collapse_results.csv")
    if os.path.isfile(deposit):
        xs, zs = [], []
        with open(deposit) as f:
            for r in csv.DictReader(f):
                xs.append(float(r["x"]) * 1000)
                zs.append(float(r["z"]) * 1000)
        fig, ax = plt.subplots(figsize=(9, 3.2))
        ax.scatter(xs, zs, s=6, color="#ff7f0e")
        ax.axvline(L0 * 1000, color="0.5", linestyle="--", linewidth=1,
                   label=r"$L_0$")
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("z [mm]")
        ax.set_title(f"Deposit profile (a = {a_t:g})")
        ax.set_aspect("equal")
        ax.legend(fontsize=9)
        fig.savefig(os.path.join(PLOT_DIR, "deposit_profile.png"), bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {PLOT_DIR}/deposit_profile.png")


def load_optional(path):
    """Load a runout CSV if it exists, else return []."""
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f))


def graph():
    rows = load_runout()
    lammps_rows = load_optional(LAMMPS_CSV)
    ok = validate(rows)            # DIRT-vs-theory only; LAMMPS never gates PASS.
    if lammps_rows:
        compare_codes(rows, lammps_rows)
    else:
        print(f"\n(no {os.path.basename(LAMMPS_CSV)} — plotting DIRT only)")
    plot(rows, lammps_rows)
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
