#!/usr/bin/env python3
"""
Kharaz, Gorham & Salman (2001) oblique-impact replication driver.

Reproduces the *experimental protocol* of Kharaz, Gorham & Salman, "An
experimental study of the elastic rebound of spheres", Powder Technology 120
(2001) 281-291: a 5 mm aluminium-oxide sphere strikes a thick soda-lime **glass**
anvil (a fully elastic response, normal restitution e_n ~= 0.98) at a **fixed
impact speed V_i = 3.85 m/s**, with the **angle of incidence swept** from near
the surface normal to near grazing. The classic Kharaz figures are the three
rebound/spin curves plotted against the incidence angle:

    (1) rebound angle              Theta_r(Theta_i)
    (2) tangential restitution     e_t = v_t' / v_t  (centre-velocity ratio)
    (3) non-dimensional rebound spin   R*omega' / V_i

DIRT geometry.  DIRT fires the sphere at a **flat, frictional glass anvil** (a
real `dirt_wall` z-plane, normal +z) — the same geometry Kharaz used (a thick
glass block).  A flat wall keeps the contact normal exactly +z throughout the
collision at every incidence angle, so the normal restitution stays
angle-independent right up to grazing (a frozen *sphere* partner would curve the
projectile at grazing and depress the apparent e_n).  The wall is infinite, so no
aiming is needed.  The normal contact model is validated separately by
`bench_hertz_rebound`; the tangential Mindlin+Coulomb model against Maw (1976)
theory + LAMMPS by `bench_oblique_impact`.

Independent reference (NOT self-consistent).  In the gross-sliding regime the
contact point slides throughout the collision and the rebound is fixed by the
exact rigid-body impulse relations (no free parameters beyond Kharaz's measured
mu and e_n):

    v_n' = e_n v_n
    v_t' = v_t - mu (1 + e_n) v_n           => e_t = 1 - mu(1+e_n)/tan(Theta_i)
    omega' = 5 mu (1 + e_n) v_n / (2 R)      => R omega'/V_i = 5 mu(1+e_n) v_n /(2 V_i)
    tan(Theta_r) = v_t'/v_n'

valid while  tan(Theta_i) > (7/2) mu (1 + e_n)   (Theta measured from the normal;
here the boundary is Theta_i ~= 32.5 deg).  This is exactly the sliding-branch
behaviour Kharaz's glass-anvil data confirmed.  Below that boundary the contact
sticks/micro-slips (Maw's S-curve); there DIRT is checked qualitatively (spin
below the sliding line, contact-point tangential-velocity reversal).

Kharaz's reported experimental scalars used as anchors (Powder Tech. 120, 2001):
    normal restitution   e_n = 0.98   (glass, ~angle-independent)
    sliding friction     mu  = 0.092
    impact speed         V_i = 3.85 m/s ,   5 mm alumina spheres.

NOTE on data.  The paper's per-point experimental *scatter* for the glass anvil
lives only in the paywalled figures (no open-access copy); this driver therefore
validates against the exact sliding kinematics + Maw theory that those points
confirmed, plus Kharaz's reported e_n and mu.  If the digitised points become
available, drop them into `data/kharaz_experiment.csv` and they will be overlaid.

Commands (run from anywhere):
    python3 examples/bench_kharaz_oblique/sweep.py generate  # write per-angle configs
    python3 examples/bench_kharaz_oblique/sweep.py start     # build + run -> CSV
    python3 examples/bench_kharaz_oblique/sweep.py graph     # validate + plot
    python3 examples/bench_kharaz_oblique/sweep.py           # all three, in order
"""

import os
import sys
import csv
import math
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_kharaz_oblique"   # flat-wall anvil recorder (main.rs)

SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
SWEEP_CSV = os.path.join(DATA_DIR, "kharaz_sweep.csv")
# Optional digitised experimental points; kept at the (tracked) example root so a
# future digitisation lands in git (data/ is gitignored). See README.
EXP_CSV = os.path.join(SCRIPT_DIR, "kharaz_experiment.csv")

# ── Kharaz (2001) glass-anvil conditions ─────────────────────────────────────
V_I = 3.85            # m/s   fixed impact speed
E_N = 0.98            # -     normal restitution (Kharaz, glass)
MU = 0.092            # -     sliding friction  (Kharaz)
NU = 0.23             # -     Poisson ratio (alumina)
YOUNGS_MOD = 380.0e9  # Pa    alumina
DENSITY = 4000.0      # kg/m^3
RADIUS = 0.005        # m     (rebound curves are radius-independent; alumina R)
DT = 1.0e-7
GAP = 3.0e-4          # m     initial clearance between sphere surface and wall
STEPS = 40000         # cover the slowest (largest-angle => smallest v_n) descent

# Incidence angles from the surface normal (deg). 0 = normal impact, 90 = grazing.
THETA_DEG = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80]

# Sliding/sticking boundary (from the normal): tan(Theta_i) = 7/2 mu (1+e_n).
THETA_SLIDE_DEG = math.degrees(math.atan(3.5 * MU * (1.0 + E_N)))

# ── DIRT config template (obliquely-launched sphere + flat frictional wall) ──
TOML_TEMPLATE = """[comm]
processors_x = 1
processors_y = 1
processors_z = 1
[domain]
x_low = -0.05
x_high = 0.05
y_low = -0.02
y_high = 0.02
z_low = 0.0
z_high = 0.05
boundary_x = "fixed"
boundary_y = "fixed"
boundary_z = "fixed"
[neighbor]
skin_fraction = 1.1
bin_size = 0.012
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
velocity_x = {vx:.6e}
velocity_z = -{vz:.6e}
region = {{ type = "block", min = [-1.0e-6, -1.0e-6, {zlo:.8e}], max = [1.0e-6, 1.0e-6, {zhi:.8e}] }}
[[wall]]
point_x = 0.0
point_y = 0.0
point_z = 0.0
normal_x = 0.0
normal_y = 0.0
normal_z = 1.0
material = "alumina"
[output]
dir = "{outdir}"
[run]
steps = {steps}
thermo = {steps}
dt = {dt:.6e}
"""


def vel_components(theta_deg):
    """Return (v_t, v_n) for a fixed impact speed V_I at incidence theta (from normal)."""
    th = math.radians(theta_deg)
    return V_I * math.sin(th), V_I * math.cos(th)


def case_tag(theta_deg):
    return f"th_{theta_deg:g}"


def case_dir(theta_deg):
    return os.path.join(SWEEP_DIR, case_tag(theta_deg))


def _config(theta_deg, outdir):
    v_t, v_n = vel_components(theta_deg)
    # Sphere centre starts a clearance GAP above the wall (contact at z = R).
    z0 = RADIUS + GAP
    return TOML_TEMPLATE.format(
        youngs=YOUNGS_MOD, nu=NU, e_n=E_N, mu=MU, radius=RADIUS, density=DENSITY,
        vx=v_t, vz=v_n, zlo=z0 - 1e-9, zhi=z0 + 1e-9, outdir=outdir, steps=STEPS, dt=DT,
    )


# ── generate ─────────────────────────────────────────────────────────────────
def generate():
    n = 0
    for th in THETA_DEG:
        cdir = case_dir(th)
        os.makedirs(cdir, exist_ok=True)
        with open(os.path.join(cdir, "config.toml"), "w") as f:
            f.write(_config(th, cdir))
        n += 1
    print(f"Generated {n} Kharaz per-angle configs under {SWEEP_DIR}")


# ── start (build once, run each angle) ───────────────────────────────────────
SWEEP_FIELDS = [
    "theta_deg", "theta_meas_deg", "v_t", "v_n",
    "e_n", "e_t", "beta_cp", "spin_nd", "theta_r_deg",
    "e_t_slide", "spin_nd_slide", "theta_r_slide_deg", "sliding",
]


def _run_dirt(cdir):
    config = os.path.join(cdir, "config.toml")
    res = os.path.join(cdir, "data", "oblique_results.csv")
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
    with open(res) as f:
        return next(csv.DictReader(f))


def _row(theta_deg, r):
    """Reduce one DIRT oblique_results.csv row to Kharaz rebound/spin quantities."""
    v_n = float(r["vn_impact"]); v_t = float(r["vt_impact"])
    vn2 = float(r["vn_rebound"]); vt2 = float(r["vt_rebound"])
    wy = float(r["omega_y_rebound"]); R = float(r["radius"])

    theta_meas = math.degrees(math.atan2(v_t, v_n))
    e_n = vn2 / v_n
    e_t = vt2 / v_t                             # centre tangential restitution
    beta_cp = -(vt2 - R * wy) / v_t             # contact-point (Maw beta)
    spin_nd = abs(R * wy) / V_I
    theta_r = math.degrees(math.atan2(vt2, vn2))

    # Exact rigid-body sliding reference (independent of DIRT).
    vt_gs = v_t - MU * (1.0 + E_N) * v_n
    vn_gs = E_N * v_n
    om_gs = 5.0 * MU * (1.0 + E_N) * v_n / (2.0 * R)
    e_t_gs = vt_gs / v_t
    spin_gs = abs(R * om_gs) / V_I
    theta_r_gs = math.degrees(math.atan2(vt_gs, vn_gs))
    sliding = math.tan(math.radians(theta_meas)) > 3.5 * MU * (1.0 + E_N)

    return {
        "theta_deg": theta_deg, "theta_meas_deg": theta_meas, "v_t": v_t, "v_n": v_n,
        "e_n": e_n, "e_t": e_t, "beta_cp": beta_cp, "spin_nd": spin_nd, "theta_r_deg": theta_r,
        "e_t_slide": e_t_gs, "spin_nd_slide": spin_gs, "theta_r_slide_deg": theta_r_gs,
        "sliding": 1 if sliding else 0,
    }


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
                    "--no-default-features", "--features", "precision-double"],
                   cwd=REPO_ROOT, check=True)
    rows = []
    n = len(THETA_DEG)
    for i, th in enumerate(THETA_DEG, 1):
        cdir = case_dir(th)
        if not os.path.isfile(os.path.join(cdir, "config.toml")):
            print(f"  [{i:2d}/{n}] missing config for theta={th} — run 'generate' first.")
            continue
        print(f"  [{i:2d}/{n}] theta={th:>3}deg", end="  ", flush=True)
        r = _run_dirt(cdir)
        if r is None:
            print("DIRT FAILED")
            continue
        row = _row(th, r)
        rows.append(row)
        print(f"e_n={row['e_n']:.3f}  e_t={row['e_t']:+.3f}  R w/V={row['spin_nd']:.3f}  "
              f"Theta_r={row['theta_r_deg']:.1f}deg  ({'slide' if row['sliding'] else 'stick'})")
    if not rows:
        print("\nERROR: no DIRT results collected.")
        sys.exit(1)
    _write_csv(SWEEP_CSV, SWEEP_FIELDS, rows)
    print(f"\n{len(rows)}/{n} cases -> {SWEEP_CSV}")


# ── graph (validate + plot) ──────────────────────────────────────────────────
# Tolerances. Kharaz report agreement of their glass-anvil data with the
# theory to within experimental scatter of a few percent; we require DIRT to
# reproduce the exact sliding kinematics well inside that.
TOL_EN = 0.02        # |e_n - 0.98|  (Kharaz measured value)
TOL_EN_SPREAD = 0.01  # e_n must be ~angle-independent (elastic glass)
TOL_ET = 0.03        # |e_t - sliding|      in the sliding regime
TOL_SPIN = 0.03      # |R w/V - sliding|    in the sliding regime
TOL_THR = 2.0        # |Theta_r - sliding| (deg) in the sliding regime
SQRT_5_6 = 0.9128709291752768


def _tsuji_alpha(e):
    """Tsuji-Tanaka-Ishida Hertz damping polynomial used by DIRT/LAMMPS."""
    return (1.2728 - 4.2783 * e + 11.087 * e**2 - 22.348 * e**3
            + 27.467 * e**4 - 18.022 * e**5 + 4.8218 * e**6)


def _hertz_beta_for_cor(e):
    if e >= 0.9999:
        return 0.0
    return _tsuji_alpha(max(1.0e-3, min(0.9999, e))) / math.sqrt(5.0)


PSI_PREF = 2.0 * (1.0 - NU) / (MU * (2.0 - NU))


def maw_beta_for_theta(theta_deg, dt=2.0e-8):
    """Independent Hertz-Mindlin oblique-impact reference in Maw variables."""
    v_t0, v_n0 = vel_components(theta_deg)
    if v_t0 <= 0.0 or v_n0 <= 0.0:
        return float("nan")

    mass = 4.0 / 3.0 * math.pi * RADIUS**3 * DENSITY
    r_eff = RADIUS
    e_eff = 1.0 / ((1.0 - NU * NU) / YOUNGS_MOD)
    g_eff = 1.0 / (2.0 * (2.0 - NU) * (1.0 + NU) / YOUNGS_MOD)
    beta_n = _hertz_beta_for_cor(E_N)

    overlap = 0.0
    overlap_rate = v_n0
    v_s = v_t0
    spring = 0.0
    for step in range(2_000_000):
        if step > 0 and overlap <= 0.0 and overlap_rate < 0.0:
            break
        if overlap > 0.0:
            sdr = math.sqrt(overlap * r_eff)
            k_n = 4.0 / 3.0 * e_eff * sdr
            s_n = 2.0 * e_eff * sdr
            k_t = 8.0 * g_eff * sdr
            f_n = k_n * overlap + 2.0 * beta_n * SQRT_5_6 * math.sqrt(s_n * mass) * overlap_rate
            f_n = max(0.0, f_n)

            spring += v_s * dt
            gamma_t = 2.0 * SQRT_5_6 * beta_n * math.sqrt(k_t * mass)
            f_t = k_t * spring + gamma_t * v_s
            f_t_max = MU * f_n
            if abs(f_t) > f_t_max:
                f_t = math.copysign(f_t_max, f_t)
                spring = (f_t - gamma_t * v_s) / k_t if k_t > 0.0 else 0.0
        else:
            f_n = 0.0
            f_t = 0.0

        overlap_rate += (-f_n / mass) * dt
        v_s += (-3.5 * f_t / mass) * dt
        overlap += overlap_rate * dt
    return -v_s / v_t0


def maw_kharaz_row(theta_deg):
    beta = maw_beta_for_theta(theta_deg)
    v_t, v_n = vel_components(theta_deg)
    e_t = 1.0 - (1.0 + beta) / 3.5
    spin_nd = (5.0 / 7.0) * (1.0 + beta) * v_t / V_I
    theta_r = math.degrees(math.atan2(e_t * v_t, E_N * v_n))
    return {"theta_deg": theta_deg, "beta_cp": beta, "e_t": e_t,
            "spin_nd": spin_nd, "theta_r_deg": theta_r}


def maw_kharaz_curve(theta_min=2.0, theta_max=82.0, n=180):
    return [maw_kharaz_row(theta_min + (theta_max - theta_min) * i / (n - 1))
            for i in range(n)]


def _load(path):
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        out = []
        for r in csv.DictReader(f):
            out.append({k: (float(v) if v not in ("", None) else float("nan"))
                        for k, v in r.items()})
        return out


def validate(rows):
    print("\n=== Kharaz (2001) oblique-impact validation ===")
    print(f"  V_i={V_I} m/s  e_n={E_N}  mu={MU}  nu={NU}   "
          f"sliding boundary Theta_i>{THETA_SLIDE_DEG:.1f}deg")
    print(f"  {'Th_i':>6}{'e_n':>8}{'e_t':>8}{'e_t*':>8}{'Rw/V':>8}{'Rw/V*':>8}"
          f"{'Th_r':>8}{'Th_r*':>8}  regime  note")
    ok = True
    e_ns = [r["e_n"] for r in rows]
    for r in sorted(rows, key=lambda x: x["theta_meas_deg"]):
        note = ""
        # (1) normal restitution matches Kharaz's measured 0.98 and is constant.
        if abs(r["e_n"] - E_N) > TOL_EN:
            note = "e_n OFF"; ok = False
        # (2) in the sliding regime the rebound/spin curves must equal the exact
        #     rigid-body kinematics Kharaz's data confirmed.
        if r["sliding"]:
            if abs(r["e_t"] - r["e_t_slide"]) > TOL_ET:
                note = (note + " e_t").strip(); ok = False
            if abs(r["spin_nd"] - r["spin_nd_slide"]) > TOL_SPIN:
                note = (note + " spin").strip(); ok = False
            if abs(r["theta_r_deg"] - r["theta_r_slide_deg"]) > TOL_THR:
                note = (note + " Theta_r").strip(); ok = False
            regime = "slide"
        else:
            # (3) sticking/micro-slip: spin must stay at or below the sliding line
            #     (micro-slip develops less spin than full sliding).
            if r["spin_nd"] > r["spin_nd_slide"] + TOL_SPIN:
                note = (note + " spin>slide").strip(); ok = False
            regime = "stick"
        print(f"  {r['theta_meas_deg']:>6.1f}{r['e_n']:>8.3f}{r['e_t']:>+8.3f}"
              f"{r['e_t_slide']:>+8.3f}{r['spin_nd']:>8.3f}{r['spin_nd_slide']:>8.3f}"
              f"{r['theta_r_deg']:>8.1f}{r['theta_r_slide_deg']:>8.1f}  {regime:>6}  {note}")
    e_spread = max(e_ns) - min(e_ns)
    if e_spread > TOL_EN_SPREAD:
        print(f"  e_n spread {e_spread:.4f} > {TOL_EN_SPREAD} — not angle-independent"); ok = False
    print(f"\n  e_n spread across sweep: {e_spread:.4f}  (Kharaz glass: ~angle-independent, 0.98)")
    n_slide = sum(r["sliding"] for r in rows)
    print(f"  sliding-regime cases checked against exact rigid-body kinematics: {n_slide}/{len(rows)}")
    deltas = [abs(r["beta_cp"] - maw_beta_for_theta(r["theta_meas_deg"])) for r in rows]
    max_delta = max(deltas) if deltas else float("nan")
    if max_delta > 0.04:
        print(f"  max |delta beta| vs Maw/Hertz-Mindlin = {max_delta:.4f} > 0.04"); ok = False
    else:
        print(f"  max |delta beta| vs Maw/Hertz-Mindlin = {max_delta:.4f}")
    print("RESULT:", "ALL CHECKS PASSED" if ok else "CHECKS FAILED")
    return ok


def plot(rows, exp):
    os.makedirs(PLOT_DIR, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150, "font.size": 10})

    d = sorted(rows, key=lambda x: x["theta_meas_deg"])
    th = [r["theta_meas_deg"] for r in d]
    maw = maw_kharaz_curve()
    th_m = [r["theta_deg"] for r in maw]
    # The rigid-body sliding relations only hold in the sliding regime; draw the
    # reference there (plus the boundary point) rather than extrapolating it into
    # the sticking/micro-slip region where it is unphysical.
    ds = [r for r in d if r["sliding"]]
    th_s = [r["theta_meas_deg"] for r in ds]

    def exp_pts(key):
        return ([r["theta_deg"] for r in exp if key in r and not math.isnan(r[key])],
                [r[key] for r in exp if key in r and not math.isnan(r[key])])

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    # (a) rebound angle
    ax = axes[0, 0]
    ax.plot(th, [r["theta_r_deg"] for r in d], "o-", label="DIRT")
    ax.plot(th_m, [r["theta_r_deg"] for r in maw], "-", color="tab:green",
            lw=1.7, label="Maw (1976)")
    ax.plot(th_s, [r["theta_r_slide_deg"] for r in ds], "k:", label="rigid sliding")
    if exp:
        x, y = exp_pts("theta_r_deg")
        if x:
            ax.plot(x, y, "r^", label="Kharaz exp.")
    ax.axvline(THETA_SLIDE_DEG, color="gray", lw=0.6, ls="--")
    ax.set_xlabel(r"incidence angle $\Theta_i$ (deg, from normal)")
    ax.set_ylabel(r"rebound angle $\Theta_r$ (deg)")
    ax.set_title("(a) Rebound angle"); ax.legend()

    # (b) tangential restitution
    ax = axes[0, 1]
    ax.plot(th, [r["e_t"] for r in d], "o-", label=r"DIRT $e_t=v_t'/v_t$")
    ax.plot(th_m, [r["e_t"] for r in maw], "-", color="tab:green",
            lw=1.7, label="Maw (1976)")
    ax.plot(th_s, [r["e_t_slide"] for r in ds], "k:", label="rigid sliding")
    if exp:
        x, y = exp_pts("e_t")
        if x:
            ax.plot(x, y, "r^", label="Kharaz exp.")
    ax.axhline(0, color="gray", lw=0.5)
    ax.axvline(THETA_SLIDE_DEG, color="gray", lw=0.6, ls="--")
    ax.set_xlabel(r"incidence angle $\Theta_i$ (deg)")
    ax.set_ylabel(r"tangential restitution $e_t$")
    ax.set_title("(b) Tangential restitution"); ax.legend()

    # (c) non-dimensional rebound spin
    ax = axes[1, 0]
    ax.plot(th, [r["spin_nd"] for r in d], "o-", label=r"DIRT $R\omega'/V_i$")
    ax.plot(th_m, [r["spin_nd"] for r in maw], "-", color="tab:green",
            lw=1.7, label="Maw (1976)")
    ax.plot(th_s, [r["spin_nd_slide"] for r in ds], "k:", label="rigid sliding")
    if exp:
        x, y = exp_pts("spin_nd")
        if x:
            ax.plot(x, y, "r^", label="Kharaz exp.")
    ax.axvline(THETA_SLIDE_DEG, color="gray", lw=0.6, ls="--")
    ax.set_xlabel(r"incidence angle $\Theta_i$ (deg)")
    ax.set_ylabel(r"non-dim rebound spin $R\omega'/V_i$")
    ax.set_title("(c) Rebound angular velocity"); ax.legend()

    # (d) normal restitution vs Kharaz's measured 0.98
    ax = axes[1, 1]
    ax.plot(th, [r["e_n"] for r in d], "o-", label="DIRT")
    ax.axhline(E_N, color="r", lw=1.0, ls="--", label=f"Kharaz exp. $e_n$={E_N}")
    ax.set_ylim(0.9, 1.02)
    ax.set_xlabel(r"incidence angle $\Theta_i$ (deg)")
    ax.set_ylabel(r"normal restitution $e_n$")
    ax.set_title("(d) Normal restitution"); ax.legend()

    fig.suptitle("Kharaz, Gorham & Salman (2001) oblique impact — 5 mm alumina on glass, "
                 r"$V_i$=3.85 m/s", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(PLOT_DIR, "kharaz_rebound_spin.png")
    fig.savefig(out)
    plt.close(fig)
    print(f"\nFigure -> {out}")


def graph():
    rows = _load(SWEEP_CSV)
    if not rows:
        print(f"No {SWEEP_CSV} — run 'start' first.")
        return False
    exp = _load(EXP_CSV)
    if exp:
        print(f"(overlaying {len(exp)} digitised experimental points from {EXP_CSV})")
    ok = validate(rows)
    plot(rows, exp)
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
