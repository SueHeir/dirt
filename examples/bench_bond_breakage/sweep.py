#!/usr/bin/env python3
"""
Bond-breakage / plasticity sweep benchmark driver.

This benchmark closes the "no breakage bench" gap: it drives the bonded-particle
model (BPM) through two quantitatively-gated experiments that exercise the
`dirt_bond` breakage + plasticity machinery end-to-end (config parse -> per-bond
threshold sampling at the actual bond length -> plastic envelope -> in-loop
breakage check -> recorded failure), reusing the existing `fiber_bond` example
binary (no new core code, no new example binary).

Two groups, each with a quantitative theory gate:

  GROUP A — Guo 2018 trilinear bending plasticity, fully-plastic moment cap
  -------------------------------------------------------------------------
  A pinned 11-bead fiber is loaded with the built-in three-step transverse tip
  schedule (fiber_bond/main.rs::apply_three_step_load, activated by the
  "bending_plastic_guo" output-dir tag). The bending channel uses the literal
  `guo_trilinear` envelope, which caps the bond moment at the fully-plastic
  moment (Guo et al. 2018, Chem. Eng. Sci. 175, 118-129, Eq. 31):

        M^p = (4/3) * sigma_0 * r_b^3

  Sweeping the yield stress sigma_0 (0.5, 1.0, 1.5 MPa) scales M^p linearly
  (0.67, 1.33, 2.0 mN*m; all below the moment the fixed tip schedule can drive
  into the middle bond, so every case reaches its cap). The gate reconstructs the peak bond moment from the
  recorded kinematics (M = K_bend * (theta_bend - theta_p_bend)) and checks it
  plateaus at M^p for every sigma_0.

  GROUP B — seeded Weibull weakest-link CDF
  ------------------------------------------------------------
  Sixty independently seeded axial-stress Weibull breakage realizations run the
  same 10-bond fiber. Each run still checks the deterministic weakest-bond
  prediction from the sampled `bond_thresholds.csv`; the ensemble then compares
  the measured first-break strain distribution against the analytical
  weakest-link Weibull CDF:

        F_min(eps) = 1 - exp[-N_bonds * (E*eps/lambda)^m]
        lambda = mean / Gamma(1 + 1/m)

  The gate requires every per-seed first break to match its sampled weakest-link
  prediction within 5% and the ensemble Kolmogorov-Smirnov statistic to stay
  below 0.18.

Commands (run from anywhere):
    python3 examples/bench_bond_breakage/sweep.py generate   # write per-case configs + geometry
    python3 examples/bench_bond_breakage/sweep.py start      # build + run all sims -> CSV
    python3 examples/bench_bond_breakage/sweep.py graph      # validate + plot (exit 1 on FAIL)
    python3 examples/bench_bond_breakage/sweep.py            # all three phases, in order

No LAMMPS needed: both gates are closed-form (Guo Eq. 31 and the weakest-link
Weibull CDF). Reference values come from theory, not back-fit.
"""

import os
import sys
import csv
import math
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE_BIN = "fiber_bond"   # reuse the existing BPM validation binary

SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
GUO_CSV = os.path.join(DATA_DIR, "guo_trilinear.csv")
WEIBULL_CSV = os.path.join(DATA_DIR, "weibull_cdf.csv")

# ── GROUP A: Guo trilinear bending, fully-plastic moment cap ─────────────────
N_B = 11                  # 11-bead pinned fiber (matches the built-in tip schedule, TIP_TAG=10)
SP_B = 4.0e-3             # bead spacing (m); spaced fiber, bonds still form (bond_tolerance=2.1)
R_B = 1.0e-3              # bead radius (m)
E_B = 1.0e9              # Pa
POISSON_B = 0.25
G_B = 4.0e8              # Pa
RHO_B = 2500.0
SIGMA0_LIST = [0.5e6, 1.0e6, 1.5e6]   # yield stress -> M^p = 0.67, 1.33, 2.0 mN*m
DT_B = 1.0e-7
STEPS_B = 120000          # run through the first load cycle's hold window


def guo_Mp(sigma_0, r_b):
    """Guo 2018 Eq. 31 fully-plastic bending moment cap."""
    return (4.0 / 3.0) * sigma_0 * r_b ** 3


# ── Validation tolerances (theory gates) ─────────────────────────────────────
# Guo trilinear peak moment: relative tolerance about M^p.
TOL_GUO_REL = 0.06
# Weibull weakest-link: Kolmogorov-Smirnov gate against the closed-form CDF.
WEIBULL_SEEDS = list(range(1001, 1061))
WEIBULL_N = 11
WEIBULL_LF = 2.0e-2
WEIBULL_RADIUS = 1.0e-3
WEIBULL_E = 1.0e9
WEIBULL_G = 4.0e8
WEIBULL_MEAN = 5.0e6
WEIBULL_MODULUS = 5.0
WEIBULL_L_CALIB = WEIBULL_LF / (WEIBULL_N - 1)
WEIBULL_PULL_V = 0.1
WEIBULL_DT = 1.0e-7
WEIBULL_STEPS = 30000
TOL_WEIBULL_PER_SEED_REL = 0.05
TOL_WEIBULL_KS = 0.18


# ── geometry / config writers ────────────────────────────────────────────────
def _straight_fiber_csv(path, n, spacing):
    with open(path, "w") as f:
        f.write("# straight fiber along x: x, y, z\n")
        for i in range(n):
            f.write(f"{i * spacing:.10e}, 0.0, 0.0\n")


GUO_CONFIG = """# Guo 2018 trilinear bending, sigma_0 = {s0:.3e} Pa, M^p = {mp:.4e} N*m
# NOTE: the output dir MUST contain "bending_plastic_guo" to activate the
# built-in three-step tip-load schedule in fiber_bond/main.rs.
[comm]
processors_x = 1
processors_y = 1
processors_z = 1
[domain]
x_low  = -0.005
x_high =  0.045
y_low  = -0.005
y_high =  0.005
z_low  = -0.040
z_high =  0.010
boundary_x = "fixed"
boundary_y = "fixed"
boundary_z = "fixed"
[neighbor]
skin_fraction = 1.1
bin_size = 0.008
every = 1
[dem]
contact_model = "hertz"
[[dem.materials]]
name = "bpm"
youngs_mod = {e:.6e}
poisson_ratio = {nu}
restitution = 0.9
friction = 0.0
[[particles.insert]]
source = "file"
file = "{csv}"
format = "csv"
material = "bpm"
radius = {radius:.6e}
density = {rho}
columns = {{ x = 0, y = 1, z = 2 }}
[bonds]
auto_bond = true
bond_tolerance = 2.1
ghost_cutoff_multiplier = 4.0
bond_radius_ratio = 1.0
youngs_modulus = {e:.6e}
shear_modulus  = {g:.6e}
beta_normal  = 1.0
beta_shear   = 1.0
beta_twist   = 1.0
beta_bending = 1.0
[bonds.plasticity.bending]
kind = "guo_trilinear"
yield_stress = {s0:.6e}
[[group]]
name = "anchor"
region = {{ type = "block", min = [-5e-4, -5e-4, -5e-4], max = [5e-4, 5e-4, 5e-4] }}
dynamic = false
[[freeze]]
group = "anchor"
[[viscous]]
group = "all"
gamma = 0.1
[output]
dir = "{outdir}"
[run]
steps = {steps}
thermo = {thermo}
dt = {dt:.6e}
"""


WEIBULL_CONFIG = """# Seeded Weibull axial-stress breakage realization.
# Seed {seed}; analytical weakest-link CDF:
# F_min(eps) = 1 - exp[-N_bonds * (E*eps/lambda)^m],
# lambda = mean / Gamma(1 + 1/m).
[comm]
processors_x = 1
processors_y = 1
processors_z = 1
[domain]
x_low  = -0.005
x_high =  0.025
y_low  = -0.005
y_high =  0.005
z_low  = -0.005
z_high =  0.005
boundary_x = "fixed"
boundary_y = "fixed"
boundary_z = "fixed"
[neighbor]
skin_fraction = 1.1
bin_size = 0.005
every = 1
[dem]
contact_model = "hertz"
[[dem.materials]]
name = "bpm"
youngs_mod = {e:.6e}
poisson_ratio = 0.25
restitution = 0.9
friction = 0.0
[[particles.insert]]
source = "file"
file = "{csv}"
format = "csv"
material = "bpm"
radius = {radius:.6e}
density = 2500.0
columns = {{ x = 0, y = 1, z = 2 }}
[bonds]
auto_bond = true
bond_tolerance = 1.001
bond_radius_ratio = 1.0
youngs_modulus = {e:.6e}
shear_modulus  = {g:.6e}
beta_normal  = 1.0
beta_shear   = 1.0
beta_twist   = 1.0
beta_bending = 1.0
seed = {seed}
[bonds.breakage]
kind = "axial_stress"
tensile = {{ kind = "weibull", mean = {mean:.6e}, m = {modulus}, l_calib = {l_calib:.6e} }}
[[group]]
name = "left_end"
region = {{ type = "block", min = [-5e-4, -5e-4, -5e-4], max = [5e-4, 5e-4, 5e-4] }}
dynamic = false
[[group]]
name = "right_end"
region = {{ type = "block", min = [0.0195, -5e-4, -5e-4], max = [0.0205, 5e-4, 5e-4] }}
dynamic = false
[[freeze]]
group = "left_end"
[[move_linear]]
group = "right_end"
vx = {vx:.6e}
[output]
dir = "{outdir}"
[run]
steps = {steps}
thermo = {thermo}
dt = {dt:.6e}
"""


def guo_case_dir(s0):
    # Dir name embeds the "bending_plastic_guo" tag so the tip-load schedule fires.
    return os.path.join(SWEEP_DIR, f"bending_plastic_guo_s{s0:g}")


def _write_guo_case(s0):
    cdir = guo_case_dir(s0)
    os.makedirs(cdir, exist_ok=True)
    csv_path = os.path.join(cdir, "fiber.csv")
    _straight_fiber_csv(csv_path, N_B, SP_B)
    cfg = GUO_CONFIG.format(
        s0=s0, mp=guo_Mp(s0, R_B), e=E_B, nu=POISSON_B, g=G_B, rho=RHO_B,
        radius=R_B, csv=csv_path, outdir=cdir, steps=STEPS_B, thermo=STEPS_B, dt=DT_B,
    )
    with open(os.path.join(cdir, "config.toml"), "w") as f:
        f.write(cfg)
    return cdir


def weibull_case_dir(seed):
    return os.path.join(SWEEP_DIR, f"weibull_seed_{seed}")


def _write_weibull_case(seed):
    cdir = weibull_case_dir(seed)
    os.makedirs(cdir, exist_ok=True)
    csv_path = os.path.join(cdir, "fiber.csv")
    _straight_fiber_csv(csv_path, WEIBULL_N, WEIBULL_LF / (WEIBULL_N - 1))
    cfg = WEIBULL_CONFIG.format(
        seed=seed, e=WEIBULL_E, g=WEIBULL_G, radius=WEIBULL_RADIUS,
        csv=csv_path, mean=WEIBULL_MEAN, modulus=WEIBULL_MODULUS,
        l_calib=WEIBULL_L_CALIB, vx=WEIBULL_PULL_V, outdir=cdir,
        steps=WEIBULL_STEPS, thermo=WEIBULL_STEPS, dt=WEIBULL_DT,
    )
    with open(os.path.join(cdir, "config.toml"), "w") as f:
        f.write(cfg)
    return cdir


# ── generate ─────────────────────────────────────────────────────────────────
def generate():
    n = 0
    for s0 in SIGMA0_LIST:
        _write_guo_case(s0)
        n += 1
    for seed in WEIBULL_SEEDS:
        _write_weibull_case(seed)
        n += 1
    print(f"Generated {n} configs under {SWEEP_DIR}")


# ── start (build + run) ──────────────────────────────────────────────────────
def _run(cdir):
    config = os.path.join(cdir, "config.toml")
    res = os.path.join(cdir, "data", "fiber_bond.csv")
    if os.path.exists(res):
        os.remove(res)
    log = os.path.join(cdir, "run.log")
    with open(log, "w") as lf:
        proc = subprocess.run(
            ["cargo", "run", "--release", "--example", EXAMPLE_BIN,
             "--no-default-features", "--features", "precision-double", "--", config],
            cwd=REPO_ROOT, stdout=lf, stderr=subprocess.STDOUT,
        )
    if proc.returncode != 0 or not os.path.isfile(res):
        return None
    with open(res) as f:
        return list(csv.DictReader(f))


def _fval(row, key):
    try:
        return float(row[key])
    except (KeyError, ValueError):
        return float("nan")


def _measure_guo_peak_moment(rows):
    """Peak reconstructed bond moment M = K_bend*(theta_bend - theta_p_bend)."""
    if not rows:
        return None
    k_bend = _fval(rows[-1], "k_bend")
    peak = 0.0
    for r in rows:
        th_e = _fval(r, "dth_bend_y_mid") - _fval(r, "theta_p_bend_y_mid")
        if math.isnan(th_e):
            continue
        peak = max(peak, abs(k_bend * th_e))
    return peak


def _load_thresholds(cdir):
    path = os.path.join(cdir, "data", "bond_thresholds.csv")
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]


def _measure_weibull(rows, cdir):
    if not rows:
        return None
    thresholds = _load_thresholds(cdir)
    if not thresholds:
        return None
    initial_bonds = int(_fval(rows[0], "bond_count"))
    if initial_bonds <= 0:
        return None
    length0 = _fval(rows[0], "length0")
    idx = next((i for i, r in enumerate(rows)
                if int(_fval(r, "bond_count")) < initial_bonds), None)
    if idx is None:
        return None
    eps_at = (_fval(rows[idx], "length_global") - length0) / length0
    eps_prev = eps_at if idx == 0 else (_fval(rows[idx - 1], "length_global") - length0) / length0
    eps_meas = 0.5 * (eps_prev + eps_at)
    sigma_min = min(t["thr0"] for t in thresholds)
    eps_pred = sigma_min / WEIBULL_E
    return {
        "n_bonds": initial_bonds,
        "eps_break_pred": eps_pred,
        "eps_break_meas": eps_meas,
        "sigma_min": sigma_min,
    }


def weibull_cdf(eps, n_bonds):
    scale = WEIBULL_MEAN / math.gamma(1.0 + 1.0 / WEIBULL_MODULUS)
    x = max(0.0, WEIBULL_E * eps / scale)
    return 1.0 - math.exp(-n_bonds * x ** WEIBULL_MODULUS)


def weibull_quantile(prob, n_bonds):
    scale = WEIBULL_MEAN / math.gamma(1.0 + 1.0 / WEIBULL_MODULUS)
    p = min(max(prob, 1.0e-12), 1.0 - 1.0e-12)
    return (scale / WEIBULL_E) * (-math.log(1.0 - p) / n_bonds) ** (1.0 / WEIBULL_MODULUS)


def _write_csv(path, fields, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})


def start():
    os.makedirs(DATA_DIR, exist_ok=True)
    print(f"Building {EXAMPLE_BIN} (release, precision-double)...", flush=True)
    subprocess.run(["cargo", "build", "--release", "--example", EXAMPLE_BIN,
                    "--no-default-features", "--features", "precision-double"],
                   cwd=REPO_ROOT, check=True)

    # GROUP A — Guo trilinear bending, M^p plateau.
    guo_rows = []
    print("\n[Group A] Guo 2018 trilinear bending, fully-plastic moment cap")
    for s0 in SIGMA0_LIST:
        cdir = guo_case_dir(s0)
        if not os.path.isfile(os.path.join(cdir, "config.toml")):
            print(f"  sigma_0={s0:g}: missing config — run 'generate' first.")
            continue
        print(f"  sigma_0={s0/1e6:4.1f} MPa", end="  ", flush=True)
        data = _run(cdir)
        peak = _measure_guo_peak_moment(data) if data is not None else None
        if peak is None:
            print("RUN FAILED")
            continue
        mp = guo_Mp(s0, R_B)
        guo_rows.append({"sigma_0": s0, "Mp_theory": mp, "peak_moment": peak})
        print(f"peak_M={peak:.4e}  M^p={mp:.4e}  ratio={peak/mp:.3f}")

    if guo_rows:
        _write_csv(GUO_CSV, ["sigma_0", "Mp_theory", "peak_moment"], guo_rows)

    # GROUP B — seeded Weibull weakest-link CDF.
    weibull_rows = []
    print("\n[Group B] seeded Weibull breakage realizations, weakest-link CDF")
    for seed in WEIBULL_SEEDS:
        cdir = weibull_case_dir(seed)
        if not os.path.isfile(os.path.join(cdir, "config.toml")):
            print(f"  seed={seed}: missing config — run 'generate' first.")
            continue
        print(f"  seed={seed}", end="  ", flush=True)
        data = _run(cdir)
        meas = _measure_weibull(data, cdir) if data is not None else None
        if meas is None:
            print("RUN/BREAK FAILED")
            continue
        rel = abs(meas["eps_break_meas"] - meas["eps_break_pred"]) / meas["eps_break_pred"]
        row = {"seed": seed, **meas, "rel_err": rel}
        weibull_rows.append(row)
        print(f"eps_meas={meas['eps_break_meas']:.5f}  eps_pred={meas['eps_break_pred']:.5f}  err={rel:.1%}")

    if weibull_rows:
        _write_csv(WEIBULL_CSV, ["seed", "n_bonds", "eps_break_pred", "eps_break_meas",
                                 "sigma_min", "rel_err"], weibull_rows)
    print(f"\nWrote {GUO_CSV} ({len(guo_rows)} rows), "
          f"{WEIBULL_CSV} ({len(weibull_rows)} rows)")


# ── graph (validate + plot) ──────────────────────────────────────────────────
def _load(path):
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        return [{k: (float(v) if k not in () else v) for k, v in r.items()}
                for r in csv.DictReader(f)]


def validate(guo_rows, weibull_rows):
    ok = True
    print("\n=== Bond-breakage / plasticity sweep validation ===")

    # ── Group A gate: Guo trilinear M^p = (4/3) sigma_0 r_b^3 ──
    print(f"\n[A] Guo trilinear peak moment vs M^p = (4/3) sigma_0 r_b^3  (rel tol {TOL_GUO_REL:.0%})")
    print(f"  {'sigma_0(MPa)':>12}{'M^p(mN*m)':>12}{'peak(mN*m)':>12}{'ratio':>8}  note")
    if len(guo_rows) != len(SIGMA0_LIST):
        print(f"  MISSING CASES: got {len(guo_rows)}/{len(SIGMA0_LIST)} Guo runs")
        ok = False
    for r in sorted(guo_rows, key=lambda x: x["sigma_0"]):
        mp = r["Mp_theory"]
        peak = r["peak_moment"]
        ratio = peak / mp if mp else 0.0
        note = "" if abs(ratio - 1.0) <= TOL_GUO_REL else "MOMENT CAP MISMATCH"
        if note:
            ok = False
        print(f"  {r['sigma_0']/1e6:>12.1f}{mp*1e3:>12.4f}{peak*1e3:>12.4f}{ratio:>8.3f}  {note}")

    # ── Group B gate: per-run weakest-link + distribution CDF ──
    print(f"\n[B] Weibull weakest-link break strain CDF  "
          f"(per-seed tol {TOL_WEIBULL_PER_SEED_REL:.0%}, KS tol {TOL_WEIBULL_KS:.2f})")
    if len(weibull_rows) != len(WEIBULL_SEEDS):
        print(f"  MISSING CASES: got {len(weibull_rows)}/{len(WEIBULL_SEEDS)} seeded Weibull runs")
        ok = False
    if weibull_rows:
        per_seed_max = max(r["rel_err"] for r in weibull_rows)
        n_bonds = int(round(weibull_rows[0]["n_bonds"]))
        eps_sorted = sorted(r["eps_break_meas"] for r in weibull_rows)
        ks = 0.0
        for i, eps in enumerate(eps_sorted, start=1):
            f = weibull_cdf(eps, n_bonds)
            ks = max(ks, abs(i / len(eps_sorted) - f),
                     abs((i - 1) / len(eps_sorted) - f))
        note = []
        if per_seed_max > TOL_WEIBULL_PER_SEED_REL:
            note.append("PER-SEED BREAK MISMATCH")
            ok = False
        if ks > TOL_WEIBULL_KS:
            note.append("KS CDF MISMATCH")
            ok = False
        print(f"  realizations={len(weibull_rows)}, bonds/run={n_bonds}, "
              f"max per-seed err={per_seed_max:.1%}, KS D={ks:.3f}  {'; '.join(note)}")
        print(f"  empirical eps range: {eps_sorted[0]:.5f} .. {eps_sorted[-1]:.5f}; "
              f"analytical median={weibull_quantile(0.5, n_bonds):.5f}")
    else:
        ok = False

    print("\nRESULT:", "ALL CHECKS PASSED" if ok else "CHECKS FAILED")
    return ok


def plot(guo_rows, weibull_rows):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"(matplotlib unavailable, skipping plots: {e})")
        return
    os.makedirs(PLOT_DIR, exist_ok=True)
    plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150, "font.size": 11})

    # Fig 1: Guo peak moment vs sigma_0.
    if guo_rows:
        rs = sorted(guo_rows, key=lambda x: x["sigma_0"])
        s0 = [r["sigma_0"] / 1e6 for r in rs]
        peak = [r["peak_moment"] * 1e3 for r in rs]
        xline = [0.0, max(s0) * 1.05]
        yline = [(4.0 / 3.0) * (x * 1e6) * R_B ** 3 * 1e3 for x in xline]
        fig, ax = plt.subplots(figsize=(6.2, 4.5))
        ax.plot(xline, yline, "k--", label=r"$M^p=\frac{4}{3}\sigma_0 r_b^3$ (Guo Eq. 31)")
        ax.plot(s0, peak, "s", ms=8, color="C3", label="DIRT (peak bond moment)")
        ax.set_xlabel(r"yield stress $\sigma_0$ (MPa)")
        ax.set_ylabel(r"peak bond moment (mN$\cdot$m)")
        ax.set_title("Guo trilinear fully-plastic moment cap")
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(PLOT_DIR, "guo_trilinear_moment.png"))
        plt.close(fig)

    # Fig 2: Weibull CDF and QQ plot.
    if weibull_rows:
        rs = sorted(weibull_rows, key=lambda x: x["eps_break_meas"])
        n_bonds = int(round(rs[0]["n_bonds"]))
        eps = [r["eps_break_meas"] for r in rs]
        empirical = [(i - 0.5) / len(rs) for i in range(1, len(rs) + 1)]
        x_max = max(max(eps) * 1.08, weibull_quantile(0.995, n_bonds))
        xline = [x_max * i / 240 for i in range(241)]
        cdf = [weibull_cdf(x, n_bonds) for x in xline]
        q_theory = [weibull_quantile(p, n_bonds) for p in empirical]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.5))
        ax1.plot(xline, cdf, "k--", label="analytical weakest-link CDF")
        ax1.step(eps, [i / len(rs) for i in range(1, len(rs) + 1)],
                 where="post", color="C0", label="DIRT empirical CDF")
        ax1.set_xlabel(r"first-break strain $\varepsilon$")
        ax1.set_ylabel("CDF")
        ax1.set_title("Seeded Weibull breakage CDF")
        ax1.legend()
        ax2.plot(q_theory, eps, "o", ms=5, label="seeded realizations")
        lim = [0.0, max(max(q_theory), max(eps)) * 1.08]
        ax2.plot(lim, lim, "k--", label="1:1")
        ax2.set_xlim(lim)
        ax2.set_ylim(lim)
        ax2.set_xlabel("analytical quantile")
        ax2.set_ylabel("measured break strain")
        ax2.set_title("Weibull QQ comparison")
        ax2.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(PLOT_DIR, "weibull_cdf_qq.png"))
        plt.close(fig)

    print(f"Figures -> {PLOT_DIR}/")


def graph():
    guo_rows = _load(GUO_CSV)
    weibull_rows = _load(WEIBULL_CSV)
    if not guo_rows and not weibull_rows:
        print(f"No results in {DATA_DIR} — run 'start' first.")
        return False
    ok = validate(guo_rows, weibull_rows)
    plot(guo_rows, weibull_rows)
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
