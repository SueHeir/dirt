#!/usr/bin/env python3
"""
Bond-breakage / plasticity sweep benchmark driver.

This benchmark closes the "no breakage bench" gap: it drives the bonded-particle
model (BPM) through two quantitatively-gated experiments that exercise the
`dirt_bond` breakage + plasticity machinery end-to-end (config parse -> per-bond
threshold sampling at the actual bond length -> plastic envelope -> in-loop
breakage check -> recorded failure), reusing the existing `fiber_bond` example
binary (no new core code, no new example binary).

Two groups, each with a closed-form theory gate:

  GROUP A — crack-band breakage across a 4x mesh-refinement sweep
  -------------------------------------------------------------------------
  A straight fiber is pulled axially, quasi-statically and symmetrically, with a
  piecewise (trilinear-shaped: elastic -> hardening -> perfectly-plastic) axial
  plasticity envelope and an `axial_strain` breakage criterion whose threshold
  is drawn from the *crack-band* length-rescaling law (Bazant 1976; Hillerborg-
  Modeer-Petersson 1976):

        eps_break(L_bond) = eps_yield + (value_ref - eps_yield) * l_ref / L_bond

  Because a uniform axial pull produces a spatially uniform strain, every bond
  reaches its (identical, deterministic) crack-band threshold at the same global
  strain, so the fiber's first-break global strain equals eps_break(L_bond)
  exactly. Refining the mesh (N = 11, 21, 41 beads over a fixed 2 mm fiber ->
  L_bond = 200, 100, 50 um) sweeps L_bond by 4x and the predicted break strain
  moves 0.04 -> 0.06 -> 0.10 (dt = 1e-7 s keeps the finest, stiffest bonds
  stable through the startup transient). The gate checks the *measured* first-break global
  strain against the crack-band closed form for every mesh, and that each bond
  goes plastic (eps_p > 0) before it breaks.

  This validates that the crack-band regularization is applied with the correct
  per-bond length scaling inside the running solver — not just in the unit test
  of `ThresholdDistribution::sample`.

  GROUP B — Guo 2018 trilinear bending plasticity, fully-plastic moment cap
  -------------------------------------------------------------------------
  A pinned 11-bead fiber is loaded with the built-in three-step transverse tip
  schedule (fiber_bond/main.rs::apply_three_step_load, activated by the
  "bending_plastic_guo" output-dir tag). The bending channel uses the literal
  `guo_trilinear` envelope, which caps the bond moment at the fully-plastic
  moment (Guo et al. 2018, Chem. Eng. Sci. 175, 118-129, Eq. 31):

        M^p = (4/3) * sigma_0 * r_b^3

  Sweeping the yield stress sigma_0 (2.0, 2.5, 3.0 MPa) scales M^p linearly
  (2.67, 3.33, 4.0 mN*m; all below the moment the fixed tip schedule can drive
  into the middle bond, so every case reaches its cap). The gate reconstructs the peak bond moment from the
  recorded kinematics (M = K_bend * (theta_bend - theta_p_bend)) and checks it
  plateaus at M^p for every sigma_0.

Commands (run from anywhere):
    python3 examples/bench_bond_breakage/sweep.py generate   # write per-case configs + geometry
    python3 examples/bench_bond_breakage/sweep.py start      # build + run all sims -> CSV
    python3 examples/bench_bond_breakage/sweep.py graph      # validate + plot (exit 1 on FAIL)
    python3 examples/bench_bond_breakage/sweep.py            # all three, in order

No LAMMPS needed: both gates are closed-form (crack-band rescaling law; Guo
Eq. 31). Reference values come from theory, not back-fit.
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
CB_CSV = os.path.join(DATA_DIR, "crackband.csv")     # one row per mesh case
GUO_CSV = os.path.join(DATA_DIR, "guo_trilinear.csv")  # one row per sigma_0 case

# ── GROUP A: crack-band axial breakage across mesh refinement ────────────────
LF = 2.0e-3               # fiber end-to-end length (m), fixed across the sweep
N_LIST = [11, 21, 41]     # bead counts -> L_bond = LF/(N-1) = 200, 100, 50 um
E_AX = 1.0e7              # Pa (soft polymer fiber)
POISSON_AX = 0.3
G_AX = E_AX / (2.0 * (1.0 + POISSON_AX))
RHO_AX = 1200.0
BETA_AX = 0.02            # light bond damping (beta*N << 1 so the strain field equilibrates)
GRIP_V = 0.02            # m/s per grip (symmetric); global strain rate = 2 v / LF = 20 /s
DT_AX = 1.0e-7            # small enough that the finest mesh (L_bond=50um) is stable at startup
STEPS_AX = 64000          # reach ~0.128 global strain (past the 0.10 break of the finest mesh)

# Crack-band threshold parameters (axial strain criterion).
CB_VALUE_REF = 0.06       # threshold at l_ref
CB_L_REF = 1.0e-4         # reference bond length (m)
CB_EPS_YIELD = 0.02       # elastic anchor: only the part above rescales with L_bond
# Piecewise (trilinear-shaped) axial plasticity envelope.
AX_BREAKPOINTS = [0.02, 0.04]
AX_SLOPES = [0.3, 0.0]


def cb_predicted_break_strain(l_bond):
    """Crack-band closed form: eps_yield + (value_ref - eps_yield) * l_ref / L_bond."""
    return CB_EPS_YIELD + (CB_VALUE_REF - CB_EPS_YIELD) * CB_L_REF / l_bond


# ── GROUP B: Guo trilinear bending, fully-plastic moment cap ─────────────────
N_B = 11                  # 11-bead pinned fiber (matches the built-in tip schedule, TIP_TAG=10)
SP_B = 4.0e-3             # bead spacing (m); spaced fiber, bonds still form (bond_tolerance=2.1)
R_B = 1.0e-3              # bead radius (m)
E_B = 1.0e9              # Pa
POISSON_B = 0.25
G_B = 4.0e8              # Pa
RHO_B = 2500.0
SIGMA0_LIST = [2.0e6, 2.5e6, 3.0e6]   # yield stress -> M^p = 2.67, 3.33, 4.0 mN*m (all reach the cap)
DT_B = 1.0e-7
STEPS_B = 120000          # through the first load cycle's hold window (~12 ms), reaches the plateau


def guo_Mp(sigma_0, r_b):
    """Guo 2018 Eq. 31 fully-plastic bending moment cap."""
    return (4.0 / 3.0) * sigma_0 * r_b ** 3


# ── Validation tolerances (theory gates) ─────────────────────────────────────
# Crack-band: relative tolerance, with an absolute floor of ~1.2 sample-strain
# steps to account for finite recording cadence (record_every=200 steps).
TOL_CB_REL = 0.08
SAMPLE_DSTRAIN = (2.0 * GRIP_V / LF) * 200 * DT_AX   # global strain advance per recorded sample
TOL_CB_ABS_FLOOR = 1.2 * SAMPLE_DSTRAIN
# Guo trilinear peak moment: relative tolerance about M^p.
TOL_GUO_REL = 0.06


# ── geometry / config writers ────────────────────────────────────────────────
def _straight_fiber_csv(path, n, spacing):
    with open(path, "w") as f:
        f.write("# straight fiber along x: x, y, z\n")
        for i in range(n):
            f.write(f"{i * spacing:.10e}, 0.0, 0.0\n")


CB_CONFIG = """# Crack-band axial breakage, {n} beads, L_bond = {lbond:.3e} m
# Predicted first-break global strain = {pred:.5f}  (crack-band closed form)
[comm]
processors_x = 1
processors_y = 1
processors_z = 1
[domain]
x_low  = {xlo:.6e}
x_high = {xhi:.6e}
y_low  = -0.0003
y_high =  0.0003
z_low  = -0.0003
z_high =  0.0003
boundary_x = "fixed"
boundary_y = "fixed"
boundary_z = "fixed"
[neighbor]
skin_fraction = 1.2
bin_size = {binsz:.6e}
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
bond_tolerance = 1.05
bond_radius_ratio = 1.0
youngs_modulus = {e:.6e}
shear_modulus  = {g:.6e}
beta_normal  = {beta}
beta_shear   = {beta}
beta_twist   = {beta}
beta_bending = {beta}
[bonds.plasticity.axial]
kind = "piecewise"
breakpoint_strains = [{bp0}, {bp1}]
slope_multipliers  = [{sl0}, {sl1}]
[bonds.breakage]
kind = "axial_strain"
tensile = {{ kind = "crack_band", value_ref = {vref}, l_ref = {lref:.6e}, eps_yield = {epsy} }}
[[group]]
name = "left_grip"
region = {{ type = "block", min = [{lgx0:.6e}, -0.0003, -0.0003], max = [{lgx1:.6e}, 0.0003, 0.0003] }}
dynamic = false
[[group]]
name = "right_grip"
region = {{ type = "block", min = [{rgx0:.6e}, -0.0003, -0.0003], max = [{rgx1:.6e}, 0.0003, 0.0003] }}
dynamic = false
[[move_linear]]
group = "left_grip"
vx = {vneg:.6e}
[[move_linear]]
group = "right_grip"
vx = {vpos:.6e}
[output]
dir = "{outdir}"
[run]
steps = {steps}
thermo = {thermo}
dt = {dt:.6e}
"""


def cb_case_dir(n):
    return os.path.join(SWEEP_DIR, f"crackband_N{n}")


def _write_cb_case(n):
    cdir = cb_case_dir(n)
    os.makedirs(cdir, exist_ok=True)
    spacing = LF / (n - 1)
    radius = spacing / 2.0            # touching beads
    l_bond = spacing
    csv_path = os.path.join(cdir, "fiber.csv")
    _straight_fiber_csv(csv_path, n, spacing)
    # Grip boxes grab exactly the end bead at each side (half-width < spacing).
    hw = 0.4 * spacing
    cfg = CB_CONFIG.format(
        n=n, lbond=l_bond, pred=cb_predicted_break_strain(l_bond),
        xlo=-0.3 * LF, xhi=1.3 * LF, binsz=max(3.0 * spacing, 5.0e-5),
        e=E_AX, nu=POISSON_AX, g=G_AX, rho=RHO_AX, beta=BETA_AX,
        radius=radius, csv=csv_path,
        bp0=AX_BREAKPOINTS[0], bp1=AX_BREAKPOINTS[1], sl0=AX_SLOPES[0], sl1=AX_SLOPES[1],
        vref=CB_VALUE_REF, lref=CB_L_REF, epsy=CB_EPS_YIELD,
        lgx0=-hw, lgx1=hw, rgx0=LF - hw, rgx1=LF + hw,
        vneg=-GRIP_V, vpos=GRIP_V,
        outdir=cdir, steps=STEPS_AX, thermo=STEPS_AX, dt=DT_AX,
    )
    with open(os.path.join(cdir, "config.toml"), "w") as f:
        f.write(cfg)
    return cdir


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


# ── generate ─────────────────────────────────────────────────────────────────
def generate():
    n = 0
    for nb in N_LIST:
        _write_cb_case(nb)
        n += 1
    for s0 in SIGMA0_LIST:
        _write_guo_case(s0)
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


def _measure_crackband(rows):
    """Return (measured_break_strain, plastic_before_break, prebreak_strain)."""
    if not rows:
        return None
    length0 = _fval(rows[0], "length0")
    if not (length0 > 0):
        return None
    gstrain = [(_fval(r, "length_global") / length0) - 1.0 for r in rows]
    broken = [int(_fval(r, "bonds_broken")) for r in rows]
    idx = next((i for i, b in enumerate(broken) if b > 0), None)
    if idx is None or idx == 0:
        return None
    # Threshold crossing is bracketed by [prev sample, first-break sample].
    prev = gstrain[idx - 1]
    meas = 0.5 * (prev + gstrain[idx])
    plastic = _fval(rows[idx - 1], "eps_p_axial_mid") > 0.0
    return meas, plastic, prev


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

    # GROUP A — crack-band mesh sweep.
    cb_rows = []
    print("\n[Group A] crack-band axial breakage, mesh refinement")
    for nb in N_LIST:
        cdir = cb_case_dir(nb)
        if not os.path.isfile(os.path.join(cdir, "config.toml")):
            print(f"  N={nb}: missing config — run 'generate' first.")
            continue
        l_bond = LF / (nb - 1)
        print(f"  N={nb:<3} L_bond={l_bond*1e6:6.1f} um", end="  ", flush=True)
        data = _run(cdir)
        m = _measure_crackband(data) if data is not None else None
        if m is None:
            print("RUN/BREAK FAILED")
            continue
        meas, plastic, prev = m
        pred = cb_predicted_break_strain(l_bond)
        cb_rows.append({"N": nb, "L_bond": l_bond, "eps_break_pred": pred,
                        "eps_break_meas": meas, "prebreak_strain": prev,
                        "plastic_before_break": int(plastic)})
        print(f"break_meas={meas:.4f}  pred={pred:.4f}  plastic={plastic}")

    # GROUP B — Guo trilinear bending, M^p plateau.
    guo_rows = []
    print("\n[Group B] Guo 2018 trilinear bending, fully-plastic moment cap")
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

    if not cb_rows and not guo_rows:
        print("\nERROR: no results collected.")
        sys.exit(1)
    if cb_rows:
        _write_csv(CB_CSV, ["N", "L_bond", "eps_break_pred", "eps_break_meas",
                            "prebreak_strain", "plastic_before_break"], cb_rows)
    if guo_rows:
        _write_csv(GUO_CSV, ["sigma_0", "Mp_theory", "peak_moment"], guo_rows)
    print(f"\nWrote {CB_CSV} ({len(cb_rows)} rows), {GUO_CSV} ({len(guo_rows)} rows)")


# ── graph (validate + plot) ──────────────────────────────────────────────────
def _load(path):
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        return [{k: (float(v) if k not in () else v) for k, v in r.items()}
                for r in csv.DictReader(f)]


def validate(cb_rows, guo_rows):
    ok = True
    print("\n=== Bond-breakage / plasticity sweep validation ===")

    # ── Group A gate: crack-band closed form + plastic-before-break ──
    print("\n[A] Crack-band breakage vs closed form "
          f"eps_yield + (value_ref - eps_yield)*l_ref/L_bond  (rel tol {TOL_CB_REL:.0%}, "
          f"abs floor {TOL_CB_ABS_FLOOR:.4f})")
    print(f"  {'N':>4}{'L_bond(um)':>12}{'pred':>9}{'meas':>9}{'rel.err':>9}  plastic  note")
    if len(cb_rows) != len(N_LIST):
        print(f"  MISSING CASES: got {len(cb_rows)}/{len(N_LIST)} crack-band runs")
        ok = False
    prev_meas = None
    for r in sorted(cb_rows, key=lambda x: x["N"]):
        pred = r["eps_break_pred"]
        meas = r["eps_break_meas"]
        rel = abs(meas - pred) / pred
        tol = max(TOL_CB_REL * pred, TOL_CB_ABS_FLOOR)
        note = ""
        if abs(meas - pred) > tol:
            note = "BREAK-STRAIN MISMATCH"
            ok = False
        if not int(r["plastic_before_break"]):
            note = (note + "; " if note else "") + "NO PLASTIC YIELD BEFORE BREAK"
            ok = False
        # Monotonicity: finer mesh (smaller L_bond, larger N) must break later.
        if prev_meas is not None and meas <= prev_meas:
            note = (note + "; " if note else "") + "NOT MONOTONE"
            ok = False
        prev_meas = meas
        print(f"  {int(r['N']):>4}{r['L_bond']*1e6:>12.1f}{pred:>9.4f}{meas:>9.4f}"
              f"{rel:>8.1%}  {'yes' if int(r['plastic_before_break']) else 'NO ':>7}  {note}")

    # ── Group B gate: Guo trilinear M^p = (4/3) sigma_0 r_b^3 ──
    print(f"\n[B] Guo trilinear peak moment vs M^p = (4/3) sigma_0 r_b^3  (rel tol {TOL_GUO_REL:.0%})")
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

    print("\nRESULT:", "ALL CHECKS PASSED" if ok else "CHECKS FAILED")
    return ok


def plot(cb_rows, guo_rows):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"(matplotlib unavailable, skipping plots: {e})")
        return
    os.makedirs(PLOT_DIR, exist_ok=True)
    plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150, "font.size": 11})

    # Fig 1: crack-band break strain vs 1/L_bond.
    if cb_rows:
        rs = sorted(cb_rows, key=lambda x: x["L_bond"])
        inv_l = [1.0 / r["L_bond"] for r in rs]
        meas = [r["eps_break_meas"] for r in rs]
        xline = [0.0, max(inv_l) * 1.05]
        yline = [CB_EPS_YIELD + (CB_VALUE_REF - CB_EPS_YIELD) * CB_L_REF * x for x in xline]
        fig, ax = plt.subplots(figsize=(6.2, 4.5))
        ax.plot(xline, yline, "k--",
                label=r"crack-band: $\varepsilon_y+(\varepsilon_{ref}-\varepsilon_y)\,l_{ref}/L$")
        ax.plot(inv_l, meas, "o", ms=8, label="DIRT (measured first-break)")
        ax.set_xlabel(r"$1/L_{bond}$  (m$^{-1}$)")
        ax.set_ylabel(r"global break strain $\varepsilon_{break}$")
        ax.set_title("Crack-band regularized breaking strain vs mesh refinement")
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(PLOT_DIR, "crackband_break_strain.png"))
        plt.close(fig)

    # Fig 2: Guo peak moment vs sigma_0.
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

    print(f"Figures -> {PLOT_DIR}/")


def graph():
    cb_rows = _load(CB_CSV)
    guo_rows = _load(GUO_CSV)
    if not cb_rows and not guo_rows:
        print(f"No results in {DATA_DIR} — run 'start' first.")
        return False
    ok = validate(cb_rows, guo_rows)
    plot(cb_rows, guo_rows)
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
