#!/usr/bin/env python3
"""
Polydisperse / multi-material mixing benchmark driver.

Every case is a single binary collision between two spheres that are (in general)
of unequal radius and/or different material. The benchmark checks that DIRT's
per-pair mixing rules enter the contact physics with the correct values:

  * reduced radius        R*   = r1·r2 / (r1 + r2)
  * effective modulus     E*   = 1 / ((1−ν1²)/E1 + (1−ν2²)/E2)   (= `e_eff_ij`)
  * restitution mixing    e_ij = √(e1·e2)                        (geometric mean)
  * sliding friction      μ_ij = √(μ1·μ2)                        (= `friction_ij`)

Two families of cases:

  HEAD-ON (free–free).  Two free spheres collide along the line of centres.
    * Elastic (e=1, μ=0): the peak overlap and contact duration follow the
      undamped Hertz collision *exactly*, evaluated with the mixed R*, E* and the
      reduced mass m*. Because those formulas depend on R* and E* only through
      the mixing rules, matching them for unequal-radius / cross-material pairs
      pins `r_eff = r1 r2/(r1+r2)` and `e_eff_ij`.
        δ_max = (15 m* v0² / (16 E* √R*))^(2/5)
        t_c   = 2.868 (m*² / (R* E*² v0))^(1/5)
      Reference: K.L. Johnson, *Contact Mechanics*, Cambridge Univ. Press, 1985.
    * Inelastic (cross restitution): the realized normal COR equals the geometric
      mean e_ij = √(e1 e2) — DIRT maps input restitution to realized COR via
      `hertz_beta_for_cor`, so this is a direct check of restitution mixing.

  OBLIQUE (projectile onto a FROZEN target).  The projectile strikes an immovable
    target sphere with a large tangential velocity, chosen in the *gross-sliding*
    regime (v_t / v_n well above the (7/2)·μ_ij·(1+e) stick threshold). There the
    tangential contact force sits on the Coulomb cap for the whole contact, so the
    ratio of tangential to normal impulse delivered to the projectile equals the
    pair friction:
        |J_t| / |J_n| = μ_ij = √(μ1·μ2).
    This isolates `friction_ij` and, crucially, distinguishes the geometric mean
    from the arithmetic mean (reported alongside).

Commands (run from anywhere):
    python3 examples/bench_polydisperse_mixing/sweep.py generate
    python3 examples/bench_polydisperse_mixing/sweep.py start
    python3 examples/bench_polydisperse_mixing/sweep.py graph
    python3 examples/bench_polydisperse_mixing/sweep.py            # all three

Outputs:
    sweep/<case>/config.toml         DIRT configs                    (gitignored)
    sweep/<case>/data/...csv         per-case measured results       (gitignored)
    data/sweep_results.csv           collected results               (gitignored)
    plots/mixing_validation.png      summary figure                  (tracked)
"""

import os
import sys
import csv
import math
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_polydisperse_mixing"

SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
SWEEP_CSV = os.path.join(DATA_DIR, "sweep_results.csv")

DENSITY = 2500.0        # kg/m³ — shared by every particle (mass set by radius)

# Two reference materials used by the head-on R*/E* family.
#   name: (Young's modulus E [Pa], Poisson ratio ν)
MAT_SOFT = (8.7e9, 0.30)     # e.g. soda-lime-ish soft glass
MAT_STIFF = (70.0e9, 0.22)   # stiff glass

# ── Tolerances ───────────────────────────────────────────────────────────────
OVERLAP_TOL = 0.01       # elastic peak overlap vs Hertz theory (1%)
DURATION_TOL = 0.03      # elastic contact duration vs Hertz theory (3%)
COR_ELASTIC_TOL = 0.003  # measured COR vs 1.0 for elastic head-on
COR_REF_TOL = 0.01       # cross-restitution COR vs same-e_ij reference COR
FRICTION_TOL = 0.05      # measured |Jt|/|Jn| vs μ_ij = √(μ1 μ2) (5%)

# Oblique impact: projectile velocity (gross-sliding regime, v_t/v_n = 3).
OBLIQUE_VX = 6.0
OBLIQUE_VZ = -2.0
OBLIQUE_OFFX = -0.015    # start x-offset so impact normal ≈ +z (radii sum = 0.01)


def estar(m1, m2):
    (e1, n1), (e2, n2) = m1, m2
    return 1.0 / ((1.0 - n1 * n1) / e1 + (1.0 - n2 * n2) / e2)


# ── Case matrix ──────────────────────────────────────────────────────────────
# Each case is a dict. Common keys: id, kind ("headon" | "oblique").
# Materials are (E, nu, restitution, friction) per side.
def _headon(cid, mat1, mat2, r1, r2, e1=1.0, e2=1.0, **flags):
    (E1, n1), (E2, n2) = mat1, mat2
    c = {
        "id": cid, "kind": "headon",
        "m1": (E1, n1, e1, 0.0), "m2": (E2, n2, e2, 0.0),
        "r1": r1, "r2": r2,
    }
    c.update(flags)  # cor_mix / cor_ref / ref_id
    return c


def _oblique(cid, mu1, mu2, r1, r2, e=0.9, mat=MAT_STIFF):
    (E, n) = mat
    return {
        "id": cid, "kind": "oblique",
        "m1": (E, n, e, mu1), "m2": (E, n, e, mu2),
        "r1": r1, "r2": r2,
    }


def build_cases():
    cases = []
    # Isolate R* (same material, unequal radii; radii always sum to 0.01 m).
    cases.append(_headon("N_R_equal", MAT_SOFT, MAT_SOFT, 0.005, 0.005))
    cases.append(_headon("N_R_46", MAT_SOFT, MAT_SOFT, 0.004, 0.006))
    cases.append(_headon("N_R_37", MAT_SOFT, MAT_SOFT, 0.003, 0.007))
    cases.append(_headon("N_R_28", MAT_SOFT, MAT_SOFT, 0.002, 0.008))
    # Isolate E* (equal radii, cross material).
    cases.append(_headon("N_E_AB", MAT_SOFT, MAT_STIFF, 0.005, 0.005))
    cases.append(_headon("N_E_BB", MAT_STIFF, MAT_STIFF, 0.005, 0.005))
    # Combined R* + E*.
    cases.append(_headon("N_RE_AB46", MAT_SOFT, MAT_STIFF, 0.004, 0.006))
    cases.append(_headon("N_RE_BA37", MAT_STIFF, MAT_SOFT, 0.003, 0.007))
    # Restitution mixing (cross restitution, equal radii). The realized COR carries
    # a KNOWN velocity/viscoelastic offset from the nominal input (see
    # bench_hertz_rebound), so we do NOT compare it to √(e1 e2) directly — that
    # would fold the calibration offset into a mixing check. Instead each cross pair
    # (e1,e2) is paired with a same-material REFERENCE at e_ref = √(e1 e2): equal
    # e_ij ⇒ equal beta_ij ⇒ identical realized COR. Matching the two isolates the
    # geometric-mean mixing rule and still fails if the code mixed arithmetically.
    for cid, ea, eb in [("COR_A", 1.0, 0.64), ("COR_B", 0.9, 0.5)]:
        e_ref = math.sqrt(ea * eb)
        cases.append(_headon(cid, MAT_STIFF, MAT_STIFF, 0.005, 0.005,
                             e1=ea, e2=eb, cor_mix=True, ref_id=cid + "_ref"))
        cases.append(_headon(cid + "_ref", MAT_STIFF, MAT_STIFF, 0.005, 0.005,
                             e1=e_ref, e2=e_ref, cor_ref=True))
    # Friction mixing (frozen-target oblique gross sliding): |Jt|/|Jn| = √(μ1 μ2).
    cases.append(_oblique("F_02_08", 0.2, 0.8, 0.005, 0.005))     # geo 0.400 vs arith 0.500
    cases.append(_oblique("F_01_09", 0.1, 0.9, 0.005, 0.005))     # geo 0.300 vs arith 0.500
    cases.append(_oblique("F_04_04", 0.4, 0.4, 0.005, 0.005))     # same-material sanity: 0.400
    cases.append(_oblique("F_016_064", 0.16, 0.64, 0.005, 0.005))  # geo 0.320 vs arith 0.400
    cases.append(_oblique("F_02_08_46", 0.2, 0.8, 0.004, 0.006))  # R* + friction combined
    return cases


CSV_FIELDS = [
    "id", "kind", "E1", "nu1", "e1", "mu1", "E2", "nu2", "e2", "mu2",
    "v_n_impact", "v_n_rebound", "cor", "contact_time", "max_overlap",
    "jt", "jn", "m1", "m2", "r1", "r2",
]


# ── config templates ─────────────────────────────────────────────────────────
HEADON_TOML = """\
# Auto-generated head-on collision — case {cid}
[comm]
processors_x = 1
processors_y = 1
processors_z = 1

[domain]
x_low = 0.0
x_high = 0.06
y_low = -0.02
y_high = 0.02
z_low = -0.02
z_high = 0.02
boundary_x = "fixed"
boundary_y = "fixed"
boundary_z = "fixed"

[neighbor]
skin_fraction = 1.2
bin_size = 0.017

[dem]
contact_model = "hertz"

[[dem.materials]]
name = "m1"
youngs_mod = {E1}
poisson_ratio = {nu1}
restitution = {e1}
friction = {mu1}

[[dem.materials]]
name = "m2"
youngs_mod = {E2}
poisson_ratio = {nu2}
restitution = {e2}
friction = {mu2}

[[particles.insert]]
source = "file"
format = "csv"
file = "{csv}"
material = "m1"
density = {rho}
columns = {{ x = 0, y = 1, z = 2, vx = 3, vy = 4, vz = 5, radius = 6, atom_type = 7 }}
type_map = {{ 1 = "m1", 2 = "m2" }}

[output]
dir = "{out}"

[[run]]
name = "collide"
dt = 2.0e-8
steps = 40000
thermo = 20000
"""

OBLIQUE_TOML = """\
# Auto-generated oblique (frozen-target) collision — case {cid}
[comm]
processors_x = 1
processors_y = 1
processors_z = 1

[domain]
x_low = -0.2
x_high = 0.2
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
name = "m1"
youngs_mod = {E1}
poisson_ratio = {nu1}
restitution = {e1}
friction = {mu1}

[[dem.materials]]
name = "m2"
youngs_mod = {E2}
poisson_ratio = {nu2}
restitution = {e2}
friction = {mu2}

# Projectile FIRST (tag 0 → "particle 1"): material m1, radius r1.
[[particles.insert]]
material = "m1"
count = 1
radius = {r1}
density = {rho}
velocity_x = {vx}
velocity_z = {vz}
region = {{ type = "block", min = [{offx_lo}, -1.0e-6, 0.064999], max = [{offx_hi}, 1.0e-6, 0.065001] }}

# Frozen target: material m2, radius r2, centred at (0, 0, 0.05).
[[particles.insert]]
material = "m2"
count = 1
radius = {r2}
density = {rho}
region = {{ type = "block", min = [-1.0e-6, -1.0e-6, 0.049999], max = [1.0e-6, 1.0e-6, 0.050001] }}

[[group]]
name = "target"
region = {{ type = "block", min = [-0.003, -0.003, 0.043], max = [0.003, 0.003, 0.057] }}

[[freeze]]
group = "target"

[output]
dir = "{out}"

[[run]]
name = "collide"
dt = 2.0e-8
steps = 200000
thermo = 100000
"""


def case_dir(cid):
    return os.path.join(SWEEP_DIR, cid)


def generate():
    n = 0
    for c in build_cases():
        cdir = case_dir(c["id"])
        os.makedirs(cdir, exist_ok=True)
        (E1, n1, e1, mu1) = c["m1"]
        (E2, n2, e2, mu2) = c["m2"]
        if c["kind"] == "headon":
            csv_path = os.path.join(cdir, "particles.csv")
            x1 = 0.010
            x2 = x1 + c["r1"] + c["r2"] + 2.0e-4   # small pre-contact gap
            with open(csv_path, "w") as f:
                f.write("x,y,z,vx,vy,vz,radius,atype\n")
                f.write(f"{x1:.6f},0.0,0.0,1.0,0.0,0.0,{c['r1']:.6f},1\n")
                f.write(f"{x2:.6f},0.0,0.0,0.0,0.0,0.0,{c['r2']:.6f},2\n")
            with open(os.path.join(cdir, "config.toml"), "w") as f:
                f.write(HEADON_TOML.format(
                    cid=c["id"], E1=E1, nu1=n1, e1=e1, mu1=mu1,
                    E2=E2, nu2=n2, e2=e2, mu2=mu2,
                    csv=csv_path, rho=DENSITY, out=cdir,
                ))
        else:  # oblique
            with open(os.path.join(cdir, "config.toml"), "w") as f:
                f.write(OBLIQUE_TOML.format(
                    cid=c["id"], E1=E1, nu1=n1, e1=e1, mu1=mu1,
                    E2=E2, nu2=n2, e2=e2, mu2=mu2,
                    r1=c["r1"], r2=c["r2"], rho=DENSITY,
                    vx=OBLIQUE_VX, vz=OBLIQUE_VZ,
                    offx_lo=OBLIQUE_OFFX - 1.0e-6, offx_hi=OBLIQUE_OFFX + 1.0e-6,
                    out=cdir,
                ))
        n += 1
    print(f"Generated {n} case configs under {SWEEP_DIR}")


def start():
    os.makedirs(DATA_DIR, exist_ok=True)
    print(f"Building {EXAMPLE} (release)...", flush=True)
    subprocess.run(
        ["cargo", "build", "--release", "--example", EXAMPLE,
         "--no-default-features", "--features", "precision-double"],
        cwd=REPO_ROOT, check=True,
    )

    cases = build_cases()
    results = []
    for i, c in enumerate(cases, 1):
        cdir = case_dir(c["id"])
        config = os.path.join(cdir, "config.toml")
        if not os.path.isfile(config):
            print(f"  [{i:2d}/{len(cases)}] {c['id']}: missing config — run 'generate' first.")
            continue
        out_csv = os.path.join(cdir, "data", "collision_results.csv")
        if os.path.isfile(out_csv):
            os.remove(out_csv)
        print(f"  [{i:2d}/{len(cases)}] {c['id']:<12} ({c['kind']})", end="  ", flush=True)
        log = os.path.join(cdir, "run.log")
        with open(log, "w") as lf:
            proc = subprocess.run(
                ["cargo", "run", "--release", "--example", EXAMPLE,
                 "--no-default-features", "--features", "precision-double", "--", config],
                cwd=REPO_ROOT, stdout=lf, stderr=subprocess.STDOUT,
            )
        if proc.returncode != 0 or not os.path.isfile(out_csv):
            print(f"FAILED (see {log})")
            continue
        with open(out_csv) as f:
            row = next(csv.DictReader(f))
        (E1, n1, e1, mu1) = c["m1"]
        (E2, n2, e2, mu2) = c["m2"]
        row.update({
            "id": c["id"], "kind": c["kind"],
            "E1": E1, "nu1": n1, "e1": e1, "mu1": mu1,
            "E2": E2, "nu2": n2, "e2": e2, "mu2": mu2,
        })
        results.append(row)
        if c["kind"] == "oblique":
            jt, jn = float(row["jt"]), float(row["jn"])
            print(f"Jt/Jn={jt/jn:.4f}")
        else:
            print(f"COR={float(row['cor']):.4f}")

    if not results:
        print("\nERROR: no results collected.")
        sys.exit(1)
    with open(SWEEP_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()
        for r in results:
            w.writerow({k: r.get(k, "") for k in CSV_FIELDS})
    print(f"\nCollected {len(results)}/{len(cases)} cases -> {SWEEP_CSV}")


# ── theory ───────────────────────────────────────────────────────────────────
def hertz_overlap(mstar, v0, Estar, Rstar):
    return (15.0 * mstar * v0 * v0 / (16.0 * Estar * math.sqrt(Rstar))) ** 0.4


def hertz_duration(mstar, v0, Estar, Rstar):
    return 2.868 * (mstar * mstar / (Rstar * Estar * Estar * v0)) ** 0.2


def load_rows():
    if not os.path.isfile(SWEEP_CSV):
        print(f"ERROR: {SWEEP_CSV} not found. Run 'start' first.")
        sys.exit(1)
    with open(SWEEP_CSV) as f:
        return list(csv.DictReader(f))


def validate(rows):
    print("=" * 74)
    print("Polydisperse / multi-material mixing validation")
    print("=" * 74)
    checks = []  # (label, passed)
    by_id = {r["id"]: r for r in rows}
    cases = build_cases()

    print("\n── Head-on: R* and E* mixing (elastic peak overlap & contact time) ──")
    for c in cases:
        if c["kind"] != "headon" or c.get("cor_mix") or c.get("cor_ref"):
            continue
        r = by_id.get(c["id"])
        if r is None:
            continue
        E1, n1 = float(r["E1"]), float(r["nu1"])
        E2, n2 = float(r["E2"]), float(r["nu2"])
        m1, m2 = float(r["m1"]), float(r["m2"])
        r1, r2 = float(r["r1"]), float(r["r2"])
        v0 = float(r["v_n_impact"])
        cor = float(r["cor"])
        Est = estar((E1, n1), (E2, n2))
        Rst = r1 * r2 / (r1 + r2)
        mst = m1 * m2 / (m1 + m2)
        om, tm = float(r["max_overlap"]), float(r["contact_time"])
        d_th = hertz_overlap(mst, v0, Est, Rst)
        t_th = hertz_duration(mst, v0, Est, Rst)
        oerr, terr = abs(om - d_th) / d_th, abs(tm - t_th) / t_th
        po, pt = oerr <= OVERLAP_TOL, terr <= DURATION_TOL
        pc = abs(cor - 1.0) <= COR_ELASTIC_TOL
        checks += [(f"{c['id']} overlap", po), (f"{c['id']} duration", pt),
                   (f"{c['id']} COR=1", pc)]
        print(f"  {c['id']:<11} R*={Rst*1e3:.3f}mm E*={Est:.3e}  "
              f"δ {om*1e6:6.2f}/{d_th*1e6:6.2f}µm ({oerr*100:4.1f}%) [{ok(po)}]  "
              f"t_c {tm*1e6:6.2f}/{t_th*1e6:6.2f}µs ({terr*100:4.1f}%) [{ok(pt)}]  "
              f"COR {cor:.4f} [{ok(pc)}]")

    print("\n── Restitution mixing: cross (e1,e2) vs same-material ref at e=√(e1 e2) ──")
    for c in cases:
        if not c.get("cor_mix"):
            continue
        r = by_id.get(c["id"])
        ref = by_id.get(c["ref_id"])
        if r is None or ref is None:
            checks.append((f"{c['id']} restitution-mix", False))
            print(f"  {c['id']:<11} MISSING measurement or reference [FAIL]")
            continue
        (_, _, e1, _), (_, _, e2, _) = c["m1"], c["m2"]
        e_ij = math.sqrt(e1 * e2)
        arith = 0.5 * (e1 + e2)
        cor, cref = float(r["cor"]), float(ref["cor"])
        d = abs(cor - cref)
        p = d <= COR_REF_TOL
        # Wrong (arithmetic) mixing would target e_ij=arith; the reference is at the
        # geometric e_ij, so a code using arithmetic mixing would separate the two by
        # ~|geo−arith| ≫ tol. Report that gap so the discrimination is explicit.
        checks.append((f"{c['id']} restitution-mix (geo)", p))
        print(f"  {c['id']:<11} e_ij=√({e1}·{e2})={e_ij:.4f} (arith {arith:.4f}):  "
              f"COR cross {cor:.4f} vs ref {cref:.4f} (Δ={d:.4f}) [{ok(p)}]")

    print("\n── Oblique (frozen target, gross sliding): friction_ij = √(μ1 μ2) ──")
    for c in cases:
        if c["kind"] != "oblique":
            continue
        r = by_id.get(c["id"])
        if r is None:
            continue
        mu1, mu2 = float(r["mu1"]), float(r["mu2"])
        jt, jn = float(r["jt"]), float(r["jn"])
        ratio = jt / jn
        geo = math.sqrt(mu1 * mu2)
        arith = 0.5 * (mu1 + mu2)
        gerr = abs(ratio - geo) / geo
        pg = gerr <= FRICTION_TOL
        # Discrimination: measured must be closer to geometric than arithmetic mean
        # (skip when the two coincide, i.e. μ1 == μ2).
        discriminates = (abs(mu1 - mu2) < 1e-12) or (abs(ratio - geo) < abs(ratio - arith))
        checks.append((f"{c['id']} Jt/Jn=√(μ1 μ2)", pg))
        checks.append((f"{c['id']} geo<arith", discriminates))
        print(f"  {c['id']:<11} Jt/Jn meas {ratio:.4f}  geo √({mu1}·{mu2})={geo:.4f} "
              f"({gerr*100:4.1f}%) [{ok(pg)}]  arith={arith:.4f}  "
              f"closer-to-geo [{ok(discriminates)}]")

    npass = sum(1 for _, p in checks if p)
    ntot = len(checks)
    print("\n" + "-" * 74)
    print(f"Overall: {npass}/{ntot} checks passed")
    allok = npass == ntot
    print("ALL CHECKS PASSED" if allok else f"CHECKS FAILED ({ntot - npass} failed)")
    return allok


def ok(b):
    return "PASS" if b else "FAIL"


def plot(rows):
    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"\n(matplotlib/numpy unavailable, skipped plot: {e})")
        return
    os.makedirs(PLOT_DIR, exist_ok=True)
    plt.rcParams.update({"font.size": 11, "figure.dpi": 150, "savefig.dpi": 150})

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Left: measured vs theory peak overlap (elastic head-on cases).
    xs, ys, labels = [], [], []
    for r in rows:
        if r["kind"] != "headon" or float(r["e1"]) != 1.0 or float(r["e2"]) != 1.0:
            continue
        E1, n1 = float(r["E1"]), float(r["nu1"])
        E2, n2 = float(r["E2"]), float(r["nu2"])
        m1, m2 = float(r["m1"]), float(r["m2"])
        r1, r2 = float(r["r1"]), float(r["r2"])
        v0 = float(r["v_n_impact"])
        Est = estar((E1, n1), (E2, n2))
        Rst = r1 * r2 / (r1 + r2)
        mst = m1 * m2 / (m1 + m2)
        xs.append(hertz_overlap(mst, v0, Est, Rst) * 1e6)
        ys.append(float(r["max_overlap"]) * 1e6)
        labels.append(r["id"])
    if xs:
        lim = [0, max(xs + ys) * 1.15]
        ax1.plot(lim, lim, "k--", lw=1, label="theory = measured")
        ax1.plot(xs, ys, "o", color="#1f77b4", ms=8, label="cases")
        ax1.set_xlabel("Hertz theory peak overlap [µm]  (mixed R*, E*)")
        ax1.set_ylabel("DIRT measured peak overlap [µm]")
        ax1.set_title("Head-on: R* & E* mixing")
        ax1.legend(loc="upper left", fontsize=9)
        ax1.set_aspect("equal")

    # Right: measured Jt/Jn vs geometric-mean friction (oblique cases).
    gm, meas, am, labs = [], [], [], []
    for r in rows:
        if r["kind"] != "oblique":
            continue
        mu1, mu2 = float(r["mu1"]), float(r["mu2"])
        gm.append(math.sqrt(mu1 * mu2))
        am.append(0.5 * (mu1 + mu2))
        meas.append(float(r["jt"]) / float(r["jn"]))
        labs.append(r["id"])
    if gm:
        order = np.argsort(gm)
        gm = np.array(gm)[order]; meas = np.array(meas)[order]
        am = np.array(am)[order]
        idx = np.arange(len(gm))
        ax2.plot(idx, gm, "s-", color="#2ca02c", label="geometric mean √(μ1 μ2)")
        ax2.plot(idx, am, "^--", color="#d62728", alpha=0.6, label="arithmetic mean (½(μ1+μ2))")
        ax2.plot(idx, meas, "o", color="#1f77b4", ms=9, label="DIRT measured Jt/Jn")
        ax2.set_xticks(idx)
        ax2.set_xticklabels([labs[i] for i in order], rotation=30, ha="right", fontsize=8)
        ax2.set_ylabel("tangential/normal impulse ratio")
        ax2.set_title("Oblique gross sliding: friction_ij mixing")
        ax2.legend(loc="best", fontsize=9)

    fig.tight_layout()
    out = os.path.join(PLOT_DIR, "mixing_validation.png")
    fig.savefig(out, bbox_inches="tight")
    import matplotlib.pyplot as plt2  # noqa
    print(f"Saved: {out}")


def graph():
    rows = load_rows()
    allok = validate(rows)
    plot(rows)
    return allok


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
        print("Usage: sweep.py [generate|start|graph]")
        sys.exit(2)


if __name__ == "__main__":
    main()
