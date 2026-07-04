#!/usr/bin/env python3
"""
SDS (spring-dashpot-slider) rolling-resistance benchmark driver.

Two identical spheres are stacked along +z. The lower sphere is a frozen anchor
([[freeze]] — it can neither translate nor spin), and the upper sphere rests on
it under gravity, seated at the static Hertz overlap so the normal contact force
is F_n = m*g. The upper sphere is given a pure ROLLING spin omega = (omega0,0,0)
about a horizontal axis (perpendicular to the contact normal n = +z), so its
projection onto n is zero => this is rolling, not twisting. Sliding friction is
turned OFF (friction = 0) so the tangential slip a horizontal spin would drive
cannot add a couple: the ONLY torque perpendicular to n is the SDS rolling
resistance under test (contact.rs). With F_n = m*g, equal-sphere r_eff = R/2 and
sphere inertia I = (2/5) m R^2, the model is validated in its two regimes:

  * ELASTIC (Coulomb cap disengaged): the rolling "displacement" delta integrates
    the rolling velocity and the couple is tau = -k_r*delta - gamma_r*omega, so
    the spin obeys the EXACT damped linear oscillator

        I*delta'' + gamma_r*delta' + k_r*delta = 0 ,  omega = delta' .

    We compare omega(t) point-by-point against the closed-form solution of that
    ODE (over-damped => (near-)exponential decay set by the larger eigenvalue
    |s1|; under-damped => a decaying oscillation whose spring restoring force
    reverses the spin). PASS requires max|omega_sim - omega_exact| <= tol*omega0.
    As a discriminating control we also report the error of the springless
    (k_r = 0) pure-dashpot curve — the spring must matter, or the test is empty.

  * COULOMB CAP (slider saturated): under a large sustained spin the spring+dashpot
    torque exceeds the cap, so the slider holds tau = tau_max = mu_r*F_n*r_eff and
    the spin decays at the EXACT constant rate

        alpha = domega/dt = -(5/4) * mu_r * g / R

    (same closed form as the constant-rolling / twisting benches). We fit the
    saturated slope past the brief spring wind-up and compare to alpha.

Every case also checks the roll stayed pure: off-axis spin |omega_perp| ~ 0 and no
lateral drift.

Commands (from the dirt repo root):
    python3 examples/bench_sds_rolling/sweep.py generate   # write per-case configs
    python3 examples/bench_sds_rolling/sweep.py start      # build + run all sims -> CSV
    python3 examples/bench_sds_rolling/sweep.py graph      # validate + plot
    python3 examples/bench_sds_rolling/sweep.py            # all three, in order

Outputs:
    sweep/<case>/config.toml       DIRT configs                       (gitignored)
    data/roll_<case>.csv           per-case DIRT time series          (gitignored)
    data/sweep.csv                 per-case validation summary        (gitignored)
    plots/*.png                    final figures                      (tracked)

Reference: the SDS rolling model is Model C ("elastic-plastic spring-dashpot") of
J. Ai, J.-F. Chen, J.M. Rotter, J.Y. Ooi, "Assessment of rolling resistance models
in discrete element simulations", Powder Technology 206 (2011) 269-282; the same
spring-dashpot-slider rolling model is implemented in LAMMPS `pair_granular`
(rolling `sds`, doc/src/pair_granular.rst). Both the elastic damped-oscillator
dynamics and the mu_r*F_n*r_eff Coulomb cap are the model's own defining behaviour,
which is what this benchmark validates.
"""

import os
import sys
import csv
import math
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_sds_rolling"

SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
SWEEP_CSV = os.path.join(DATA_DIR, "sweep.csv")

# ── Material / geometry (fixed across the sweep) ──────────────────────────────
NU = 0.3              # Poisson ratio
E_N = 0.3             # normal restitution (damps the tiny settling transient)
YOUNGS_MOD = 1.0e8    # Pa — soft (keeps dt large; the decay rate is E-independent)
DENSITY = 2500.0      # kg/m^3
RADIUS = 0.005        # m — sphere radius R
GRAVITY = 9.81        # m/s^2
OMEGA0 = 8.0          # rad/s — initial rolling spin (passed to main.rs via SDS_OMEGA0)
DT = 1.0e-5           # s

# Derived constants used by the analytics.
MASS = DENSITY * (4.0 / 3.0) * math.pi * RADIUS**3
INERTIA = 0.4 * MASS * RADIUS**2          # I = (2/5) m R^2
F_N = MASS * GRAVITY                       # seated normal force F_n = m g
R_EFF = RADIUS / 2.0                        # equal spheres

# ── Cases ─────────────────────────────────────────────────────────────────────
# Each case is a dict. kind="elastic" validates the damped-oscillator trajectory;
# kind="cap" validates the saturated Coulomb slope. k_r / gamma_r are the rolling
# spring stiffness (N*m/rad) and dashpot (N*m*s/rad); mu_r the rolling friction.
#
# Elastic cases use a LARGE mu_r so the cap never engages (pure spring-dashpot).
# The over-damped case gives a clean exponential omega decay; the under-damped
# case oscillates (the spring reverses the spin) — an unambiguous spring signature.
# Cap cases use gamma_r = 0 and a stiff spring that saturates within a few steps,
# then the slider holds tau = mu_r*F_n*r_eff; mu_r is swept.
CASES = [
    {"tag": "elastic_overdamped",  "kind": "elastic",
     "k_r": 1.45e-5,  "gamma_r": 2.618e-6, "mu_r": 5.0,  "steps": 40000},
    {"tag": "elastic_underdamped", "kind": "elastic",
     "k_r": 4.7124e-3, "gamma_r": 1.5708e-6, "mu_r": 50.0, "steps": 6000},
    {"tag": "cap_mu_0.05", "kind": "cap",
     "k_r": 1.0e-2, "gamma_r": 0.0, "mu_r": 0.05, "steps": 60000},
    {"tag": "cap_mu_0.10", "kind": "cap",
     "k_r": 1.0e-2, "gamma_r": 0.0, "mu_r": 0.10, "steps": 40000},
    {"tag": "cap_mu_0.20", "kind": "cap",
     "k_r": 1.0e-2, "gamma_r": 0.0, "mu_r": 0.20, "steps": 20000},
]

# ── Theory ────────────────────────────────────────────────────────────────────
def cap_alpha(mu_r):
    """Exact saturated (Coulomb-cap) rolling spin-down rate:
    alpha = tau_max/I = mu_r*F_n*r_eff / ((2/5) m R^2) = (5/4) mu_r g / R."""
    return (5.0 / 4.0) * mu_r * GRAVITY / RADIUS


def elastic_omega(t, k_r, gamma_r, w0):
    """Exact omega(t) of I*delta'' + gamma_r*delta' + k_r*delta = 0 with
    delta(0)=0, omega(0)=w0 (omega=delta'). Handles over/critical/under-damped."""
    I = INERTIA
    disc = gamma_r * gamma_r - 4.0 * k_r * I
    if disc > 0.0:                     # over-damped: two real roots
        rt = math.sqrt(disc)
        s1 = (-gamma_r - rt) / (2.0 * I)
        s2 = (-gamma_r + rt) / (2.0 * I)
        return w0 * (s1 * math.exp(s1 * t) - s2 * math.exp(s2 * t)) / (s1 - s2)
    elif disc < 0.0:                   # under-damped: decaying oscillation
        wn = math.sqrt(k_r / I)
        zeta = gamma_r / (2.0 * math.sqrt(k_r * I))
        wd = wn * math.sqrt(1.0 - zeta * zeta)
        return w0 * math.exp(-zeta * wn * t) * (
            math.cos(wd * t) - (zeta * wn / wd) * math.sin(wd * t))
    else:                              # critically damped
        wn = math.sqrt(k_r / I)
        return w0 * math.exp(-wn * t) * (1.0 - wn * t)


def elastic_springless(t, gamma_r, w0):
    """Discriminating control: k_r = 0 pure-dashpot decay, omega = w0 exp(-(gamma/I) t)."""
    return w0 * math.exp(-(gamma_r / INERTIA) * t)


def elastic_eigen(k_r, gamma_r):
    """Return (label, dominant_rate, zeta, wn) describing the elastic regime."""
    I = INERTIA
    wn = math.sqrt(k_r / I)
    zeta = gamma_r / (2.0 * math.sqrt(k_r * I))
    disc = gamma_r * gamma_r - 4.0 * k_r * I
    if disc > 0.0:
        rt = math.sqrt(disc)
        s1 = (gamma_r + rt) / (2.0 * I)   # |larger eigenvalue| (dominant decay)
        return "overdamped", s1, zeta, wn
    elif disc < 0.0:
        return "underdamped", zeta * wn, zeta, wn
    return "critical", wn, zeta, wn


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
rolling_model = "sds"
[[dem.materials]]
name = "grain"
youngs_mod = {youngs:.6e}
poisson_ratio = {nu}
restitution = {e_n}
friction = 0.0
rolling_friction = {mu_r:.6e}
rolling_stiffness = {k_r:.6e}
rolling_damping = {gamma_r:.6e}
[[particles.insert]]
material = "grain"
count = 1
radius = {radius}
density = {density}
region = {{ type = "block", min = [-1.0e-6, -1.0e-6, 0.0099990], max = [1.0e-6, 1.0e-6, 0.0100010] }}
[[particles.insert]]
material = "grain"
count = 1
radius = {radius}
density = {density}
region = {{ type = "block", min = [-1.0e-6, -1.0e-6, 0.0199970], max = [1.0e-6, 1.0e-6, 0.0199990] }}
[[group]]
name = "anchor"
region = {{ type = "block", min = [-0.002, -0.002, 0.008], max = [0.002, 0.002, 0.012] }}
[[freeze]]
group = "anchor"
[output]
dir = "{outdir}"
[run]
steps = {steps}
thermo = {steps}
dt = {dt:.6e}
"""


# ── helpers ───────────────────────────────────────────────────────────────────
def case_dir(c):
    return os.path.join(SWEEP_DIR, c["tag"])


def roll_csv(c):
    return os.path.join(DATA_DIR, f"roll_{c['tag']}.csv")


def _dirt_config(c, outdir):
    return TOML_TEMPLATE.format(
        g=GRAVITY, youngs=YOUNGS_MOD, nu=NU, e_n=E_N, mu_r=c["mu_r"],
        k_r=c["k_r"], gamma_r=c["gamma_r"], radius=RADIUS, density=DENSITY,
        outdir=outdir, steps=c["steps"], dt=DT,
    )


# ── generate ──────────────────────────────────────────────────────────────────
def generate():
    for c in CASES:
        cdir = case_dir(c)
        os.makedirs(cdir, exist_ok=True)
        with open(os.path.join(cdir, "config.toml"), "w") as f:
            f.write(_dirt_config(c, cdir))
    print(f"Generated {len(CASES)} DIRT sweep configs under {SWEEP_DIR}")


# ── start ─────────────────────────────────────────────────────────────────────
SWEEP_FIELDS = ["tag", "kind", "mu_r", "k_r", "gamma_r", "metric", "value",
                "target", "rel_err", "null_err", "max_perp", "max_drift", "npts"]


def _read_timeseries(path):
    """Read t,omega_roll,omega_perp,drift from an sds_rolling_results.csv."""
    rows = []
    with open(path) as f:
        for parts in csv.reader(f):
            if not parts or parts[0].startswith("#") or parts[0] == "t":
                continue
            if len(parts) >= 4:
                rows.append(tuple(float(p) for p in parts[:4]))
    return rows  # (t, omega_roll, omega_perp, drift)


def _run_dirt(c):
    cdir = case_dir(c)
    config = os.path.join(cdir, "config.toml")
    res = os.path.join(cdir, "data", "sds_rolling_results.csv")
    if os.path.exists(res):
        os.remove(res)
    env = dict(os.environ, SDS_OMEGA0=f"{OMEGA0}")
    log = os.path.join(cdir, "run.log")
    with open(log, "w") as lf:
        proc = subprocess.run(
            ["cargo", "run", "--release", "--example", EXAMPLE,
             "--no-default-features", "--features", "precision-double", "--", config],
            cwd=REPO_ROOT, stdout=lf, stderr=subprocess.STDOUT, env=env,
        )
    if proc.returncode != 0 or not os.path.isfile(res):
        return None
    rows = _read_timeseries(res)
    with open(roll_csv(c), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t", "omega_roll", "omega_perp", "drift"])
        w.writerows(rows)
    return rows


def _analyze_elastic(c, rows):
    """Compare omega(t) to the exact damped-oscillator solution; report the
    max |omega_sim - omega_exact| / omega0 (metric) and the springless-control
    error (null_err). Returns a summary dict or None."""
    if len(rows) < 20:
        return None
    err = 0.0
    null = 0.0
    for t, w, _, _ in rows:
        err = max(err, abs(w - elastic_omega(t, c["k_r"], c["gamma_r"], OMEGA0)))
        null = max(null, abs(w - elastic_springless(t, c["gamma_r"], OMEGA0)))
    err /= OMEGA0
    null /= OMEGA0
    max_perp = max(r[2] for r in rows) / OMEGA0
    max_drift = max(r[3] for r in rows)
    _, rate, _, _ = elastic_eigen(c["k_r"], c["gamma_r"])
    return {"metric": "max_abs_err_over_w0", "value": err, "target": rate,
            "rel_err": err, "null_err": null, "max_perp": max_perp,
            "max_drift": max_drift, "npts": len(rows)}


def _analyze_cap(c, rows):
    """Fit the saturated spin-down slope past the wind-up and compare to
    alpha = (5/4) mu_r g / R. Returns a summary dict or None."""
    lo, hi = 0.15 * OMEGA0, 0.85 * OMEGA0
    prefix = []
    for r in rows:               # leading monotone-decay prefix
        if r[1] < lo:
            break
        prefix.append(r)
    win = [r for r in prefix if lo <= r[1] <= hi]
    if len(win) < 10:
        return None
    ts = [r[0] for r in win]
    ws = [r[1] for r in win]
    n = len(ts)
    tbar = sum(ts) / n
    wbar = sum(ws) / n
    sxx = sum((t - tbar) ** 2 for t in ts)
    sxy = sum((t - tbar) * (w - wbar) for t, w in zip(ts, ws))
    a_fit = -(sxy / sxx) if sxx > 0 else 0.0
    a_pred = cap_alpha(c["mu_r"])
    rel = abs(a_fit - a_pred) / a_pred
    max_perp = max(r[2] for r in win) / OMEGA0
    max_drift = max(r[3] for r in win)
    return {"metric": "cap_slope_alpha", "value": a_fit, "target": a_pred,
            "rel_err": rel, "null_err": 0.0, "max_perp": max_perp,
            "max_drift": max_drift, "npts": n}


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

    if os.path.exists(SWEEP_CSV):
        os.remove(SWEEP_CSV)

    out_rows = []
    n = len(CASES)
    for i, c in enumerate(CASES, 1):
        cdir = case_dir(c)
        if not os.path.isfile(os.path.join(cdir, "config.toml")):
            print(f"  [{i}/{n}] missing config for {c['tag']} — run 'generate' first.")
            continue
        print(f"  [{i}/{n}] {c['kind']:<8} {c['tag']:<20}", end="  ", flush=True)
        rows = _run_dirt(c)
        if rows is None:
            print("DIRT FAILED")
            continue
        summary = (_analyze_elastic(c, rows) if c["kind"] == "elastic"
                   else _analyze_cap(c, rows))
        if summary is None:
            print("no fit")
            continue
        rec = {"tag": c["tag"], "kind": c["kind"], "mu_r": c["mu_r"],
               "k_r": c["k_r"], "gamma_r": c["gamma_r"], **summary}
        out_rows.append(rec)
        if c["kind"] == "elastic":
            print(f"err={summary['rel_err']*100:.2f}% (null={summary['null_err']*100:.2f}%) "
                  f"perp={summary['max_perp']:.1e} drift={summary['max_drift']:.1e}")
        else:
            print(f"a_fit={summary['value']:.2f} a_pred={summary['target']:.2f} "
                  f"rel={summary['rel_err']*100:.2f}% perp={summary['max_perp']:.1e}")

    if not out_rows:
        print("\nERROR: no DIRT results collected.")
        sys.exit(1)
    _write_csv(SWEEP_CSV, SWEEP_FIELDS, out_rows)
    print(f"\nDIRT: {len(out_rows)}/{n} cases -> {SWEEP_CSV}")


# ── graph (validate + plot) ───────────────────────────────────────────────────
ELASTIC_TOL = 0.015   # elastic omega(t) must match exact analytical within 1.5% of omega0
NULL_MARGIN = 3.0     # ... and the springless control must be >= 3x worse (spring must matter)
SLOPE_TOL = 0.03      # 3% relative error on the saturated Coulomb slope (theory is exact)
PERP_TOL = 1.0e-3     # |omega_perp| must stay < 0.1% of omega0 (pure roll)
DRIFT_TOL = 1.0e-5    # lateral drift must stay < 10 um (no translation)


def _load(path):
    if not os.path.isfile(path):
        return []
    out = []
    with open(path) as f:
        for r in csv.DictReader(f):
            row = dict(r)
            for k in ("mu_r", "k_r", "gamma_r", "value", "target", "rel_err",
                      "null_err", "max_perp", "max_drift", "npts"):
                row[k] = float(row[k]) if row[k] not in ("", None) else 0.0
            out.append(row)
    return out


def validate(rows):
    print("\n=== SDS rolling-resistance validation ===")
    print(f"  R={RADIUS} m  equal spheres (r_eff=R/2)  g={GRAVITY}  omega0={OMEGA0}")
    print(f"  elastic: omega(t) vs exact  I d2delta + gamma d(delta) + k delta = 0")
    print(f"  cap:     alpha = (5/4) mu_r g / R   (saturated slider)")
    print(f"  {'case':>20}{'kind':>10}{'metric':>10}{'target':>11}"
          f"{'err/rel':>10}{'null':>9}{'perp':>10}{'drift':>10}  note")
    ok = True
    for r in sorted(rows, key=lambda x: (x["kind"], x["tag"])):
        note = ""
        if r["kind"] == "elastic":
            if r["rel_err"] > ELASTIC_TOL:
                note = "TRAJ MISMATCH"; ok = False
            if r["null_err"] < NULL_MARGIN * max(r["rel_err"], 1e-9):
                note = (note + " SPRING-NULL").strip(); ok = False
            metric = f"{r['rel_err']*100:.2f}%"
            tgt = f"|s1|={r['target']:.1f}"
        else:
            if r["rel_err"] > SLOPE_TOL:
                note = "SLOPE MISMATCH"; ok = False
            metric = f"{r['value']:.1f}"
            tgt = f"{r['target']:.1f}"
        if r["max_perp"] > PERP_TOL:
            note = (note + " OFF-AXIS").strip(); ok = False
        if r["max_drift"] > DRIFT_TOL:
            note = (note + " DRIFTED").strip(); ok = False
        print(f"  {r['tag']:>20}{r['kind']:>10}{metric:>10}{tgt:>11}"
              f"{r['rel_err']*100:>9.2f}%{r['null_err']*100:>8.2f}%"
              f"{r['max_perp']:>10.2e}{r['max_drift']:>10.2e}  {note}")
    print(f"\n  tolerances: elastic traj <= {ELASTIC_TOL*100:.1f}% of omega0 "
          f"(& springless-null >= {NULL_MARGIN:.0f}x worse), "
          f"cap slope <= {SLOPE_TOL*100:.0f}% rel,")
    print(f"              omega_perp <= {PERP_TOL*100:.1f}% of omega0, "
          f"drift <= {DRIFT_TOL*1e6:.0f} um")
    print("RESULT:", "PASS" if ok else "FAIL")
    return ok


def plot(rows):
    os.makedirs(PLOT_DIR, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150, "font.size": 11})

    # ── elastic: omega(t) measured vs exact damped-oscillator ──
    elastic = [c for c in CASES if c["kind"] == "elastic"]
    if elastic:
        fig, axes = plt.subplots(1, len(elastic), figsize=(6.4 * len(elastic), 4.6),
                                 squeeze=False)
        for ax, c in zip(axes[0], elastic):
            path = roll_csv(c)
            if not os.path.isfile(path):
                continue
            ts, ws = [], []
            with open(path) as f:
                for row in csv.DictReader(f):
                    ts.append(float(row["t"]))
                    ws.append(float(row["omega_roll"]))
            ax.plot(ts, ws, "o", ms=2.5, mfc="none", color="tab:blue",
                    label="DIRT (sds)")
            te = [ts[0] + (ts[-1] - ts[0]) * k / 400 for k in range(401)]
            ax.plot(te, [elastic_omega(t, c["k_r"], c["gamma_r"], OMEGA0) for t in te],
                    "k-", lw=1.4, label="exact damped osc.")
            ax.plot(te, [elastic_springless(t, c["gamma_r"], OMEGA0) for t in te],
                    "r:", lw=1.2, label="springless (k=0)")
            lab, rate, zeta, wn = elastic_eigen(c["k_r"], c["gamma_r"])
            ax.axhline(0.0, color="0.6", lw=0.6)
            ax.set_title(f"{c['tag']}\n{lab}: zeta={zeta:.2f}, wn={wn:.0f}, rate={rate:.1f}/s",
                         fontsize=9)
            ax.set_xlabel("time t (s)")
            ax.set_ylabel(r"rolling spin $\omega$ (rad/s)")
            ax.legend(fontsize=8)
        fig.suptitle("SDS rolling — elastic regime: spring-dashpot decay vs exact")
        fig.tight_layout()
        fig.savefig(os.path.join(PLOT_DIR, "sds_rolling_elastic.png"))
        plt.close(fig)

    # ── cap: saturated spin-down vs alpha = (5/4) mu_r g / R ──
    caps = [c for c in CASES if c["kind"] == "cap"]
    if caps:
        fig, ax = plt.subplots(figsize=(6.6, 4.7))
        colors = plt.cm.viridis([0.15, 0.5, 0.82][: len(caps)])
        for c, col in zip(caps, colors):
            path = roll_csv(c)
            if not os.path.isfile(path):
                continue
            ts, ws = [], []
            with open(path) as f:
                for row in csv.DictReader(f):
                    ts.append(float(row["t"]))
                    ws.append(float(row["omega_roll"]))
            ax.plot(ts, ws, "-", color=col, lw=1.5,
                    label=fr"DIRT $\mu_r$={c['mu_r']}")
            ap = cap_alpha(c["mu_r"])
            tstop = OMEGA0 / ap
            ax.plot([0, tstop], [OMEGA0, 0.0], "k:", lw=1.0,
                    label="theory" if c is caps[0] else None)
        ax.set_xlabel("time t (s)")
        ax.set_ylabel(r"rolling spin $\omega$ (rad/s)")
        ax.set_title(r"SDS rolling — Coulomb cap: $\alpha=(5/4)\mu_r g/R$")
        ax.set_ylim(bottom=0)
        ax.legend(fontsize=9)
        fig.tight_layout()
        fig.savefig(os.path.join(PLOT_DIR, "sds_rolling_cap.png"))
        plt.close(fig)
    print(f"\nFigures -> {PLOT_DIR}/sds_rolling_elastic.png, sds_rolling_cap.png")


def graph():
    rows = _load(SWEEP_CSV)
    if not rows:
        print(f"No {SWEEP_CSV} — run 'start' first.")
        return False
    ok = validate(rows)
    plot(rows)
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
