#!/usr/bin/env python3
"""
Marshall twisting-model spin-down benchmark driver.

This validates DIRT's `twisting_model = "marshall"` — the twisting model whose
stiffness/damping/friction coefficients are DERIVED from the active tangential
(Mindlin) contact model, exactly as LAMMPS `pair_granular` `twisting marshall`
does (doc/src/pair_granular.rst, "twisting_marshall"; Marshall 2009 eqs 32-33):

    k_twist  = 0.5 * k_t * a^2
    gamma_tw = 0.5 * gamma_t * a^2
    mu_twist = (2/3) * a * mu_t

where a = sqrt(R* * delta) is the Hertz contact radius, k_t / gamma_t are the
tangential (Mindlin) spring/damping, and mu_t is the tangential (SLIDING)
friction coefficient. This is distinct from DIRT's constant / SDS twisting
models (validated in bench_twisting_friction), whose coefficients are direct
user inputs.

## Test geometry (same seating as bench_twisting_friction)
Two identical spheres are stacked along +z. The lower sphere is a frozen anchor
([[freeze]] — it can neither translate nor spin), and the upper sphere rests on
it under gravity, seated at the static Hertz overlap so the normal contact force
is F_n = m*g. The upper sphere is given a PURE twisting spin omega = (0,0,omega0)
about the contact normal n = +z.

A pure twist about n produces zero relative surface velocity at the contact
point, so there is no sliding and no rolling — the tangential spring stays at
zero even though mu_t > 0. The ONLY contact torque is the Marshall twisting
couple, which in the saturated (sliding) regime is

    tau_tw = mu_twist * F_n = (2/3) * a * mu_t * F_n     (opposing the twist).

With F_n = m*g, equal-sphere r_eff = R/2 (so a = sqrt(R*delta/2)), and sphere
inertia I = (2/5) m R^2, the spin decays at the EXACT constant rate

    alpha = domega/dt = (5/3) * a * mu_t * g / R^2.

The tangential (sliding) friction mu_t is the SWEPT variable: because
mu_twist = (2/3) a mu_t, doubling mu_t doubles the twisting cap and hence the
spin-down rate. This is the essential Marshall behaviour — the twist couple is
tied to the sliding coefficient, not to independent twist inputs.

This script sweeps a few mu_t, fits the omega_z(t) slope, and validates it
against alpha; PASS requires the fitted deceleration to match within tolerance
AND the twist to have stayed pure (omega_perp ~= 0, no lateral drift).

Commands (from anywhere):
    python3 examples/bench_marshall_twisting/sweep.py generate   # write per-case configs
    python3 examples/bench_marshall_twisting/sweep.py start      # build + run all sims -> CSV
    python3 examples/bench_marshall_twisting/sweep.py graph      # validate + plot
    python3 examples/bench_marshall_twisting/sweep.py            # all three, in order

Outputs:
    sweep/<case>/config.toml       DIRT configs                       (gitignored)
    data/decay_<case>.csv          per-case DIRT time series          (gitignored)
    data/sweep.csv                 fitted-slope summary               (gitignored)
    plots/*.png                    final figures                      (tracked)

References:
  * LAMMPS pair_granular documentation, "twisting_marshall" derived coefficients
    (~/projects/reference/lammps/doc/src/pair_granular.rst).
  * W.R. Marshall, "Discrete-element modeling of particulate aerosol flows",
    J. Comput. Phys. 228 (2009) 1541-1561, eqs 32-33 (the k_twist, eta_twist,
    mu_twist expressions in terms of the tangential coefficients).
"""

import math
import os
import sys
import csv
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_marshall_twisting"

SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
SWEEP_CSV = os.path.join(DATA_DIR, "sweep.csv")

# ── Material / geometry (fixed across the sweep) ──────────────────────────────
NU = 0.3              # Poisson ratio
E_N = 0.3             # normal restitution (damps the tiny settling transient)
YOUNGS_MOD = 1.0e8    # Pa — soft (keeps dt large; a and the rate follow from the seating)
DENSITY = 2500.0      # kg/m^3
RADIUS = 0.005        # m — sphere radius R
GRAVITY = 9.81        # m/s^2
OMEGA0 = 8.0          # rad/s — initial twisting spin (must match main.rs OMEGA0)
DT = 1.0e-5           # s
STEPS = 40000         # plenty for the slowest (smallest mu_t) case to spin down

# Static Hertz overlap that makes F_n = m g (must match main.rs SEAT_OVERLAP).
# The Marshall contact radius a = sqrt(r_eff * delta) is evaluated at this
# overlap, so it appears directly in the analytical spin-down rate below.
SEAT_OVERLAP = 2.307e-6  # m

# Swept TANGENTIAL (sliding) friction coefficients. Marshall sets
# mu_twist = (2/3) a mu_t, so this is what tunes the twisting cap. Kept modest so
# every case spins down within the step budget.
MU_T_LIST = [0.5, 0.7, 0.9]

R_EFF = RADIUS / 2.0                        # equal spheres
A_CONTACT = math.sqrt(R_EFF * SEAT_OVERLAP)  # Hertz contact radius a = sqrt(R* delta)


# ── Theory ────────────────────────────────────────────────────────────────────
def a_pred(mu_t):
    """Exact constant Marshall twisting spin-down rate. The saturated couple is
    tau = mu_twist F_n = (2/3) a mu_t (m g); with I = (2/5) m R^2:
        alpha = tau / I = (5/3) * a * mu_t * g / R^2."""
    return (5.0 / 3.0) * A_CONTACT * mu_t * GRAVITY / (RADIUS * RADIUS)


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
twisting_model = "marshall"
[[dem.materials]]
name = "grain"
youngs_mod = {youngs:.6e}
poisson_ratio = {nu}
restitution = {e_n}
friction = {mu_t}
twisting_friction = 0.0
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
def case_tag(mu_t):
    return f"marshall_mu_{mu_t:g}"


def case_dir(mu_t):
    return os.path.join(SWEEP_DIR, case_tag(mu_t))


def decay_csv(mu_t):
    return os.path.join(DATA_DIR, f"decay_{case_tag(mu_t)}.csv")


def _dirt_config(mu_t, outdir):
    return TOML_TEMPLATE.format(
        g=GRAVITY, youngs=YOUNGS_MOD, nu=NU, e_n=E_N, mu_t=mu_t,
        radius=RADIUS, density=DENSITY, outdir=outdir, steps=STEPS, dt=DT,
    )


# ── generate ──────────────────────────────────────────────────────────────────
def generate():
    n = 0
    for mu_t in MU_T_LIST:
        cdir = case_dir(mu_t)
        os.makedirs(cdir, exist_ok=True)
        with open(os.path.join(cdir, "config.toml"), "w") as f:
            f.write(_dirt_config(mu_t, cdir))
        n += 1
    print(f"Generated {n} DIRT sweep configs under {SWEEP_DIR}")


# ── start ─────────────────────────────────────────────────────────────────────
SWEEP_FIELDS = ["mu_t", "a_fit", "a_pred", "rel_err", "max_perp", "max_drift", "npts"]


def _read_timeseries(path):
    """Read t,omega_z,omega_perp,drift from a twisting_results.csv (skip '#')."""
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
    return rows  # list of (t, omega_z, omega_perp, drift)


def _fit_decay(rows):
    """Linear-fit omega_z(t) over the monotone decay window; return
    (a_fit, max_perp, max_drift, npts).

    a_fit = -slope (deceleration, positive). The window is
    0.15*omega0 <= omega_z <= 0.85*omega0 — above the near-zero discretization
    chatter and below the brief Marshall spring wind-up. max_perp is the largest
    off-axis spin |omega_perp|/omega0 and max_drift the largest lateral
    displacement (m) over that window (twist-purity checks)."""
    lo, hi = 0.15 * OMEGA0, 0.85 * OMEGA0
    prefix = []
    for r in rows:
        if r[1] < lo:
            break
        prefix.append(r)
    win = [r for r in prefix if lo <= r[1] <= hi]
    if len(win) < 10:
        win = [r for r in prefix if r[1] > 0.01 * OMEGA0] or prefix
    n = len(win)
    if n < 2:
        return None
    ts = [r[0] for r in win]
    ws = [r[1] for r in win]
    tbar = sum(ts) / n
    wbar = sum(ws) / n
    sxx = sum((t - tbar) ** 2 for t in ts)
    sxy = sum((t - tbar) * (w - wbar) for t, w in zip(ts, ws))
    slope = sxy / sxx if sxx > 0 else 0.0
    a_fit = -slope
    max_perp = max(r[2] for r in win) / OMEGA0
    max_drift = max(r[3] for r in win)
    return a_fit, max_perp, max_drift, n


def _run_dirt(mu_t):
    cdir = case_dir(mu_t)
    config = os.path.join(cdir, "config.toml")
    res = os.path.join(cdir, "data", "twisting_results.csv")
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
    rows = _read_timeseries(res)
    with open(decay_csv(mu_t), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t", "omega_z", "omega_perp", "drift"])
        w.writerows(rows)
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
    subprocess.run(["cargo", "build", "--release", "--example", EXAMPLE,
                    "--no-default-features", "--features", "precision-double"],
                   cwd=REPO_ROOT, check=True)

    if os.path.exists(SWEEP_CSV):
        os.remove(SWEEP_CSV)

    out_rows = []
    n = len(MU_T_LIST)
    for i, mu_t in enumerate(MU_T_LIST, 1):
        cdir = case_dir(mu_t)
        if not os.path.isfile(os.path.join(cdir, "config.toml")):
            print(f"  [{i}/{n}] missing config for {case_tag(mu_t)} — run 'generate' first.")
            continue
        print(f"  [{i}/{n}] marshall mu_t={mu_t:<5}", end="  ", flush=True)
        rows = _run_dirt(mu_t)
        if rows is None:
            print("DIRT FAILED")
            continue
        fit = _fit_decay(rows)
        if fit is None:
            print("DIRT (no fit)")
            continue
        a_fit, max_perp, max_drift, npts = fit
        ap = a_pred(mu_t)
        rel = abs(a_fit - ap) / ap
        out_rows.append({"mu_t": mu_t, "a_fit": a_fit, "a_pred": ap,
                         "rel_err": rel, "max_perp": max_perp,
                         "max_drift": max_drift, "npts": npts})
        print(f"a_fit={a_fit:.3f} a_pred={ap:.3f} rel={rel*100:.2f}% "
              f"perp={max_perp:.1e} drift={max_drift:.1e}")

    if not out_rows:
        print("\nERROR: no DIRT results collected.")
        sys.exit(1)
    _write_csv(SWEEP_CSV, SWEEP_FIELDS, out_rows)
    print(f"\nDIRT: {len(out_rows)}/{n} cases -> {SWEEP_CSV}")


# ── graph (validate + plot) ───────────────────────────────────────────────────
SLOPE_TOL = 0.04      # 4% relative error on the fitted deceleration (theory is exact)
PERP_TOL = 1.0e-3     # |omega_perp| must stay < 0.1% of omega0 (pure twist)
DRIFT_TOL = 1.0e-5    # lateral drift must stay < 10 um (co-axial, no translation)


def _load(path):
    if not os.path.isfile(path):
        return []
    out = []
    with open(path) as f:
        for r in csv.DictReader(f):
            row = dict(r)
            for k in ("mu_t", "a_fit", "a_pred", "rel_err", "max_perp",
                      "max_drift", "npts"):
                row[k] = float(row[k])
            out.append(row)
    return out


def validate(rows):
    print("\n=== Marshall twisting spin-down validation ===")
    print(f"  R={RADIUS} m  equal spheres (r_eff = R/2)  g={GRAVITY}  omega0={OMEGA0}")
    print(f"  Hertz contact radius a = sqrt(r_eff*delta) = {A_CONTACT:.6e} m "
          f"(delta = {SEAT_OVERLAP:.3e} m)")
    print(f"  Marshall: mu_twist = (2/3) a mu_t  =>  alpha = (5/3) a mu_t g / R^2  (exact)")
    print(f"  {'mu_t':>7}{'a_fit':>10}{'a_pred':>10}"
          f"{'rel_err':>9}{'max_perp':>10}{'drift':>10}  note")
    ok = True
    npass = 0
    for r in sorted(rows, key=lambda x: x["mu_t"]):
        note = ""
        good = True
        if r["rel_err"] > SLOPE_TOL:
            note = "SLOPE MISMATCH"; good = False
        if r["max_perp"] > PERP_TOL:
            note = (note + " OFF-AXIS").strip(); good = False
        if r["max_drift"] > DRIFT_TOL:
            note = (note + " DRIFTED").strip(); good = False
        ok = ok and good
        npass += 1 if good else 0
        print(f"  {r['mu_t']:>7.3f}{r['a_fit']:>10.3f}"
              f"{r['a_pred']:>10.3f}{r['rel_err']*100:>8.2f}%"
              f"{r['max_perp']:>10.2e}{r['max_drift']:>10.2e}  {note}")
    print(f"\n  tolerances: slope <= {SLOPE_TOL*100:.0f}% rel, "
          f"omega_perp <= {PERP_TOL*100:.1f}% of omega0, drift <= {DRIFT_TOL*1e6:.0f} um")
    print(f"  {npass}/{len(rows)} checks passed")
    print("RESULT:", "PASS" if ok else "FAIL")
    print("ALL CHECKS PASSED" if ok else "CHECKS FAILED")
    return ok


def plot(rows):
    os.makedirs(PLOT_DIR, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150, "font.size": 11})

    # ── omega_z(t) spin-down curves: measured vs theory, per case ──
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    colors = plt.cm.viridis([0.15, 0.5, 0.82])
    cmap = {mu: c for mu, c in zip(MU_T_LIST, colors)}
    for mu_t in MU_T_LIST:
        path = decay_csv(mu_t)
        if not os.path.isfile(path):
            continue
        ts, ws = [], []
        with open(path) as f:
            for row in csv.DictReader(f):
                ts.append(float(row["t"]))
                ws.append(float(row["omega_z"]))
        ax.plot(ts, ws, "-", color=cmap[mu_t], lw=1.5,
                label=fr"marshall $\mu_t$={mu_t}")
    # Theory lines omega(t) = omega0 - a_pred*t, clamped at 0.
    for mu_t in MU_T_LIST:
        ap = a_pred(mu_t)
        tstop = OMEGA0 / ap
        ax.plot([0, tstop], [OMEGA0, 0.0], "k:", lw=1.0,
                label="theory" if mu_t == MU_T_LIST[0] else None)
    ax.set_xlabel("time t (s)")
    ax.set_ylabel(r"twist rate $\omega_z$ (rad/s)")
    ax.set_title(r"Marshall twisting spin-down: measured vs $\alpha=(5/3)a\mu_t g/R^2$")
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "marshall_spindown.png"))
    plt.close(fig)

    # ── fitted deceleration vs theory ──
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    mu_line = sorted(set(r["mu_t"] for r in rows))
    ax.plot(mu_line, [a_pred(m) for m in mu_line], "k-",
            label=r"theory $(5/3)a g/R^2$")
    d = sorted(rows, key=lambda x: x["mu_t"])
    ax.plot([r["mu_t"] for r in d], [r["a_fit"] for r in d], "o", ms=8,
            color="tab:blue", label="DIRT marshall (fit)")
    ax.set_xlabel(r"tangential friction $\mu_t$")
    ax.set_ylabel(r"spin-down rate $\alpha$ (rad/s$^2$)")
    ax.set_title("Fitted Marshall twisting spin-down vs theory")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "marshall_spindown_vs_mu_t.png"))
    plt.close(fig)
    print(f"\nFigures -> {PLOT_DIR}/marshall_spindown.png, marshall_spindown_vs_mu_t.png")


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
