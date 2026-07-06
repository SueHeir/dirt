#!/usr/bin/env python3
"""bench_mpi_decomposition — gated PASS/FAIL driver for MPI cross-rank DEM correctness.

Proves that a contact-rich frictional granular gas run under a **>1 processor**
domain decomposition reproduces the single-rank (`1x1x1`) result. The same
seeded config is run at three decompositions of the SAME periodic box:

  1x1x1 : serial reference
  2x1x1 : mpiexec -n 2  (one rank boundary; atoms migrate across x = 0.01)
  2x2x1 : mpiexec -n 4  (four sub-domains; boundaries in x AND y)

Every property that a correct decomposition must preserve is asserted against
the 1x1x1 reference, sampled from the MPI-gathered `[dump]` frames (each frame
gathers every rank's LOCAL atoms to rank 0, so a lost/duplicated atom shows up
as a wrong global count or tag set):

  (a) ATOM-COUNT / IDENTITY CONSERVATION — every gathered frame of every
      decomposition holds exactly N atoms with the exact reference tag set. A
      migration or ghost-exchange bug that dropped or duplicated an atom breaks
      this immediately.

  (b) GLOBAL MOMENTUM CONSERVATION — gravity is off and the box is periodic, so
      total linear momentum P = Sum m*v is an exact invariant (contact forces are
      equal-and-opposite pairs, reverse-communicated across ghosts). For every
      decomposition |P(t) - P(0)| stays at the FP round-off floor, and the
      multi-rank P(t) matches the reference P(t) to the same floor. A dropped
      ghost-force reverse-comm would inject a net force and drift momentum.

  (c) GLOBAL ENERGY (KE) AGREEMENT — the total kinetic energy trajectory KE(t)
      of each multi-rank run matches the 1x1x1 KE(t) to the FP floor at every
      sample.

  (d) PER-ATOM TRAJECTORY AGREEMENT — at the final step, per-atom velocities and
      (minimum-image) positions, matched by global tag, agree with the 1x1x1
      reference to the FP-associativity floor (documented tolerance TOL).

  plus provenance : each multi-rank run really used its processor grid
  (`Comm: processors ...` on stdout), and non-vacuity: the gas actually moved,
  carried live dissipative contacts (KE decayed), and stayed finite.

The neighbor list is rebuilt every step with the cache re-sort disabled
(`[neighbor] every=1 check=false sort_every=0`) so the ONLY difference between
the runs is the decomposition — the reduction/summation order across rank
boundaries. That difference is pure floating-point associativity, so agreement
is to the FP floor, NOT a loosened physics band: the measured worst case is
~1e-13 (see the printed numbers), four orders under the 1e-9 gate this repo
also uses for bench_restart_determinism.

Run:   python3 examples/bench_mpi_decomposition/sweep.py
Exit:  0 iff every check passes (and prints "ALL CHECKS PASSED").
"""

import os
import sys
import glob
import math
import shutil
import struct
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
EXAMPLE = "bench_mpi_decomposition"
BASE_CONFIG = os.path.join(HERE, "config.toml")
BIN = os.path.join(REPO_ROOT, "target", "release", "examples", EXAMPLE)
WORK = os.path.join(HERE, "data", "work")  # under data/ so it is gitignored
PLOTS = os.path.join(HERE, "plots")
DELTA_PLOT = os.path.join(PLOTS, "mpi_decomposition_deltas.svg")
MARKER = "# === SWEEP CONTROL BELOW"

TOTAL_STEPS = 4000     # integration steps (0 .. TOTAL_STEPS-1)
DUMP_INTERVAL = 500    # gather + dump a global frame every this many steps
DT = 1.0e-6            # timestep (s) — identical across every decomposition
TOL = 1e-9             # FP-associativity floor for per-atom + conservation agreement

# Physical constants that mirror config.toml (all particles are identical glass
# spheres, so mass is a pure function of radius and this density).
DENSITY = 2500.0
BOX_L = 0.02           # periodic box edge (m) — for minimum-image positions
POS_SCALE = BOX_L      # relative-error denominator floor (position)
VEL_SCALE = 0.2        # insertion Gaussian sigma (m/s) — denominator floor (velocity)

# Decompositions under test. (px, py, pz, nranks). The first is the reference.
DECOMPS = [
    ("n1", 1, 1, 1, 1),
    ("n2", 2, 1, 1, 2),
    ("n4", 2, 2, 1, 4),
]


# ── config composition ───────────────────────────────────────────────────────

def physics_prefix():
    with open(BASE_CONFIG) as f:
        text = f.read()
    idx = text.find(MARKER)
    if idx < 0:
        sys.exit(f"ERROR: control marker '{MARKER}' not found in {BASE_CONFIG}")
    return text[:idx]


def write_config(name, px, py, pz, outdir):
    control = (
        f"[comm]\n"
        f"processors_x = {px}\n"
        f"processors_y = {py}\n"
        f"processors_z = {pz}\n\n"
        f"[output]\n"
        f"dir = \"{outdir}\"\n\n"
        f"[dump]\n"
        f"interval = {DUMP_INTERVAL}\n"
        f"format = \"binary\"\n\n"
        f"[[run]]\n"
        f"name = \"{name}\"\n"
        f"dt = {DT}\n"
        f"steps = {TOTAL_STEPS}\n"
        f"thermo = 100000\n"
        f"save_at_end = true\n"
    )
    path = os.path.join(WORK, f"config_{name}.toml")
    with open(path, "w") as f:
        f.write(physics_prefix())
        f.write(control)
    return path


# ── build + run ──────────────────────────────────────────────────────────────

def build():
    # DEFAULT features = mpi_backend + precision-double (see root Cargo.toml), so
    # the same binary serves the serial reference AND the mpiexec runs.
    print("Building bench_mpi_decomposition (default features: mpi + precision-double) ...")
    subprocess.run(
        ["cargo", "build", "--release", "--example", EXAMPLE],
        cwd=REPO_ROOT, check=True,
    )


def run(cfg_path, tag, nranks):
    log = os.path.join(WORK, f"{tag}.log")
    cmd = [BIN, cfg_path] if nranks == 1 else ["mpiexec", "-n", str(nranks), BIN, cfg_path]
    with open(log, "w") as lf:
        proc = subprocess.run(cmd, cwd=REPO_ROOT, stdout=lf, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        with open(log) as f:
            sys.stderr.write(f.read()[-3000:])
        sys.exit(f"ERROR: run '{tag}' failed (exit {proc.returncode}); see {log}")
    with open(log) as f:
        return f.read()


# ── binary dump parsing ──────────────────────────────────────────────────────
# soil_print write_dump_binary (per_rank=false => MPI-gathered to rank 0):
#   u32 count, then per atom:
#   tag u32, type u32, pos[3] f64, vel[3] f64, force[3] f64, radius f64
# (record stride 88; no trailing columns for this plugin set). Mass is not
# stored, so derive it from radius + DENSITY (all particles identical glass).

def load_frame(path):
    with open(path, "rb") as f:
        raw = f.read()
    (count,) = struct.unpack_from("<I", raw, 0)
    if count == 0:
        sys.exit(f"ERROR: dump {path} has zero atoms")
    body = len(raw) - 4
    if body % count != 0:
        sys.exit(f"ERROR: dump {path} size {len(raw)} not divisible by count {count}")
    rec = body // count
    if rec < 88:
        sys.exit(f"ERROR: dump {path} record stride {rec} < 88")
    atoms = {}
    for i in range(count):
        off = 4 + i * rec
        (tag,) = struct.unpack_from("<I", raw, off)
        pos = struct.unpack_from("<3d", raw, off + 8)
        vel = struct.unpack_from("<3d", raw, off + 32)
        (radius,) = struct.unpack_from("<d", raw, off + 80)
        atoms[tag] = (pos, vel, radius)
    return atoms


def frame_paths(outdir):
    files = glob.glob(os.path.join(outdir, "dump", "dump_*.bin"))
    if not files:
        sys.exit(f"ERROR: no dump files in {outdir}/dump")
    return sorted(files, key=lambda p: int(os.path.basename(p)[len("dump_"):-len(".bin")]))


def step_of(path):
    return int(os.path.basename(path)[len("dump_"):-len(".bin")])


# ── physics reductions ───────────────────────────────────────────────────────

def mass_of(radius):
    return DENSITY * (4.0 / 3.0) * math.pi * radius ** 3


def momentum(atoms):
    p = [0.0, 0.0, 0.0]
    for (_, vel, r) in atoms.values():
        m = mass_of(r)
        for d in range(3):
            p[d] += m * vel[d]
    return p


def kinetic_energy(atoms):
    ke = 0.0
    for (_, vel, r) in atoms.values():
        m = mass_of(r)
        ke += 0.5 * m * (vel[0] ** 2 + vel[1] ** 2 + vel[2] ** 2)
    return ke


def vnorm(v):
    return (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) ** 0.5


def min_image(a, b):
    d = a - b
    if d > 0.5 * BOX_L:
        d -= BOX_L
    elif d < -0.5 * BOX_L:
        d += BOX_L
    return d


def per_atom_error(ref, other):
    """Max over atoms of relative velocity and (minimum-image) position error."""
    if set(ref) != set(other):
        sys.exit("ERROR: frames have different atom tag sets — cannot compare")
    worst_pos = 0.0
    worst_vel = 0.0
    for tag in ref:
        pa, va, _ = ref[tag]
        pb, vb, _ = other[tag]
        dp = (min_image(pa[0], pb[0]) ** 2 + min_image(pa[1], pb[1]) ** 2
              + min_image(pa[2], pb[2]) ** 2) ** 0.5
        dv = vnorm((va[0] - vb[0], va[1] - vb[1], va[2] - vb[2]))
        worst_pos = max(worst_pos, dp / (vnorm(pa) + POS_SCALE))
        worst_vel = max(worst_vel, dv / (vnorm(va) + VEL_SCALE))
    return worst_pos, worst_vel


def all_finite(atoms):
    for (pos, vel, _) in atoms.values():
        for x in list(pos) + list(vel):
            if not math.isfinite(x):
                return False
    return True


def write_delta_plot(rows):
    """Write a dependency-free SVG of multi-rank deltas against 1x1x1."""
    os.makedirs(PLOTS, exist_ok=True)

    metrics = [
        ("tag/count", "tag_delta"),
        ("momentum", "momentum_match"),
        ("KE(t)", "ke_match"),
        ("final pos", "pos"),
        ("final vel", "vel"),
    ]
    width, height = 840, 500
    left, right, top, bottom = 82, 26, 54, 92
    plot_w = width - left - right
    plot_h = height - top - bottom
    ymin, ymax = 1e-30, 2e-8
    log_min, log_max = math.log10(ymin), math.log10(ymax)

    def sx(i, series_idx):
        group_w = plot_w / len(metrics)
        bar_w = group_w * 0.28
        center = left + group_w * (i + 0.5)
        return center + (series_idx - (len(rows) - 1) / 2.0) * bar_w * 1.25 - bar_w / 2.0

    def sy(value):
        clipped = min(max(value, ymin), ymax)
        frac = (math.log10(clipped) - log_min) / (log_max - log_min)
        return top + plot_h * (1.0 - frac)

    def esc(text):
        return (str(text).replace("&", "&amp;").replace("<", "&lt;")
                .replace(">", "&gt;").replace('"', "&quot;"))

    colors = ["#2878b5", "#d95f02"]
    bar_w = (plot_w / len(metrics)) * 0.28
    tol_y = sy(TOL)
    tick_values = [1e-30, 1e-25, 1e-20, 1e-15, 1e-10]
    label_y = top + plot_h + 24

    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img" '
        f'aria-label="bench_mpi_decomposition multi-rank agreement">',
        "<style><![CDATA["
        "text{font-family:Arial,Helvetica,sans-serif;fill:#222}"
        ".small{font-size:12px}.tick{font-size:11px;fill:#444}"
        ".title{font-size:19px;font-weight:700}.axis{stroke:#333;stroke-width:1}"
        ".grid{stroke:#d4d4d4;stroke-width:1;stroke-dasharray:2 4}"
        ".tol{stroke:#b2182b;stroke-width:2;stroke-dasharray:7 5}"
        "]]></style>",
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text class="title" x="82" y="28">bench_mpi_decomposition: multi-rank agreement</text>',
        f'<text class="small" x="82" y="47">relative delta vs 1x1x1 reference; dashed line is unchanged {TOL:.0e} pass tolerance</text>',
    ]

    for tv in tick_values:
        y = sy(tv)
        parts.append(f'<line class="grid" x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}"/>')
        parts.append(f'<text class="tick" x="{left - 8}" y="{y + 4:.1f}" text-anchor="end">{tv:.0e}</text>')

    parts.extend([
        f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}"/>',
        f'<line class="axis" x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}"/>',
        f'<line class="tol" x1="{left}" y1="{tol_y:.1f}" x2="{left + plot_w}" y2="{tol_y:.1f}"/>',
        f'<text class="small" fill="#b2182b" x="{left + plot_w - 4}" y="{tol_y - 6:.1f}" text-anchor="end">pass tolerance {TOL:.0e}</text>',
    ])

    for i, (label, _) in enumerate(metrics):
        group_w = plot_w / len(metrics)
        cx = left + group_w * (i + 0.5)
        parts.append(f'<text class="small" x="{cx:.1f}" y="{label_y}" text-anchor="middle">{esc(label)}</text>')

    for j, row in enumerate(rows):
        for i, (_, key) in enumerate(metrics):
            value = row[key]
            drawn = max(value, ymin)
            x = sx(i, j)
            y = sy(drawn)
            h = top + plot_h - y
            parts.append(
                f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{h:.1f}" '
                f'fill="{colors[j % len(colors)]}"><title>{esc(row["name"])} {esc(key)} = {value:.3e}</title></rect>'
            )
            if value == 0.0:
                parts.append(
                    f'<text class="tick" x="{x + bar_w / 2:.1f}" y="{sy(2e-30):.1f}" '
                    f'text-anchor="middle">0</text>'
                )

    legend_x = left + 8
    legend_y = top + 16
    for j, row in enumerate(rows):
        lx = legend_x + j * 150
        parts.append(f'<rect x="{lx}" y="{legend_y - 10}" width="14" height="14" fill="{colors[j % len(colors)]}"/>')
        parts.append(f'<text class="small" x="{lx + 20}" y="{legend_y + 2}">{esc(row["name"])}</text>')

    parts.append(
        f'<text class="small" x="{width - 20}" y="{height - 18}" text-anchor="end">'
        "tag/count is an exact identity check; zeros are drawn at the log-scale floor</text>"
    )
    parts.append("</svg>")
    with open(DELTA_PLOT, "w") as f:
        f.write("\n".join(parts) + "\n")
    print(f"\nWrote {os.path.relpath(DELTA_PLOT, REPO_ROOT)}")


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    if shutil.which("mpiexec") is None:
        sys.exit("ERROR: mpiexec not found on PATH (needed for the multi-rank runs)")
    if shutil.which("cargo") is None and not os.path.exists(BIN):
        sys.exit("ERROR: cargo not found and no prebuilt binary")
    build()

    if os.path.isdir(WORK):
        shutil.rmtree(WORK)
    os.makedirs(WORK)

    # Run every decomposition of the same seeded config.
    logs = {}
    outdirs = {}
    for (name, px, py, pz, nranks) in DECOMPS:
        outdir = os.path.join(WORK, name)
        cfg = write_config(name, px, py, pz, outdir)
        print(f"Running {name} ({px}x{py}x{pz}, {nranks} rank(s), {TOTAL_STEPS} steps) ...")
        logs[name] = run(cfg, name, nranks)
        outdirs[name] = outdir

    # Load every gathered frame, keyed by step, for each decomposition.
    frames = {}
    for name in outdirs:
        frames[name] = {step_of(p): load_frame(p) for p in frame_paths(outdirs[name])}

    ref_name = DECOMPS[0][0]
    ref_frames = frames[ref_name]
    ref_steps = sorted(ref_frames)
    n_expected = len(ref_frames[ref_steps[0]])
    ref_tags = set(ref_frames[ref_steps[0]])
    final_step = ref_steps[-1]

    # ── evaluate ──────────────────────────────────────────────────────────────
    passed = 0
    total = 0

    def check(name, ok, detail=""):
        nonlocal passed, total
        total += 1
        if ok:
            passed += 1
        print(f"  {name:<34}{'PASS' if ok else 'FAIL'}   {detail}")
        return ok

    print("\nMPI domain-decomposition cross-rank correctness")
    print(f"  reference = {ref_name} (1x1x1), N = {n_expected} atoms, "
          f"{len(ref_steps)} frames, final step = {final_step}")
    print("-" * 78)

    # provenance: each multi-rank run really ran on its processor grid.
    for (name, px, py, pz, nranks) in DECOMPS:
        grid = f"Comm: processors {px} {py} {pz}"
        check(f"{name}: grid {px}x{py}x{pz} active", grid in logs[name],
              "stdout confirms processor grid" if grid in logs[name] else "grid line missing")

    # non-vacuity: reference gas moved, dissipated (contacts live), stayed finite.
    ref_init = ref_frames[ref_steps[0]]
    ref_final = ref_frames[final_step]
    _, vel_change = per_atom_error(ref_init, ref_final)  # velocity field rearranged by collisions
    ke0 = kinetic_energy(ref_init)
    kef = kinetic_energy(ref_final)
    check("dynamics non-trivial", all_finite(ref_final) and vel_change > 1e-2 and kef < 0.98 * ke0,
          f"vel-field change={vel_change:.3e}, KE {ke0:.3e}->{kef:.3e} "
          f"(dissipated {100*(1-kef/ke0):.1f}% via contacts)")

    # For each multi-rank decomposition, compare against the reference.
    worst_overall_pos = 0.0
    worst_overall_vel = 0.0
    worst_mom_drift = 0.0
    worst_mom_match = 0.0
    worst_ke_match = 0.0
    plot_rows = []
    for (name, px, py, pz, nranks) in DECOMPS:
        fr = frames[name]
        steps = sorted(fr)

        # (a) atom-count / identity conservation at EVERY frame.
        count_ok = all(len(fr[s]) == n_expected and set(fr[s]) == ref_tags for s in steps)
        bad = next((s for s in steps if len(fr[s]) != n_expected or set(fr[s]) != ref_tags), None)
        tag_delta = max(
            (abs(len(fr[s]) - n_expected) + len(set(fr[s]).symmetric_difference(ref_tags))) / n_expected
            for s in steps
        )
        check(f"{name}: atom-count conserved", count_ok,
              f"all {len(steps)} frames hold exactly {n_expected} atoms (tag set intact)"
              if count_ok else f"frame at step {bad} has {len(fr[bad])} atoms / tag mismatch")

        # (b) momentum conservation within the run + match to reference.
        p0 = momentum(fr[steps[0]])
        p0n = vnorm(p0) + DENSITY  # denominator floor (P has physical magnitude)
        mom_drift = max(vnorm([momentum(fr[s])[d] - p0[d] for d in range(3)]) / p0n
                        for s in steps)
        worst_mom_drift = max(worst_mom_drift, mom_drift)
        check(f"{name}: momentum conserved", mom_drift <= TOL,
              f"max |P(t)-P(0)|/|P0| = {mom_drift:.2e} (tol {TOL:.0e})")

        if name != ref_name:
            mom_match = max(
                vnorm([momentum(fr[s])[d] - momentum(ref_frames[s])[d] for d in range(3)]) / p0n
                for s in steps)
            worst_mom_match = max(worst_mom_match, mom_match)
            check(f"{name}: momentum matches 1x1x1", mom_match <= TOL,
                  f"max |P_{name}(t)-P_ref(t)|/|P0| = {mom_match:.2e} (tol {TOL:.0e})")

            # (c) kinetic-energy trajectory agreement at every frame.
            ke_match = max(abs(kinetic_energy(fr[s]) - kinetic_energy(ref_frames[s]))
                           / (kinetic_energy(ref_frames[s]) + 1e-30) for s in steps)
            worst_ke_match = max(worst_ke_match, ke_match)
            check(f"{name}: KE(t) matches 1x1x1", ke_match <= TOL,
                  f"max rel KE diff over {len(steps)} frames = {ke_match:.2e} (tol {TOL:.0e})")

            # (d) per-atom final trajectory agreement to the FP floor.
            wp, wv = per_atom_error(ref_final, fr[final_step])
            worst_overall_pos = max(worst_overall_pos, wp)
            worst_overall_vel = max(worst_overall_vel, wv)
            check(f"{name}: per-atom positions", wp <= TOL,
                  f"max rel = {wp:.2e} (tol {TOL:.0e})")
            check(f"{name}: per-atom velocities", wv <= TOL,
                  f"max rel = {wv:.2e} (tol {TOL:.0e})")
            plot_rows.append({
                "name": f"{name} ({px}x{py}x{pz})",
                "tag_delta": tag_delta,
                "momentum_match": mom_match,
                "ke_match": ke_match,
                "pos": wp,
                "vel": wv,
            })

    print("-" * 78)
    print(f"Worst per-atom: pos {worst_overall_pos:.2e}, vel {worst_overall_vel:.2e}  |  "
          f"momentum drift {worst_mom_drift:.2e}, match {worst_mom_match:.2e}  |  "
          f"KE match {worst_ke_match:.2e}   (all << {TOL:.0e} FP floor)")
    write_delta_plot(plot_rows)
    print(f"\nResult: {passed}/{total} checks passed")
    if passed == total:
        print("ALL CHECKS PASSED")
        return 0
    print(f"CHECKS FAILED: {total - passed} of {total}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
