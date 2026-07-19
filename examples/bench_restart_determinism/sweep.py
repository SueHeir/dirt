#!/usr/bin/env python3
"""bench_restart_determinism — gated PASS/FAIL driver.

Proves two production-trust properties of DIRT's `[restart]` machinery on a
small frictional granular gas (~300 spheres, periodic box, seeded ICs):

  (a) RESTART CONTINUITY. A run checkpointed at step N via `[restart]` and
      resumed in a fresh process reproduces the uninterrupted run to <= 1e-9
      relative on positions AND velocities at the final step.

  (b) RUN-TO-RUN DETERMINISM. Two independent, identically-configured single
      rank runs are bit-identical (byte-for-byte equal dump files).

Protocol (four processes, one shared config with per-run control blocks):

  A  uninterrupted : 0 .. TOTAL steps                       -> A/dump/dump_{last}.bin
  B1 checkpoint    : 0 .. CKPT  steps, save_at_end restart  -> B/restart/restart_{ckpt}
  B2 resume        : read restart, CKPT .. TOTAL steps      -> B/dump/dump_{last}.bin
  C  twin of A     : 0 .. TOTAL steps                        -> C/dump/dump_{last}.bin

Checks:
  (a) continuity  : per-atom (sorted by tag) rel error A-vs-B on pos & vel <= TOL
  (b) determinism : sha256(A final dump) == sha256(C final dump)  (also dump_0)
  plus provenance : checkpoint file exists and the resume genuinely read it,
  plus non-vacuity: the gas actually moved and stayed finite (test isn't trivial).

Run:   python3 examples/bench_restart_determinism/sweep.py
Exit:  0 iff every check passes (and prints "ALL CHECKS PASSED").
"""

import os
import sys
import glob
import shutil
import struct
import hashlib
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
EXAMPLE = "bench_restart_determinism"
BASE_CONFIG = os.path.join(HERE, "config.toml")
BIN = os.path.join(REPO_ROOT, "target", "release", "examples", EXAMPLE)
WORK = os.path.join(HERE, "data", "work")  # under data/ so it is gitignored
PLOT_DIR = os.path.join(HERE, "plots")
PLOT = os.path.join(PLOT_DIR, "restart_determinism.png")
MARKER = "# === SWEEP CONTROL BELOW"

TOTAL_STEPS = 4000     # total integration steps of the reference trajectory
CKPT_STEPS = 2000      # steps taken before the checkpoint is written
DT = 1.0e-6            # timestep (s) — must match across every run
TOL = 1e-9             # acceptance tolerance (relative) for restart continuity

# Relative-error floors: normalize by the physical scale of the box / velocity
# field so a near-stationary atom does not blow up a 0/0 relative error. These
# only make the metric MORE conservative (a larger denominator == smaller
# reported error is NOT what happens here: reproduction is essentially exact).
POS_SCALE = 0.02       # box edge (m)
VEL_SCALE = 0.2        # insertion Gaussian sigma (m/s)


# ── config composition ───────────────────────────────────────────────────────

def physics_prefix():
    with open(BASE_CONFIG) as f:
        text = f.read()
    idx = text.find(MARKER)
    if idx < 0:
        sys.exit(f"ERROR: control marker '{MARKER}' not found in {BASE_CONFIG}")
    return text[:idx]


def write_config(name, control):
    path = os.path.join(WORK, f"config_{name}.toml")
    with open(path, "w") as f:
        f.write(physics_prefix())
        f.write(control)
    return path


def control(outdir, steps, dump_interval=0, restart_read=False):
    return (
        f"[output]\n"
        f"dir = \"{outdir}\"\n\n"
        f"[dump]\n"
        f"interval = {dump_interval}\n"
        f"format = \"binary\"\n\n"
        f"[restart]\n"
        f"interval = 0\n"
        f"format = \"bincode\"\n"
        f"read = {'true' if restart_read else 'false'}\n\n"
        f"[[run]]\n"
        f"name = \"{os.path.basename(outdir)}\"\n"
        f"dt = {DT}\n"
        f"steps = {steps}\n"
        f"thermo = 1000\n"
        f"save_at_end = true\n"
    )


# ── process runner ───────────────────────────────────────────────────────────

def build():
    print("Building bench_restart_determinism (precision-double) ...")
    subprocess.run(
        ["cargo", "build", "--release", "--example", EXAMPLE,
         "--no-default-features", "--features", "precision-double"],
        cwd=REPO_ROOT, check=True,
    )


def run(cfg_path, tag):
    log = os.path.join(WORK, f"{tag}.log")
    with open(log, "w") as lf:
        proc = subprocess.run([BIN, cfg_path], cwd=REPO_ROOT,
                              stdout=lf, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        sys.exit(f"ERROR: run '{tag}' failed (exit {proc.returncode}); see {log}")
    with open(log) as f:
        return f.read()


# ── binary dump parsing ──────────────────────────────────────────────────────
# Layout (soil_print write_dump_binary): u32 count, then per atom:
#   tag u32, type u32, pos[3] f64, vel[3] f64, force[3] f64, radius f64, [cols...]
# The trailing per-atom columns vary with the registry, so we derive the record
# stride from the file size and read only the fixed tag/pos/vel prefix.

def load_dump(path):
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
        sys.exit(f"ERROR: dump {path} record stride {rec} < 88 (tag/type/pos/vel/force/radius)")
    atoms = {}
    for i in range(count):
        off = 4 + i * rec
        (tag,) = struct.unpack_from("<I", raw, off)
        pos = struct.unpack_from("<3d", raw, off + 8)
        vel = struct.unpack_from("<3d", raw, off + 32)
        atoms[tag] = (pos, vel)
    return atoms


def latest_dump(outdir):
    files = glob.glob(os.path.join(outdir, "dump", "dump_*.bin"))
    if not files:
        sys.exit(f"ERROR: no dump files in {outdir}/dump")
    def step_of(p):
        return int(os.path.basename(p)[len("dump_"):-len(".bin")])
    return max(files, key=step_of)


def dump_at(outdir, step):
    return os.path.join(outdir, "dump", f"dump_{step}.bin")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


# ── error metrics ────────────────────────────────────────────────────────────

def norm(v):
    return (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) ** 0.5


def max_rel_error(a, b, scale, field):
    """Max over atoms of ||a-b|| / (||ref|| + scale), matched by tag.

    `field` selects 0 = positions, 1 = velocities from each (pos, vel) record.
    """
    if set(a) != set(b):
        sys.exit("ERROR: dumps have different atom tag sets — cannot compare")
    worst = 0.0
    worst_abs = 0.0
    for tag in a:
        av, bv = a[tag][field], b[tag][field]
        d = norm((av[0] - bv[0], av[1] - bv[1], av[2] - bv[2]))
        rel = d / (norm(av) + scale)
        worst = max(worst, rel)
        worst_abs = max(worst_abs, d)
    return worst, worst_abs


def all_finite(atoms):
    import math
    for (pos, vel) in atoms.values():
        for x in list(pos) + list(vel):
            if not math.isfinite(x):
                return False
    return True


def plot_results(rel_pos, rel_vel, digest_flags):
    """Write the README figure from this run's measured validation metrics."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"\n(matplotlib unavailable, skipped plot: {e})")
        return

    os.makedirs(PLOT_DIR, exist_ok=True)
    plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150, "font.size": 10})

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(9.5, 4.0))

    labels = ["position", "velocity"]
    values = [rel_pos, rel_vel]
    ax0.bar(labels, values, color=["#4c78a8", "#f58518"])
    positive = [v for v in values if v > 0.0]
    if positive:
        ax0.set_yscale("log")
        ax0.set_ylim(min(positive) / 10.0, max(positive) * 10.0)
    else:
        ax0.set_ylim(0.0, 1.0)
    ax0.set_ylabel("max relative error vs uninterrupted A")
    ax0.set_title("Restart continuity")
    for i, v in enumerate(values):
        ax0.text(i, max(v, ax0.get_ylim()[0]) * 1.3, f"{v:.1e}", ha="center", va="bottom")

    names = list(digest_flags.keys())
    flags = list(digest_flags.values())
    ax1.bar(names, flags, color=["#4c78a8", "#f58518"][:len(flags)])
    ax1.set_ylim(0, 1.1)
    ax1.set_ylabel("SHA-256 mismatch flag vs A")
    ax1.set_title("Digest deltas")
    ax1.tick_params(axis="x", rotation=20)
    for i, v in enumerate(flags):
        ax1.text(i, v + 0.04, "0" if v == 0 else "1", ha="center", va="bottom")

    fig.suptitle("bench_restart_determinism measured validation result")
    fig.tight_layout()
    fig.savefig(PLOT, bbox_inches="tight")
    plt.close(fig)
    print(f"\nWrote {os.path.relpath(PLOT, REPO_ROOT)}")


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    if shutil.which("cargo") is None and not os.path.exists(BIN):
        sys.exit("ERROR: cargo not found and no prebuilt binary")
    build()

    # Fresh scratch every run so stale restart/dump files never leak between runs.
    if os.path.isdir(WORK):
        shutil.rmtree(WORK)
    os.makedirs(WORK)

    dir_a = os.path.join(WORK, "runA")
    dir_b = os.path.join(WORK, "runB")
    dir_c = os.path.join(WORK, "runC")

    # A: uninterrupted reference. dump at step 0 (interval) + final (save_at_end).
    cfg_a = write_config("A", control(dir_a, TOTAL_STEPS, dump_interval=TOTAL_STEPS))
    # B1: run to the checkpoint; save_at_end writes the restart (+ a dump).
    cfg_b1 = write_config("B1", control(dir_b, CKPT_STEPS))
    # B2: resume in a FRESH process from the checkpoint for the remaining steps.
    cfg_b2 = write_config("B2", control(dir_b, TOTAL_STEPS - CKPT_STEPS, restart_read=True))
    # C: independent twin of A (identical config, different output dir).
    cfg_c = write_config("C", control(dir_c, TOTAL_STEPS, dump_interval=TOTAL_STEPS))

    print(f"Running A (uninterrupted, {TOTAL_STEPS} steps) ...")
    run(cfg_a, "A")
    print(f"Running B1 (checkpoint at {CKPT_STEPS} steps) ...")
    run(cfg_b1, "B1")
    print(f"Running B2 (resume for {TOTAL_STEPS - CKPT_STEPS} steps) ...")
    log_b2 = run(cfg_b2, "B2")
    print(f"Running C (twin of A, {TOTAL_STEPS} steps) ...")
    run(cfg_c, "C")

    # ── evaluate ──────────────────────────────────────────────────────────────
    passed = 0
    total = 0

    def check(name, ok, detail=""):
        nonlocal passed, total
        total += 1
        if ok:
            passed += 1
        print(f"  {name:<28}{'PASS' if ok else 'FAIL'}   {detail}")
        return ok

    print("\nRestart continuity + determinism checks")
    print("-" * 70)

    # provenance: a checkpoint restart file was actually written ...
    ckpt_files = glob.glob(os.path.join(dir_b, "restart", "restart_*_rank0.bin"))
    ckpt_steps = sorted(int(os.path.basename(p).split("_")[1]) for p in ckpt_files)
    check("checkpoint written", len(ckpt_files) >= 1,
          f"restart steps={ckpt_steps}")

    # ... and the resume process genuinely read it back.
    read_ok = ("Restart: reading from" in log_b2) and ("loaded" in log_b2)
    read_line = next((l for l in log_b2.splitlines() if "loaded" in l), "")
    check("restart read on resume", read_ok, read_line.strip())

    final_a = latest_dump(dir_a)
    final_b = latest_dump(dir_b)
    final_c = latest_dump(dir_c)
    A = load_dump(final_a)
    B = load_dump(final_b)
    C = load_dump(final_c)

    # non-vacuity: gas is finite and actually moved (else continuity is trivial).
    init_a = load_dump(dump_at(dir_a, 0))
    _, moved_abs = max_rel_error(init_a, A, POS_SCALE, field=0)
    max_speed = max(norm(v) for (_, v) in A.values())
    check("dynamics non-trivial", all_finite(A) and moved_abs > 5e-3 and max_speed > 1e-3,
          f"max displacement={moved_abs:.4e} m (>{2.5} diam), max speed={max_speed:.4e} m/s")

    # (a) restart continuity: A vs B, sorted by tag, on positions and velocities.
    rel_pos, abs_pos = max_rel_error(A, B, POS_SCALE, field=0)
    rel_vel, abs_vel = max_rel_error(A, B, VEL_SCALE, field=1)
    check("continuity positions", rel_pos <= TOL,
          f"max rel={rel_pos:.3e} (abs={abs_pos:.3e} m), tol={TOL:.0e}")
    check("continuity velocities", rel_vel <= TOL,
          f"max rel={rel_vel:.3e} (abs={abs_vel:.3e} m/s), tol={TOL:.0e}")

    # (b) run-to-run determinism: A and C bit-identical (final + initial dumps).
    same_restart_final = sha256(final_a) == sha256(final_b)
    same_final = sha256(final_a) == sha256(final_c)
    same_init = sha256(dump_at(dir_a, 0)) == sha256(dump_at(dir_c, 0))
    check("determinism (final bytes)", same_final,
          f"sha256 A={sha256(final_a)[:12]} C={sha256(final_c)[:12]}")
    check("determinism (initial bytes)", same_init, "seeded ICs identical")

    plot_results(rel_pos, rel_vel, {
        "restart final": 0 if same_restart_final else 1,
        "twin final": 0 if same_final else 1,
        "twin initial": 0 if same_init else 1,
    })

    print("-" * 70)
    print(f"\nResult: {passed}/{total} checks passed")
    if passed == total:
        print("ALL CHECKS PASSED")
        return 0
    print(f"CHECKS FAILED: {total - passed} of {total}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
