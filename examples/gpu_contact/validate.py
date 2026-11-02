"""Is the GPU contact force within single-precision reach of the CPU?

A settled pile is chaotic, so any change in rounding eventually sends grains elsewhere.
The fair question is whether the GPU differs from double-precision CPU by about as much
as single-precision CPU does, over a dense run with springs under load. Atoms re-sort
every 200 steps, so history surviving a sort is part of what this checks.

    cargo build --release --features gpu --example gpu_contact
    cargo build --release --no-default-features --features gpu,precision-single \
        --example gpu_contact --target-dir target-single
    python3 examples/gpu_contact/validate.py
"""
import csv
import glob
import os
import re
import shutil
import subprocess
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
DOUBLE = os.path.join(ROOT, "target/release/examples/gpu_contact")
SINGLE = os.path.join(ROOT, "target-single/release/examples/gpu_contact")
CHECKPOINTS = [1000, 2000, 5000]


def config(steps, sort_every=200, window=1):
    text = open(os.path.join(HERE, "config.toml")).read()
    # Packed as closely as random insertion manages, so springs load early.
    text = re.sub(r"max = \[0.037, 0.037, 0.115\]", "max = [0.037, 0.037, 0.058]", text)
    text = re.sub(r"sort_every = \d+", f"sort_every = {sort_every}", text)
    text = re.sub(r"thermo = \d+", f"thermo = {max(steps // window // 2, 1)}", text)
    return re.sub(r"steps = \d+", f"steps = {steps // window}", text)


def run(binary, mode, steps, ranks=1, sort_every=200, window=1):
    with tempfile.TemporaryDirectory() as tmp:
        text = config(steps, sort_every, window)
        if ranks > 1:
            text = text.replace("processors_x = 1", f"processors_x = {ranks}")
        path = os.path.join(tmp, "config.toml")
        with open(path, "w") as f:
            f.write(text)
        every = "100" if mode != "resident" else "1"
        env = dict(os.environ, DIRT_CONTACT=mode, DIRT_WRITE_EVERY=every, DIRT_WINDOW=str(window))
        command = [os.path.abspath(binary), "config.toml"]
        if ranks > 1:
            command = ["mpirun", "-n", str(ranks), "--oversubscribe"] + command
        subprocess.run(command, cwd=tmp, env=env, check=True, capture_output=True)
        rows = {}
        for name in glob.glob(os.path.join(tmp, f"{mode}.rank*.csv")):
            with open(name) as f:
                for row in csv.DictReader(f):
                    rows[int(row["tag"])] = [float(row[k]) for k in row if k != "tag"]
        return rows


def worst(a, b):
    assert a.keys() == b.keys(), f"{len(a)} atoms against {len(b)}"
    return {
        "position": max(abs(a[t][i] - b[t][i]) for t in a for i in range(3)),
        "velocity": max(abs(a[t][i] - b[t][i]) for t in a for i in range(3, 6)),
        "spin": max(abs(a[t][i] - b[t][i]) for t in a for i in range(6, 9)),
    }


print(f"{'steps':>6} {'quantity':>9} {'CPU single':>12} {'GPU force':>12} {'whole step':>12}  (each against CPU double)")
for steps in CHECKPOINTS:
    reference = run(DOUBLE, "cpu", steps)
    single = worst(reference, run(SINGLE, "cpu", steps))
    gpu = worst(reference, run(DOUBLE, "gpu", steps))
    resident = worst(reference, run(DOUBLE, "resident", steps, window=100))
    for quantity in ("position", "velocity", "spin"):
        print(f"{steps:>6} {quantity:>9} {single[quantity]:>12.3e} {gpu[quantity]:>12.3e} {resident[quantity]:>12.3e}", flush=True)

print("\nWhole-step windows should change the trajectory only by rounding")
base = run(DOUBLE, "resident", 2000, window=500)
for window in (10, 100):
    diff = worst(base, run(DOUBLE, "resident", 2000, window=window))
    print(f"window {window} against 500: position {diff['position']:.3e}, spin {diff['spin']:.3e}", flush=True)

print("\nDifferent rank counts, against one rank of the same contact force, 2000 steps")
for mode in ("cpu", "gpu"):
    base = run(DOUBLE, mode, 2000)
    for ranks in (2, 4):
        diff = worst(base, run(DOUBLE, mode, 2000, ranks=ranks))
        print(f"{mode} {ranks} ranks: position {diff['position']:.3e}, velocity {diff['velocity']:.3e}, spin {diff['spin']:.3e}", flush=True)
