"""Do the GPU kernels agree with the CPU force?

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

# Several materials and walls, with no gravity. With no
# gravity and no contact at the start, the CPU's first half kick loses nothing, so the whole
# step can be compared directly. Build with --features gpu-cpu-f64 --target-dir target-f64 to
# run the f64 kernels on CubeCL's CPU runtime, and this should reach rounding error.
F64 = os.path.join(ROOT, "target-f64/release/examples/gpu_contact")


def materials_config(steps, window=1):
    text = open(os.path.join(HERE, "config.toml")).read()
    text = re.sub(r"gz = -9.81", "gz = 0.0", text)
    text = re.sub(r"\[\[particles.insert\]\].*?(?=\[\[run\]\])", """[[particles.insert]]
material = "glass"
count = 250
radius = 0.0018
density = 2500.0
velocity = 0.3
seed = 7
region = { type = "block", min = [0.003, 0.003, 0.003], max = [0.037, 0.037, 0.026] }

[[particles.insert]]
material = "steel"
count = 200
radius = 0.0012
density = 7800.0
velocity = 0.3
seed = 11
region = { type = "block", min = [0.003, 0.003, 0.03], max = [0.037, 0.037, 0.05] }

""", text, flags=re.S)
    # The blocks do not touch, because a contact at the first step would feel the CPU's first
    # half kick. Rolling and twisting friction stay off here. Their constant-torque models are singular
    # when two grains spin alike: the torque keeps full size, with a direction set by rounding
    # noise, so trajectories diverge from the first such contact on any arithmetic change.
    text = text.replace('friction = 0.5\n', """friction = 0.5

[[dem.materials]]
name = "steel"
youngs_mod = 3.0e7
poisson_ratio = 0.3
restitution = 0.7
friction = 0.2
""", 1)
    text = text.replace('z_high = 0.12', 'z_high = 0.06')
    # A lid, since grains now fly upward too. Atoms leaving the box would change the count.
    text = text.replace('[[particles.insert]]', '''[[wall]]
type = "plane"
point_x = 0.0
point_y = 0.0
point_z = 0.055
normal_x = 0.0
normal_y = 0.0
normal_z = -1.0
material = "steel"

[[particles.insert]]''', 1)
    text = re.sub(r"thermo = \d+", f"thermo = {max(steps // window // 2, 1)}", text)
    return re.sub(r"steps = \d+", f"steps = {steps // window}", text)


def materials_run(binary, mode, steps, window=1):
    global config
    original = config
    config = lambda s, sort, w: materials_config(s, w)
    try:
        return run(binary, mode, steps, window=window)
    finally:
        config = original


if os.path.exists(F64):
    print("\nSeveral materials and walls, no gravity: f64 kernels on the CPU runtime")
    for steps in (1000, 2000, 4000):
        reference = materials_run(DOUBLE, "cpu", steps)
        force = worst(reference, materials_run(F64, "gpu", steps))
        whole = worst(reference, materials_run(F64, "resident", steps, window=100))
        print(f"{steps:>6} GPU force: position {force['position']:.2e}, spin {force['spin']:.2e}; whole step: position {whole['position']:.2e}, spin {whole['spin']:.2e}", flush=True)
