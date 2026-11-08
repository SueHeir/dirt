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


# A bonded block at rest, hit by one free grain. Nothing is in contact at the first step, and
# bonds start unstretched, so the CPU's first half kick loses nothing here either.
def bonded_config(steps, window, breakage):
    text = open(os.path.join(HERE, "config.toml")).read()
    text = re.sub(r"gz = -9.81", "gz = 0.0", text)
    text = re.sub(r"\[\[particles.insert\]\].*?(?=\[\[run\]\])", """[[wall]]
type = "plane"
point_x = 0.0
point_y = 0.0
point_z = 0.11
normal_x = 0.0
normal_y = 0.0
normal_z = -1.0
material = "glass"

[[particles.insert]]
source = "file"
format = "csv"
file = "init.csv"
material = "glass"
radius = 0.0015
density = 2500.0
columns = { x = 0, y = 1, z = 2, radius = 3, vx = 4, vy = 5, vz = 6 }

[bonds]
auto_bond = true
bond_tolerance = 1.001
bond_radius_ratio = 0.5
youngs_modulus = 1.0e7
shear_modulus = 4.0e6
beta_normal = 0.2
beta_shear = 0.2
beta_twist = 0.2
beta_bending = 0.2

""" + breakage + "\n", text, flags=re.S)
    text = re.sub(r"thermo = \d+", f"thermo = {max(steps // window // 2, 1)}", text)
    return re.sub(r"steps = \d+", f"steps = {steps // window}", text)


def bonded_init(path):
    r = 0.0015
    with open(path, "w") as f:
        f.write("x,y,z,radius,vx,vy,vz\n")
        for i in range(6):
            for j in range(6):
                for k in range(6):
                    f.write(f"{0.008 + 2 * r * i},{0.008 + 2 * r * j},{0.03 + 2 * r * k},{r},0,0,0\n")
        # The free grain, aimed off-centre so that it twists and bends the block.
        f.write(f"{0.0165},{0.0175},{0.02},{r},0.1,0.2,3.0\n")


def bonded_run(binary, mode, steps, window, breakage):
    with tempfile.TemporaryDirectory() as tmp:
        with open(os.path.join(tmp, "config.toml"), "w") as f:
            f.write(bonded_config(steps, window, breakage))
        bonded_init(os.path.join(tmp, "init.csv"))
        env = dict(os.environ, DIRT_CONTACT=mode, DIRT_BONDS="1", DIRT_WINDOW=str(window), DIRT_WRITE_EVERY="1")
        done = subprocess.run([os.path.abspath(binary), "config.toml"], cwd=tmp, env=env, capture_output=True, text=True)
        if done.returncode != 0:
            at = done.stderr.find("panicked")
            raise RuntimeError(f"{mode}: {done.stderr[at:at + 600]}")
        rows = {}
        for name in glob.glob(os.path.join(tmp, f"{mode}.rank*.csv")):
            with open(name) as f:
                for row in csv.DictReader(f):
                    rows[int(row["tag"])] = [float(row[k]) for k in row if k != "tag"]
        return rows


if os.path.exists(F64):
    cases = [("elastic bonds", ""), ("breaking bonds", """[bonds.breakage]
kind = "combined_stress"
tensile = { kind = "constant", value = 2.0e3 }
shear = { kind = "constant", value = 2.0e3 }
""")]
    print("\nA bonded block hit by a free grain: f64 kernels on the CPU runtime against CPU double")
    moved = worst(bonded_run(DOUBLE, "cpu", 2000, 1, ""), bonded_run(DOUBLE, "cpu", 2000, 1, cases[1][1]))
    print(f"Breakage changes the CPU run by {moved['position']:.2e} in position at 2000 steps, so it did happen")
    for name, breakage in cases:
        for steps in (1000, 2000, 4000):
            reference = bonded_run(DOUBLE, "cpu", steps, 1, breakage)
            force = worst(reference, bonded_run(F64, "gpu", steps, 1, breakage))
            whole = worst(reference, bonded_run(F64, "resident", steps, 50, breakage))
            print(f"{name:>15} {steps:>6} GPU force: position {force['position']:.2e}, spin {force['spin']:.2e}; whole step: position {whole['position']:.2e}, spin {whole['spin']:.2e}", flush=True)


# Spring-based rolling and twisting, which are well conditioned, so trajectories can be compared.
def rolling_run(binary, mode, steps, window, rolling, twisting, mutate=lambda text: text):
    global materials_config
    plain = materials_config
    def config_with_models(s, w=1):
        text = plain(s, w)
        text = text.replace('contact_model = "hertz"', f'contact_model = "hertz"\nrolling_model = "{rolling}"\ntwisting_model = "{twisting}"')
        sds_extra = """rolling_friction = {r}
rolling_stiffness = 1.0e-3
rolling_damping = 1.0e-6
twisting_friction = {t}
twisting_stiffness = 1.0e-3
twisting_damping = 1.0e-6
"""
        text = text.replace("friction = 0.5\n", "friction = 0.5\n" + sds_extra.format(r=0.1, t=0.05), 1)
        return mutate(text.replace("friction = 0.2\n", "friction = 0.2\n" + sds_extra.format(r=0.05, t=0.02), 1))
    materials_config = config_with_models
    try:
        return materials_run(binary, mode, steps, window=window)
    finally:
        materials_config = plain


if os.path.exists(F64):
    print("\nSpring-based rolling and twisting, two materials, walls: f64 kernels against CPU double")
    for rolling, twisting in (("sds", "sds"), ("sds", "marshall")):
        for steps in (1000, 2000, 4000):
            reference = rolling_run(DOUBLE, "cpu", steps, 1, rolling, twisting)
            force = worst(reference, rolling_run(F64, "gpu", steps, 1, rolling, twisting))
            whole = worst(reference, rolling_run(F64, "resident", steps, 100, rolling, twisting))
            print(f"{rolling}/{twisting} {steps:>6} GPU force: position {force['position']:.2e}, spin {force['spin']:.2e}; whole step: position {whole['position']:.2e}, spin {whole['spin']:.2e}", flush=True)
