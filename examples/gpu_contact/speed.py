"""Steps per second for the CPU and GPU contact force as the system grows.

A dense periodic granular gas (solid fraction 0.5, the regime in perf_mpi_scaling), where
contact is most of the work. Cases alternate within each size to even out background load.

    cargo build --release --features gpu --example gpu_contact
    python3 examples/gpu_contact/speed.py
"""
import json
import os
import random
import re
import subprocess
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
BINARY = os.path.abspath(os.path.join(HERE, "..", "..", "target/release/examples/gpu_contact"))
SIDES = [20, 30, 40, 60]
# (contact force, MPI ranks)
CASES = [("cpu", 1), ("cpu", 16), ("gpu", 1), ("resident", 1)]
WINDOW = 100
RADIUS, DENSITY, PHI, SIGMA = 0.0008, 2500.0, 0.50, 0.5
EVERY, STEPS, REPEATS = 10, 2000, 2


def setup(tmp, side, ranks, ticks):
    diameter = 2 * RADIUS
    spacing = diameter * (0.5235987756 / PHI) ** (1 / 3)
    box = side * spacing
    rng = random.Random(42)
    with open(os.path.join(tmp, "init.csv"), "w") as f:
        f.write("x,y,z,radius,vx,vy,vz\n")
        for i in range(side):
            for j in range(side):
                for k in range(side):
                    v = [rng.gauss(0, SIGMA) for _ in range(3)]
                    f.write(f"{(i + .5) * spacing:.8e},{(j + .5) * spacing:.8e},{(k + .5) * spacing:.8e},"
                            f"{RADIUS},{v[0]:.6e},{v[1]:.6e},{v[2]:.6e}\n")
    px = {1: (1, 1, 1), 16: (4, 2, 2)}[ranks]
    axes = "".join(f"{a}_low=0.0\n{a}_high={box:.8e}\nboundary_{a}=\"periodic\"\n" for a in "xyz")
    with open(os.path.join(tmp, "config.toml"), "w") as f:
        f.write(f"""[comm]
processors_x={px[0]}
processors_y={px[1]}
processors_z={px[2]}
[domain]
{axes}[neighbor]
skin_fraction=1.1
every={EVERY}
check=false
bin_size={2 * diameter:.6e}
[gravity]
gx=0.0
gy=0.0
gz=0.0
[dem]
contact_model="hertz"
[[dem.materials]]
name="glass"
youngs_mod=5.0e7
poisson_ratio=0.3
restitution=0.9
friction=0.3
[[particles.insert]]
source="file"
format="csv"
file="init.csv"
material="glass"
radius={RADIUS}
columns={{ x=0, y=1, z=2, radius=3, vx=4, vy=5, vz=6 }}
density={DENSITY}
[output]
dir="."
[[run]]
name="run"
steps={ticks}
thermo={ticks // 2}
""")


def steps_per_second(mode, ranks, side):
    # The whole-step path counts one tick per window, and steps by window.
    per_tick = WINDOW if mode == "resident" else 1
    ticks = max(STEPS // per_tick, 20)
    with tempfile.TemporaryDirectory() as tmp:
        setup(tmp, side, ranks, ticks)
        command = [BINARY, "config.toml"]
        if ranks > 1:
            command = ["mpirun", "-n", str(ranks), "--oversubscribe"] + command
        out = subprocess.run(command, cwd=tmp, capture_output=True, text=True, check=True,
                             env=dict(os.environ, DIRT_CONTACT=mode, DIRT_WINDOW=str(WINDOW))).stdout
        # Step, atoms, KE, neighbors, walltime, steps per second, from the timed window.
        rows = re.findall(rf"^{ticks // 2}\s+\d+\s+\S+\s+\d+\s+\S+\s+(\S+)\s*$", out, flags=re.M)
        return float(rows[-1]) * per_tick


results = []
print(f"{'atoms':>8} " + " ".join(f"{m} x{r}".rjust(13) for m, r in CASES) + "   whole/1 core   whole/16 ranks", flush=True)
for side in SIDES:
    best = {c: 0.0 for c in CASES}
    for _ in range(REPEATS):
        for case in CASES:
            best[case] = max(best[case], steps_per_second(*case, side))
    whole, one, many = best[("resident", 1)], best[("cpu", 1)], best[("cpu", 16)]
    results.append({"atoms": side ** 3, **{f"{m}_{r}": best[(m, r)] for m, r in CASES}})
    print(f"{side ** 3:>8} " + " ".join(f"{best[c]:>13.1f}" for c in CASES) + f"   {whole / one:>12.2f}   {whole / many:>14.2f}", flush=True)
with open("gpu_speed.json", "w") as f:
    json.dump(results, f, indent=1)
