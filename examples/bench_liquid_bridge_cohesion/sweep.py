#!/usr/bin/env python3
"""Validate DIRT's opt-in Willett pendular liquid-bridge cohesion.

Checks:
  1. single bridge force vs the Willett et al. (2000) closed form,
  2. a small lifted-cylinder heap has a larger static repose angle as bridge
     volume increases, matching the established wet-granular trend.
"""

import csv
import math
import os
import shutil
import subprocess
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")

RADIUS = 0.002
R_EFF = RADIUS / 2.0
VOLUME = 1.0e-11
GAMMA = 0.072
THETA = 0.0
RUPTURE = 1.5e-4
FORCE_TOL = 1.0e-9

REPOSE_CASES = [
    ("dry", 0.0),
    ("wet_low", 3.0e-12),
    ("wet_high", 1.2e-11),
]
REPOSE_MIN_INCREASE_DEG = 2.0


def run(cmd):
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=REPO_ROOT, check=True)


def bridge_ref(separation, volume=VOLUME):
    if separation < 0.0 or separation > RUPTURE or volume <= 0.0:
        return 0.0
    s_hat = separation * math.sqrt(R_EFF / volume)
    return -(2.0 * math.pi * R_EFF * GAMMA * math.cos(THETA)) / (
        1.0 + 1.05 * s_hat + 2.5 * s_hat * s_hat
    )


def start_bridge():
    run([
        "cargo",
        "run",
        "--release",
        "--example",
        "bench_liquid_bridge_cohesion",
        "--no-default-features",
        "--features",
        "precision-double",
        "--",
        os.path.join(SCRIPT_DIR, "config.toml"),
    ])


def start_dry_identity():
    case_root = os.path.join(SWEEP_DIR, "dry_identity")
    os.makedirs(case_root, exist_ok=True)
    base = open(os.path.join(SCRIPT_DIR, "config.toml")).read()
    no_liquid = "\n".join(
        line
        for line in base.splitlines()
        if not line.startswith("liquid_")
        and not line.startswith("velocity_x")
        and "liquid_bridge_model" not in line
    )
    no_liquid = no_liquid.replace(
        'dir = "examples/bench_liquid_bridge_cohesion"',
        f'dir = "{os.path.join(case_root, "default")}"',
    )
    no_liquid = no_liquid.replace("[run]\nsteps = 200", "[run]\nsteps = 50")
    no_liquid = no_liquid.replace(
        'region = { type = "block", min = [0.0040009',
        'velocity_x = 0.02\nregion = { type = "block", min = [0.0040009',
    )
    off = base.replace('liquid_bridge_model = "willett2000"', 'liquid_bridge_model = "off"')
    off = off.replace(
        'dir = "examples/bench_liquid_bridge_cohesion"',
        f'dir = "{os.path.join(case_root, "off")}"',
    )
    off = off.replace("[run]\nsteps = 200", "[run]\nsteps = 50")
    off = off.replace("velocity_x = 0.06", "velocity_x = 0.02")
    configs = [
        (os.path.join(case_root, "default.toml"), no_liquid),
        (os.path.join(case_root, "off.toml"), off),
    ]
    for cfg, text in configs:
        with open(cfg, "w") as f:
            f.write(text)
        run([
            "cargo",
            "run",
            "--release",
            "--example",
            "bench_liquid_bridge_cohesion",
            "--no-default-features",
            "--features",
            "precision-double",
            "--",
            cfg,
        ])


def start_invalid_model_rejected():
    case_root = os.path.join(SWEEP_DIR, "invalid_model")
    os.makedirs(case_root, exist_ok=True)
    base = open(os.path.join(SCRIPT_DIR, "config.toml")).read()
    bad = base.replace('liquid_bridge_model = "willett2000"', 'liquid_bridge_model = "willet2000"')
    bad = bad.replace(
        'dir = "examples/bench_liquid_bridge_cohesion"',
        f'dir = "{case_root}"',
    )
    path = os.path.join(case_root, "config.toml")
    with open(path, "w") as f:
        f.write(bad)
    cmd = [
        "cargo",
        "run",
        "--release",
        "--example",
        "bench_liquid_bridge_cohesion",
        "--no-default-features",
        "--features",
        "precision-double",
        "--",
        path,
    ]
    print("+", " ".join(cmd), flush=True)
    result = subprocess.run(cmd, cwd=REPO_ROOT, text=True, capture_output=True)
    combined = result.stdout + result.stderr
    if result.returncode == 0:
        raise SystemExit("invalid liquid_bridge_model typo unexpectedly exited 0")
    if "liquid_bridge_model" not in combined or "willet2000" not in combined:
        raise SystemExit(
            "invalid liquid_bridge_model error must name the field and bad value; "
            f"got returncode={result.returncode}, output:\n{combined}"
        )
    print("invalid_model: rejected bad liquid_bridge_model -> PASS")


REPOSE_TEMPLATE = """\
[comm]
processors_x = 1
processors_y = 1
processors_z = 1

[domain]
x_low = -0.06
x_high = 0.06
y_low = -0.06
y_high = 0.06
z_low = 0.0
z_high = 0.12
boundary_x = "fixed"
boundary_y = "fixed"
boundary_z = "fixed"

[neighbor]
skin_fraction = 1.2
bin_size = 0.007
every = 1

[gravity]
gz = -9.81

[dem]
contact_model = "hertz"
rolling_model = "sds"
liquid_bridge_model = "{liquid_model}"

[[dem.materials]]
name = "glass"
youngs_mod = 1.0e7
poisson_ratio = 0.25
restitution = 0.4
friction = 0.3
rolling_friction = 0.1
rolling_stiffness = 1.0e-2
rolling_damping = 1.0e-6
liquid_bridge_volume = {volume:.8e}
liquid_surface_tension = 0.072
liquid_contact_angle = 0.0
liquid_rupture_distance = {rupture:.8e}

[[wall]]
type = "cylinder"
axis = "z"
center = [0.0, 0.0]
radius = 0.018
lo = 0.0
hi = 0.10
inside = true
material = "glass"
name = "cylinder"

[[wall]]
type = "plane"
point_x = 0.0
point_y = 0.0
point_z = 0.0
normal_x = 0.0
normal_y = 0.0
normal_z = 1.0
material = "glass"

[[wall]]
type = "cylinder"
axis = "z"
center = [0.0, 0.0]
radius = 0.055
lo = 0.0
hi = 0.12
inside = true
material = "glass"

[[particles.insert]]
material = "glass"
count = 260
radius = 0.0025
density = 2500.0
velocity_z = -0.1
region = {{ type = "cylinder", center = [0.0, 0.0], radius = 0.016, axis = "z", lo = 0.004, hi = 0.085 }}

[output]
dir = "{output_dir}"

[[run]]
name = "fill"
steps = 30000
thermo = 10000

[[run]]
name = "lift"
steps = 70000
thermo = 10000
"""


def write_repose_config(name, volume):
    case_dir = os.path.join(SWEEP_DIR, name)
    os.makedirs(case_dir, exist_ok=True)
    model = "off" if volume <= 0.0 else "willett2000"
    text = REPOSE_TEMPLATE.format(
        liquid_model=model,
        volume=volume,
        rupture=1.5 * volume ** (1.0 / 3.0) if volume > 0.0 else 0.0,
        output_dir=case_dir,
    )
    path = os.path.join(case_dir, "config.toml")
    with open(path, "w") as f:
        f.write(text)
    return path


def start_repose():
    os.makedirs(SWEEP_DIR, exist_ok=True)
    for name, volume in REPOSE_CASES:
        cfg = write_repose_config(name, volume)
        run([
            "cargo",
            "run",
            "--release",
            "--example",
            "bench_angle_of_repose",
            "--no-default-features",
            "--features",
            "precision-double",
            "--",
            cfg,
        ])


def fit_repose_angle(path):
    rows = []
    with open(path) as f:
        for row in csv.DictReader(f):
            x = float(row["x"])
            y = float(row["y"])
            z = float(row["z"])
            radius = float(row["radius"])
            rows.append((math.hypot(x, y), z + radius))
    if len(rows) < 50:
        raise RuntimeError(f"too few particles in {path}")
    rmax = max(r for r, _ in rows)
    bins = 12
    surface = []
    for b in range(bins):
        lo = rmax * b / bins
        hi = rmax * (b + 1) / bins
        zs = [z for r, z in rows if lo <= r < hi]
        if len(zs) >= 3:
            surface.append((0.5 * (lo + hi), max(zs)))
    fit = [(r, z) for r, z in surface if 0.15 * rmax <= r <= 0.85 * rmax]
    n = len(fit)
    sx = sum(r for r, _ in fit)
    sy = sum(z for _, z in fit)
    sxx = sum(r * r for r, _ in fit)
    sxy = sum(r * z for r, z in fit)
    slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    return math.degrees(math.atan(max(0.0, -slope))), surface


def graph():
    os.makedirs(PLOT_DIR, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    trace_path = os.path.join(DATA_DIR, "bridge_trace.csv")
    samples = []
    with open(trace_path) as f:
        for row in csv.DictReader(f):
            s = float(row["separation"])
            f_n = float(row["normal_force"])
            if 0.0 <= s <= RUPTURE and f_n < 0.0:
                samples.append((s, f_n, bridge_ref(s)))
    if not samples:
        raise SystemExit(
            "no tensile liquid-bridge samples found; the DEM cutoff must include "
            "bridge-only separations with skin_fraction=1.0"
        )
    max_rel = max(abs((f - ref) / ref) for s, f, ref in samples)
    force_pass = max_rel <= FORCE_TOL and max(s for s, _, _ in samples) > 0.95 * RUPTURE

    with open(os.path.join(DATA_DIR, "bridge_force_summary.csv"), "w") as f:
        f.write("max_relative_error,tolerance,pass\n")
        f.write(f"{max_rel:.12e},{FORCE_TOL:.12e},{force_pass}\n")

    dry_a = os.path.join(SWEEP_DIR, "dry_identity", "default", "data", "bridge_trace.csv")
    dry_b = os.path.join(SWEEP_DIR, "dry_identity", "off", "data", "bridge_trace.csv")
    dry_identical = open(dry_a, "rb").read() == open(dry_b, "rb").read()
    with open(os.path.join(DATA_DIR, "dry_identity_summary.csv"), "w") as f:
        f.write("default_trace,off_trace,byte_identical\n")
        f.write(f"{dry_a},{dry_b},{dry_identical}\n")

    s_vals = [s * 1e6 for s, _, _ in samples]
    measured = [-f * 1e3 for _, f, _ in samples]
    ref = [-r * 1e3 for _, _, r in samples]
    plt.figure(figsize=(6.0, 4.0))
    plt.plot(s_vals, ref, "k-", label="Willett 2000 closed form")
    plt.plot(s_vals, measured, "o", ms=3, label="DIRT")
    plt.axvline(RUPTURE * 1e6, color="0.5", ls="--", label="rupture")
    plt.xlabel("surface separation (micron)")
    plt.ylabel("tensile force (mN)")
    plt.title("Pendular liquid-bridge force")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "bridge_force.png"), dpi=180)
    plt.close()

    angles = []
    surfaces = {}
    for name, volume in REPOSE_CASES:
        path = os.path.join(SWEEP_DIR, name, "data", "repose_results.csv")
        angle, surface = fit_repose_angle(path)
        angles.append((name, volume, angle))
        surfaces[name] = surface
    with open(os.path.join(DATA_DIR, "repose_trend.csv"), "w") as f:
        f.write("case,bridge_volume,angle_deg\n")
        for name, volume, angle in angles:
            f.write(f"{name},{volume:.12e},{angle:.6f}\n")

    dry = angles[0][2]
    wet_high = angles[-1][2]
    monotone = all(angles[i + 1][2] >= angles[i][2] - 0.5 for i in range(len(angles) - 1))
    repose_pass = monotone and wet_high - dry >= REPOSE_MIN_INCREASE_DEG

    plt.figure(figsize=(6.0, 4.0))
    xs = list(range(len(angles)))
    plt.plot(xs, [a for _, _, a in angles], "o-", label="DIRT mini heap")
    plt.xticks(xs, [f"{name}\n{volume:.1e} m^3" for name, volume, _ in angles])
    plt.xlabel("case and liquid bridge volume per contact")
    plt.ylabel("static repose angle (deg)")
    plt.title("Wet heap trend")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "wet_repose_trend.png"), dpi=180)
    plt.close()

    print(f"bridge_force: max_rel={max_rel:.3e} tol={FORCE_TOL:.1e} -> {'PASS' if force_pass else 'FAIL'}")
    print(f"dry_identity: byte_identical={dry_identical} -> {'PASS' if dry_identical else 'FAIL'}")
    print("repose_angles:", ", ".join(f"{n}={a:.2f} deg" for n, _, a in angles))
    print(
        f"wet_repose_trend: high-dry={wet_high-dry:.2f} deg "
        f"gate={REPOSE_MIN_INCREASE_DEG:.2f} deg -> {'PASS' if repose_pass else 'FAIL'}"
    )
    if not (force_pass and dry_identical and repose_pass):
        raise SystemExit(1)
    print("ALL CHECKS PASSED")


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    os.makedirs(DATA_DIR, exist_ok=True)
    if cmd in ("all", "start"):
        start_bridge()
        start_dry_identity()
        start_invalid_model_rejected()
        start_repose()
    if cmd in ("all", "graph"):
        graph()


if __name__ == "__main__":
    main()
