#!/usr/bin/env python3
"""
analyze_discharge.py — post-process a hopper-discharge run and report the steady
mass-flow rate W (the "discharge rate").

This is the ANALYZE step of a config -> run -> analyze workflow. It reads the
standard per-atom dump snapshots written by the generic `run` driver (it needs
NO custom recorder in the simulation), counts the cumulative mass of grains that
have dropped below the orifice plane vs. time, fits the steady slope of that
curve, and prints the discharge rate. It also saves a figure of the discharge
curve with the fitted steady line.

Usage:
    python3 examples/run/analyze_discharge.py examples/run/hopper_discharge.toml
    python3 examples/run/analyze_discharge.py examples/run/hopper_discharge.toml --orifice-z 0.045

The single positional argument is the CONFIG you ran. The script reads the run's
output directory, timestep, and slot depth from that config, so the same command
works for any hopper config you adapt — change the config, re-run the driver,
re-run this script.

Dump CSV columns (written by the driver's `[dump]` section):
    tag,type,x,y,z,vx,vy,vz,fx,fy,fz,radius
"""

import argparse
import glob
import math
import os
import sys
import tomllib

import numpy as np


def load_config(config_path):
    """Read the run's output dir, timestep, orifice guess, and slot depth."""
    with open(config_path, "rb") as f:
        cfg = tomllib.load(f)

    # Output directory: [output].dir if present, else the config's own folder.
    out_dir = cfg.get("output", {}).get("dir")
    if out_dir is None:
        out_dir = os.path.splitext(config_path)[0]

    # Timestep from the (single) run stage.
    runs = cfg.get("run", [])
    dt = runs[0]["dt"] if runs else None

    # Slot depth = periodic y-extent (for the per-unit-depth flow rate).
    dom = cfg.get("domain", {})
    depth = None
    if "y_high" in dom and "y_low" in dom:
        depth = float(dom["y_high"]) - float(dom["y_low"])

    return out_dir, dt, depth


def read_frame(path):
    """Return (tag, z, radius) arrays for one dump CSV frame."""
    arr = np.genfromtxt(path, delimiter=",", names=True)
    # A frame with a single atom parses 0-d; force at least 1-d.
    tag = np.atleast_1d(arr["tag"])
    z = np.atleast_1d(arr["z"])
    radius = np.atleast_1d(arr["radius"])
    return tag, z, radius


def main():
    ap = argparse.ArgumentParser(description="Fit the steady hopper discharge rate W.")
    ap.add_argument("config", help="the config.toml you ran (e.g. examples/run/hopper_discharge.toml)")
    ap.add_argument("--orifice-z", type=float, default=0.045,
                    help="a grain is 'discharged' once its centre drops below this z (m). Default 0.045.")
    ap.add_argument("--density", type=float, default=2500.0,
                    help="grain density (kg/m^3) for the grain mass. Default 2500 (glass).")
    ap.add_argument("--lo-frac", type=float, default=0.1,
                    help="start of the steady-fit window as a fraction of total discharged mass. Default 0.1.")
    ap.add_argument("--hi-frac", type=float, default=0.9,
                    help="end of the steady-fit window. Default 0.9.")
    args = ap.parse_args()

    out_dir, dt, depth = load_config(args.config)
    if dt is None:
        sys.exit("error: could not read [[run]].dt from the config")

    dump_dir = os.path.join(out_dir, "dump")
    files = glob.glob(os.path.join(dump_dir, "dump_*.csv"))
    if not files:
        sys.exit(f"error: no dump/dump_*.csv frames under {out_dir} — did the run write [dump]?")

    # Order frames by step (parsed from dump_<step>.csv).
    def step_of(p):
        return int(os.path.basename(p)[len("dump_"):-len(".csv")])
    files.sort(key=step_of)

    # Grain mass per tag from the FIRST frame (robust to polydispersity and to
    # grains being deleted as they fall past the domain floor).
    tag0, z0, r0 = read_frame(files[0])
    mass_of = {int(t): args.density * (4.0 / 3.0) * math.pi * r**3
               for t, r in zip(tag0, r0)}
    total_mass = sum(mass_of.values())
    n0 = len(tag0)

    times, discharged = [], []
    for p in files:
        step = step_of(p)
        tag, z, _r = read_frame(p)
        # Mass still ABOVE the orifice = grains not yet discharged. Everything
        # else (below the plane, or already deleted at the floor) has discharged.
        above = sum(mass_of.get(int(t), 0.0) for t, zz in zip(tag, z) if zz >= args.orifice_z)
        times.append(step * dt)
        discharged.append(total_mass - above)

    times = np.array(times)
    discharged = np.array(discharged)

    # Steady window: the 10%-90% portion of the total discharged mass, excluding
    # the start-up transient and the empty-out tail.
    final = discharged[-1] if discharged[-1] > 0 else discharged.max()
    lo, hi = args.lo_frac * final, args.hi_frac * final
    mask = (discharged >= lo) & (discharged <= hi)
    if mask.sum() < 2:
        sys.exit("error: steady window has < 2 points — run longer or widen --lo/--hi-frac")

    # Linear fit; slope is the mass-flow rate W (kg/s).
    slope, intercept = np.polyfit(times[mask], discharged[mask], 1)
    fit = slope * times[mask] + intercept
    ss_res = float(np.sum((discharged[mask] - fit) ** 2))
    ss_tot = float(np.sum((discharged[mask] - discharged[mask].mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    W = slope
    print(f"frames                : {len(files)}  ({n0} grains initially, "
          f"total bed mass {total_mass*1e3:.2f} g)")
    print(f"steady-fit window     : {mask.sum()} frames "
          f"({times[mask][0]:.3f}-{times[mask][-1]:.3f} s)")
    print(f"discharge rate W      : {W*1e3:.3f} g/s   ({W:.4e} kg/s)")
    if depth:
        print(f"  per unit slot depth : {W/depth*1e3:.3f} g/(s·m)   "
              f"(slot depth {depth*1e3:.1f} mm)")
    print(f"fit quality R^2       : {r2:.4f}")

    # Figure: cumulative discharged mass vs time + the fitted steady line.
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        ax.plot(times, discharged * 1e3, "o-", ms=3, color="#1f77b4",
                label="cumulative discharged mass")
        ax.plot(times[mask], fit * 1e3, "-", lw=2.2, color="#d62728",
                label=f"steady fit: W = {W*1e3:.2f} g/s  (R²={r2:.3f})")
        ax.set_xlabel("time (s)")
        ax.set_ylabel("discharged mass (g)")
        ax.set_title("Hopper discharge — cumulative mass vs. time")
        ax.legend(loc="lower right", fontsize=9)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        out_png = os.path.join(out_dir, "discharge_curve.png")
        fig.savefig(out_png, dpi=120)
        print(f"figure                : {out_png}")
    except Exception as e:  # matplotlib optional — the W value is the result
        print(f"(figure skipped: {e})")


if __name__ == "__main__":
    main()
