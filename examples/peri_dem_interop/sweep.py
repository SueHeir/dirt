#!/usr/bin/env python3
"""
peri_dem_interop validation driver.

Builds and runs the peridynamics(pond)+DEM(dirt) same-substrate interop example
and checks the physics that the acceptance criterion demands:

  1. HARD GATE — mass & momentum conserved across the peri->DEM transition.
     The example itself prints `CONSERVATION: PASS` only when the max relative
     momentum drift and mass drift over the whole run are < 1e-9. We require it.
  2. Fracture actually happened (peri bonds broke, damage -> 1) — otherwise the
     conservation check is trivial. We require peak damage ~1 and the surviving
     bond count to collapse to a small fraction of the reference family.
  3. The fragments interact via DEM contact through the shared neighbour list —
     the DEM-contact count (overlapping, non-peri-excluded pairs) must rise from
     0 (all in-horizon pairs peri-bonded) to a non-trivial number after impact.

Only the conservation gate is a physics-tolerance check; it is NOT weakenable
here (the tolerance lives in the example and comes from round-off, not a fit).
Checks 2–3 are guards that the run is a genuine mixed peri->DEM fragmentation.

Usage:
    python3 examples/peri_dem_interop/sweep.py            # build + run + validate
    python3 examples/peri_dem_interop/sweep.py start      # build + run -> log
    python3 examples/peri_dem_interop/sweep.py graph      # validate last log
"""

import os
import re
import sys
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "peri_dem_interop"
CONFIG = os.path.join("examples", EXAMPLE, "config.toml")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
LOG = os.path.join(DATA_DIR, "run.log")

CARGO_FLAGS = ["--no-default-features", "--features", "precision-double"]

STEP_RE = re.compile(
    r"step\s+(\d+)\s+KE=(\S+)\s+J\s+\|p\|=(\S+)\s+bonds=\s*(\d+)\s+"
    r"dmax=(\S+)\s+.*DEMcontacts=\s*(\d+)"
)


def run(cmd, **kw):
    print("+", " ".join(cmd))
    return subprocess.run(cmd, cwd=REPO_ROOT, check=True, **kw)


def start():
    os.makedirs(DATA_DIR, exist_ok=True)
    run(["cargo", "build", "--release", "--example", EXAMPLE, *CARGO_FLAGS])
    with open(LOG, "w") as fh:
        run(
            ["cargo", "run", "--release", "--example", EXAMPLE, *CARGO_FLAGS, "--", CONFIG],
            stdout=fh,
            stderr=subprocess.STDOUT,
        )
    print(f"wrote {LOG}")


def graph():
    with open(LOG) as fh:
        text = fh.read()

    steps = [
        (int(m[1]), float(m[3]), int(m[4]), float(m[5]), int(m[6]))
        for m in (STEP_RE.match(line) for line in text.splitlines())
        if m
    ]
    if not steps:
        print("CHECKS FAILED: no step diagnostics found in log")
        return 1

    bonds0 = steps[0][2]
    bonds_final = steps[-1][2]
    peak_damage = max(s[3] for s in steps)
    max_dem = max(s[4] for s in steps)
    conserved = "CONSERVATION: PASS" in text

    checks = [
        ("mass & momentum conserved (example gate)", conserved),
        ("peri bonds broke (fracture occurred)", bonds_final < 0.10 * bonds0),
        ("damage reached ~1 (fragments separated)", peak_damage >= 0.99),
        ("fragments interact via DEM contact", max_dem >= 8),
    ]

    npass = sum(1 for _, ok in checks for _ in [ok] if ok)
    print("\n── peri_dem_interop validation ──")
    print(f"  reference bonds        = {bonds0}")
    print(f"  surviving bonds        = {bonds_final}  ({100*bonds_final/bonds0:.1f}% of reference)")
    print(f"  peak damage            = {peak_damage:.3f}")
    print(f"  peak DEM contacts      = {max_dem}")
    for name, ok in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")

    if npass == len(checks):
        print(f"\n{npass}/{len(checks)} checks passed")
        print("ALL CHECKS PASSED")
        return 0
    print(f"\n{npass}/{len(checks)} checks passed")
    print("CHECKS FAILED")
    return 1


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    if cmd in ("start", "run"):
        start()
        return 0
    if cmd == "graph":
        return graph()
    if cmd == "all":
        start()
        return graph()
    print(f"unknown command: {cmd}")
    return 2


if __name__ == "__main__":
    sys.exit(main())
