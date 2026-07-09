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
    python3 examples/peri_dem_interop/sweep.py graph      # validate + plot last log
"""

import os
import re
import sys
import subprocess
import math
from html import escape

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "peri_dem_interop"
CONFIG = os.path.join("examples", EXAMPLE, "config.toml")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
LOG = os.path.join(DATA_DIR, "run.log")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
PLOT = os.path.join(PLOT_DIR, "peri_dem_transition_validation.svg")

CARGO_FLAGS = ["--no-default-features", "--features", "precision-double"]

STEP_RE = re.compile(
    r"step\s+(\d+)\s+KE=(\S+)\s+J\s+\|p\|=(\S+)\s+bonds=\s*(\d+)\s+"
    r"dmax=(\S+)\s+.*DEMcontacts=\s*(\d+)"
)
REL_RE = re.compile(r"max (mass|momentum) drift\s+=\s+\S+\s+\(rel\s+(\S+)\)")
TOL_RE = re.compile(r"tolerance \(relative\)\s+=\s+(\S+)")


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


def write_plot(steps, rel_mass, rel_p, tol):
    os.makedirs(PLOT_DIR, exist_ok=True)

    xs = [s[0] for s in steps]
    bonds = [s[2] for s in steps]
    damage = [s[3] for s in steps]
    contacts = [s[4] for s in steps]

    width, height = 1200, 520
    left = (70, 80, 490, 380)
    right = (590, 80, 520, 380)

    def sx(step):
        xmin, xmax = min(xs), max(xs)
        return right[0] + (step - xmin) / (xmax - xmin) * right[2]

    def sy_count(value):
        ymax = max(max(bonds), 1.0) * 1.08
        return right[1] + right[3] - value / ymax * right[3]

    def sy_contact(value):
        ymax = max(max(contacts), 1.0) * 1.18
        return right[1] + right[3] - value / ymax * right[3]

    def sy_damage(value):
        return right[1] + right[3] - value / 1.05 * right[3]

    def poly(points):
        return " ".join(f"{x:.1f},{y:.1f}" for x, y in points)

    def line(x1, y1, x2, y2, color, width=1.4, dash=None):
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        return (
            f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="{color}" stroke-width="{width}"{dash_attr}/>'
        )

    def text(x, y, body, size=16, anchor="start", color="#1f2933", weight="400"):
        return (
            f'<text x="{x:.1f}" y="{y:.1f}" font-family="Arial, sans-serif" '
            f'font-size="{size}" text-anchor="{anchor}" fill="{color}" '
            f'font-weight="{weight}">{escape(str(body), quote=False)}</text>'
        )

    floor = max(tol * 1e-9, 1e-18)
    log_min = math.log10(floor / 2)
    log_max = math.log10(tol * 1.8)

    def sy_err(value):
        value = max(value, floor)
        frac = (math.log10(value) - log_min) / (log_max - log_min)
        return left[1] + left[3] - frac * left[3]

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        text(315, 38, "Measured conservation error vs reference", 18, "middle", weight="700"),
        text(850, 38, "Peri fracture and DEM contact handoff", 18, "middle", weight="700"),
    ]

    for plot in (left, right):
        x, y, w, h = plot
        svg.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="none" stroke="#cbd5e0"/>')
        for i in range(1, 5):
            gy = y + h * i / 5
            svg.append(line(x, gy, x + w, gy, "#e2e8f0", 1))

    bar_w = 90
    for label, value, color, cx in [
        ("mass", rel_mass, "#2b6cb0", left[0] + 165),
        ("momentum", rel_p, "#c2410c", left[0] + 325),
    ]:
        y = sy_err(value)
        svg.append(f'<rect x="{cx - bar_w/2:.1f}" y="{y:.1f}" width="{bar_w}" height="{left[1] + left[3] - y:.1f}" fill="{color}"/>')
        svg.append(text(cx, left[1] + left[3] + 28, label, 15, "middle"))
        shown = "0.0" if value == 0.0 else f"{value:.1e}"
        svg.append(text(cx, y - 12, shown, 14, "middle"))

    svg.append(text(22, 285, "max relative conservation error", 14, "middle"))
    svg.append(text(left[0] - 10, sy_err(floor), f"{floor:.0e}", 12, "end"))

    svg.append(f'<polyline fill="none" stroke="#2f855a" stroke-width="3" points="{poly((sx(x), sy_count(y)) for x, y in zip(xs, bonds))}"/>')
    svg.append(f'<polyline fill="none" stroke="#805ad5" stroke-width="3" points="{poly((sx(x), sy_contact(y)) for x, y in zip(xs, contacts))}"/>')
    svg.append(f'<polyline fill="none" stroke="#718096" stroke-width="2.4" points="{poly((sx(x), sy_damage(y)) for x, y in zip(xs, damage))}"/>')

    svg.append(text(right[0] + right[2] / 2, right[1] + right[3] + 34, "simulation step", 15, "middle"))
    svg.append(text(right[0] - 46, right[1] + 180, "surviving peri bonds", 14, "middle", "#2f855a"))
    svg.append(text(right[0] + right[2] + 56, right[1] + 170, "DEM contacts / damage", 14, "middle", "#805ad5"))

    svg.append(text(right[0], right[1] + right[3] + 18, f"{min(xs)}", 12, "middle"))
    svg.append(text(right[0] + right[2], right[1] + right[3] + 18, f"{max(xs)}", 12, "middle"))

    svg.append("</svg>\n")

    with open(PLOT, "w") as fh:
        fh.write("\n".join(svg))
    print(f"wrote {PLOT}")


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
    rel = {m[1]: float(m[2]) for m in REL_RE.finditer(text)}
    tol_match = TOL_RE.search(text)
    if "mass" not in rel or "momentum" not in rel or not tol_match:
        print("CHECKS FAILED: conservation summary not found in log")
        return 1
    tol = float(tol_match[1])

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

    write_plot(steps, rel["mass"], rel["momentum"], tol)

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
