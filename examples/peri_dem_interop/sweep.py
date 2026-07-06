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

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "peri_dem_interop"
CONFIG = os.path.join("examples", EXAMPLE, "config.toml")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
LOG = os.path.join(DATA_DIR, "run.log")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
PLOT = os.path.join(PLOT_DIR, "peri_dem_transition_validation.png")

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
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    os.makedirs(PLOT_DIR, exist_ok=True)
    plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150, "font.size": 10})

    xs = [s[0] for s in steps]
    bonds = [s[2] for s in steps]
    damage = [s[3] for s in steps]
    contacts = [s[4] for s in steps]
    bonds0 = bonds[0]
    bond_gate = 0.10 * bonds0
    contact_gate = 8

    fig, (ax_err, ax_counts) = plt.subplots(1, 2, figsize=(11.5, 4.6))

    labels = ["mass", "momentum"]
    vals = [rel_mass, rel_p]
    floor = max(tol * 1e-9, 1e-18)
    plot_vals = [max(v, floor) for v in vals]
    bars = ax_err.bar(labels, plot_vals, color=["#2b6cb0", "#c2410c"], width=0.55)
    for bar, val in zip(bars, vals):
        label = "0.0" if val == 0.0 else f"{val:.1e}"
        ax_err.text(bar.get_x() + bar.get_width() / 2, bar.get_height() * 1.25,
                    label, ha="center", va="bottom", fontsize=9)
    ax_err.axhline(tol, color="black", linestyle="--", linewidth=1.5,
                   label=f"PASS limit = {tol:.0e}")
    ax_err.set_yscale("log")
    ax_err.set_ylim(floor / 2, tol * 1.8)
    ax_err.set_ylabel("max relative conservation error")
    ax_err.set_title("Measured conservation error vs reference")
    ax_err.grid(True, which="both", axis="y", alpha=0.25)
    ax_err.legend(loc="upper left", frameon=False)

    ax_counts.plot(xs, bonds, color="#2f855a", linewidth=2.0, label="surviving peri bonds")
    ax_counts.axhline(bond_gate, color="#2f855a", linestyle="--", linewidth=1.2,
                      label="fracture PASS: <10% initial bonds")
    ax_counts.set_xlabel("simulation step")
    ax_counts.set_ylabel("surviving peri bonds")
    ax_counts.grid(True, axis="both", alpha=0.25)
    ax_counts.set_title("Peri fracture and DEM contact handoff")

    ax_contact = ax_counts.twinx()
    ax_contact.plot(xs, contacts, color="#805ad5", linewidth=2.0, label="active DEM contacts")
    ax_contact.axhline(contact_gate, color="#805ad5", linestyle=":", linewidth=1.4,
                       label="contact PASS: >=8")
    ax_contact.set_ylabel("active DEM contacts")

    ax_damage = ax_counts.twinx()
    ax_damage.spines["right"].set_position(("axes", 1.15))
    ax_damage.plot(xs, damage, color="#718096", linewidth=1.6, alpha=0.9,
                   label="peak damage")
    ax_damage.axhline(0.99, color="#718096", linestyle="-.", linewidth=1.1,
                      label="damage PASS: >=0.99")
    ax_damage.set_ylim(0, 1.05)
    ax_damage.set_ylabel("peak damage")

    handles = [
        Line2D([0], [0], color="#2f855a", lw=2, label="surviving peri bonds"),
        Line2D([0], [0], color="#2f855a", lw=1.2, ls="--",
               label="fracture PASS: <10% initial bonds"),
        Line2D([0], [0], color="#805ad5", lw=2, label="active DEM contacts"),
        Line2D([0], [0], color="#805ad5", lw=1.4, ls=":",
               label="contact PASS: >=8"),
        Line2D([0], [0], color="#718096", lw=1.6, label="peak damage"),
        Line2D([0], [0], color="#718096", lw=1.1, ls="-.",
               label="damage PASS: >=0.99"),
    ]
    ax_counts.legend(handles=handles, loc="center left", bbox_to_anchor=(0.02, 0.55),
                     frameon=False, fontsize=8)

    fig.tight_layout()
    fig.savefig(PLOT, bbox_inches="tight")
    plt.close(fig)
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
