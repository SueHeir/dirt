#!/usr/bin/env python3
"""Static bond / contact-exclusion audit of the 44-sub-sphere monkey.

Investigation artifact for goal `dirt-bpm-monkey-blowup-investigation`. Answers,
without running the engine: does the AS-BUILT bpm monkey leave any overlapping (or
near-overlapping) intra-monkey sub-sphere pair CONTACT-ACTIVE — i.e. not excluded
by the 1-2 / 1-3 bond rule — so that the Hertz law would see the ~23.6 % build
overlap directly?

It mirrors the engine exactly:
  * auto_bond  — `dirt_bond::auto_bond_touching`: bond (i,j) if
                 dist <= bond_tolerance * (R_i + R_j); r0 = dist.
  * exclusion  — `soil_core::bond::BondStore::are_excluded`: exclude a pair from
                 contact iff it is directly bonded (1-2) OR shares a bonded
                 neighbour (1-3). NOTHING beyond 1-3 is excluded — unlike a rigid
                 clump, whose sub-spheres skip ALL same-body contact.

Usage:
  python3 bond_exclusion_audit.py [monkey_Deq0.1.toml] [--tol 1.1]
"""
import argparse
import itertools
import math
import re
import sys


def parse_spheres(path):
    txt = open(path).read()
    sph = []
    for m in re.finditer(
        r"offset\s*=\s*\[([^\]]+)\].*?radius\s*=\s*([0-9.eE+-]+)", txt
    ):
        x, y, z = (float(v) for v in m.group(1).split(","))
        sph.append((x, y, z, float(m.group(2))))
    return sph


def dist(a, b):
    return math.sqrt(sum((a[k] - b[k]) ** 2 for k in range(3)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("toml", nargs="?", default="examples/monkey_shear/monkey_Deq0.1.toml")
    ap.add_argument("--tol", type=float, default=1.1, help="bond_tolerance")
    a = ap.parse_args()

    sph = parse_spheres(a.toml)
    n = len(sph)
    print(f"parsed {n} sub-spheres from {a.toml}")

    # auto_bond graph
    bonds = {i: set() for i in range(n)}
    nb = 0
    for i, j in itertools.combinations(range(n), 2):
        if dist(sph[i], sph[j]) <= a.tol * (sph[i][3] + sph[j][3]):
            bonds[i].add(j)
            bonds[j].add(i)
            nb += 1
    print(f"auto_bond pairs (bonds_per_monkey) @ tol {a.tol} = {nb}")

    def excluded(i, j):
        if j in bonds[i]:
            return "1-2"
        if bonds[i] & bonds[j]:
            return "1-3"
        return None

    active_overlaps = []  # non-excluded AND geometrically overlapping (t=0 injection)
    n_excluded = 0
    n_overlap_excluded = 0
    active_gaps = []  # (dist-(Ri+Rj))/(Ri+Rj) for non-excluded, non-overlapping pairs
    for i, j in itertools.combinations(range(n), 2):
        d = dist(sph[i], sph[j])
        sr = sph[i][3] + sph[j][3]
        exc = excluded(i, j)
        overlapping = d < sr
        if exc:
            n_excluded += 1
            if overlapping:
                n_overlap_excluded += 1
        elif overlapping:
            active_overlaps.append((i, j, (sr - d) / sr))
        else:
            active_gaps.append((d - sr) / sr)

    print("\n--- exclusion (1-2 / 1-3 only) ---")
    print(f"excluded pairs                         : {n_excluded} / {n*(n-1)//2}")
    print(f"overlapping pairs that ARE excluded    : {n_overlap_excluded}")
    print(
        f"overlapping pairs CONTACT-ACTIVE (t=0) : {len(active_overlaps)}"
        "   <-- >0 would be a direct build-overlap Hertz injection"
    )
    for i, j, r in sorted(active_overlaps, key=lambda t: -t[2])[:10]:
        print(f"    active build-overlap {i:2d}-{j:2d}  {r*100:.1f}% of (Ri+Rj)")

    active_gaps.sort()
    print(
        f"\n--- 'priming' of the {len(active_gaps)} non-excluded, non-overlapping "
        "intra pairs ---"
    )
    print("gap=(dist-(Ri+Rj))/(Ri+Rj); a small gap closes to self-contact under small strain")
    for pct in (0, 1, 2, 5, 10):
        k = min(len(active_gaps) - 1, int(pct / 100 * len(active_gaps)))
        print(f"  {pct:2d}th pctile gap = {active_gaps[k]*100:6.2f}%")
    print(f"  nearest contact-active intra pair sits at a {active_gaps[0]*100:.1f}% gap")
    print(
        "\nInterpretation: 0 active build-overlaps ⇒ no t=0 Hertz injection; a large\n"
        "nearest-gap ⇒ intra-body self-contact needs substantial deformation, so it is\n"
        "a consequence of shear heating rather than its initiator."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
