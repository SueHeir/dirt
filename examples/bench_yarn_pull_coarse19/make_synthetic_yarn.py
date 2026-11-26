"""Write a synthetic coarse-19-style yarn in LAMMPS bpm/sphere format.

The real coarse-19 data file comes from the private fibre generator and is not
committed. This makes a stand-in with the same layout so the yarn example can
be run, timed and MPI-checked anywhere:

* 19 fibres on a hexagonal 1 + 6 + 12 cross-section (or 1 + 6 + 12 + 18 = 37
  with ``--rings 3``), each a straight chain of overlapping beads twisted
  helically about the yarn axis;
* bead diameter 74, spacing 36 along z (the coarse-19 file's ~0.97 r);
* molecule id = fibre index, so intra-fibre contact is excluded;
* bond type 2 for bonds with both beads inside a grip (z < 500 or z > 9500),
  type 1 elsewhere, as in the coarse-19 file.

Usage:  python3 make_synthetic_yarn.py OUT.data [--length 10000] [--rings 2]
                                                [--plies 1]
"""
import argparse
import math

p = argparse.ArgumentParser()
p.add_argument("out")
p.add_argument("--length", type=float, default=10000.0)
p.add_argument("--rings", type=int, default=2)
p.add_argument("--plies", type=int, default=1, help="copies side by side in x/y")
p.add_argument("--twist", type=float, default=2.0 * math.pi / 20000.0, help="rad per unit z")
args = p.parse_args()

D = 74.0
R = D / 2.0
SPACING = 36.0
GRIP_LO, GRIP_HI = 500.0, args.length - 500.0

# Hexagonal fibre centres, ring by ring.
centres = [(0.0, 0.0)]
for ring in range(1, args.rings + 1):
    for side in range(6):
        a0 = math.pi / 3.0 * side
        a1 = math.pi / 3.0 * (side + 1)
        x0, y0 = ring * D * math.cos(a0), ring * D * math.sin(a0)
        x1, y1 = ring * D * math.cos(a1), ring * D * math.sin(a1)
        for k in range(ring):
            t = k / ring
            centres.append((x0 + (x1 - x0) * t, y0 + (y1 - y0) * t))

yarn_r = (args.rings + 0.5) * D
ply_offsets = [(0.0, 0.0)]
if args.plies > 1:
    ply_offsets = []
    side = math.ceil(math.sqrt(args.plies))
    for k in range(args.plies):
        ply_offsets.append(((k % side) * 2.2 * yarn_r, (k // side) * 2.2 * yarn_r))

nbead = int(args.length / SPACING) + 1
atoms, bonds = [], []
tag = 0
mol = 0
for ox, oy in ply_offsets:
    for cx, cy in centres:
        mol += 1
        prev = None
        for k in range(nbead):
            z = k * SPACING
            th = args.twist * z
            x = ox + cx * math.cos(th) - cy * math.sin(th)
            y = oy + cx * math.sin(th) + cy * math.cos(th)
            tag += 1
            atoms.append((tag, mol, 1, D, 1.0, x, y, z))
            if prev is not None:
                zp = (k - 1) * SPACING
                btype = 2 if (z < GRIP_LO and zp < GRIP_LO) or (z > GRIP_HI and zp > GRIP_HI) else 1
                bonds.append((len(bonds) + 1, btype, prev, tag))
            prev = tag

xs = [a[5] for a in atoms]
ys = [a[6] for a in atoms]
pad = 2.0 * D
with open(args.out, "w") as f:
    f.write("synthetic coarse yarn (make_synthetic_yarn.py)\n\n")
    f.write(f"{len(atoms)} atoms\n{len(bonds)} bonds\n1 atom types\n2 bond types\n\n")
    f.write(f"{min(xs) - pad:.6f} {max(xs) + pad:.6f} xlo xhi\n")
    f.write(f"{min(ys) - pad:.6f} {max(ys) + pad:.6f} ylo yhi\n")
    f.write(f"{-1500.0:.6f} {args.length + 1500.0:.6f} zlo zhi\n\n")
    f.write("Atoms # bpm/sphere\n\n")
    for a in atoms:
        f.write("%d %d %d %.6f %.6f %.6f %.6f %.6f\n" % a)
    f.write("\nBonds\n\n")
    for b in bonds:
        f.write("%d %d %d %d\n" % b)
print(f"{args.out}: {len(atoms)} atoms, {len(bonds)} bonds, {mol} fibres")
