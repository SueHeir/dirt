"""Build a 4-ply coarse-19 yarn from the relaxed 1-ply, porting the ply assembly
of ../fiber_config/fibers_all.py (BondBPMGuo, ~lines 2303-2386).

The 4-ply is four copies of the 19-fibre ply placed at the corners (+/-R, +/-R)
with R = yarn_radius * ply_percent_sep (181.8 um), each keeping its internal
fibre twist. The whole set of plies then spirals about the yarn axis at the ply
twist ply_rotation = 0.00056 rad/um (switch_rotation=True -> negative sense).

Shortcut vs regrowing helices: reuse the already-relaxed 1-ply (fibre twist baked
in), replicate x4 to the corners, and apply only the ply-level rotation
theta(z) = -ply_rotation * z about the yarn axis. Within-fibre bonds are copied x4
with renumbered tags; each of the 76 fibres gets a unique molecule id.

Output: data/yarn_gl10000_4ply.data (same LAMMPS bpm/sphere format as the 1-ply).
"""
import math, os

SRC = os.path.join(os.path.dirname(__file__), "data", "yarn_gl10000_1ply.data")
OUT = os.path.join(os.path.dirname(__file__), "data", "yarn_gl10000_4ply.data")

YARN_RADIUS = 180.0
PLY_PERCENT_SEP = 1.01
R = YARN_RADIUS * PLY_PERCENT_SEP                      # 181.8 um corner offset
PLY_ROTATION = 2.0 * math.pi / (2594 + 2976 + 2824 * 2.0)  # 0.000560 rad/um (4-ply)
# ply corner offsets (left, right, up, down) exactly as fibers_all.py
PLY_OFFSETS = [(+R, +R), (-R, +R), (+R, -R), (-R, -R)]

# ── parse the 1-ply data ──────────────────────────────────────────────────────
atoms, bonds = [], []
sec = None
box = {}
for line in open(SRC):
    s = line.rstrip("\n")
    t = s.strip()
    if t.endswith(("xlo xhi", "ylo yhi", "zlo zhi")):
        f = t.split(); box[f[-2]] = (float(f[0]), float(f[1]))
    if t.startswith("Atoms"): sec = "A"; continue
    if t.startswith("Bonds"): sec = "B"; continue
    if t.startswith("Velocities"): sec = "V"; continue
    if not t: continue
    if sec == "A":
        f = t.split()
        # id mol type diam dens x y z [img img img]
        atoms.append([int(f[0]), int(f[1]), int(f[2]), f[3], f[4],
                      float(f[5]), float(f[6]), float(f[7]), f[8:]])
    elif sec == "B":
        f = t.split()
        bonds.append([int(f[0]), int(f[1]), int(f[2]), int(f[3])])
    elif sec == "V":
        sec = None  # ignore velocities

N = len(atoms); NB = len(bonds); NFIB = max(a[1] for a in atoms)  # 19
print(f"1-ply: {N} atoms, {NB} bonds, {NFIB} fibres, ply offset R={R:.1f}, ply_rot={PLY_ROTATION:.6f}")

# ── replicate into 4 plies with corner offset + ply-level helical twist ───────
new_atoms, new_bonds = [], []
for p, (dx, dy) in enumerate(PLY_OFFSETS):
    for a in atoms:
        aid, mol, typ, diam, dens, x, y, z, img = a
        # offset to this ply's corner, then rotate about the yarn axis by theta(z)
        px, py = x + dx, y + dy
        th = -PLY_ROTATION * z                      # switch_rotation=True -> negative
        c, s = math.cos(th), math.sin(th)
        nx = px * c - py * s
        ny = px * s + py * c
        new_atoms.append([p * N + aid, p * NFIB + mol, typ, diam, dens, nx, ny, z, img])
    for b in bonds:
        bid, bt, a1, a2 = b
        new_bonds.append([p * NB + bid, bt, p * N + a1, p * N + a2])

# ── box: 4-ply spans ~ (R*sqrt2 + yarn_radius) ~ 437 um radially ──────────────
rmax = max(math.hypot(a[5], a[6]) for a in new_atoms)
half = math.ceil((rmax + 40) / 50) * 50
print(f"4-ply: {len(new_atoms)} atoms, {len(new_bonds)} bonds, max radius={rmax:.0f} -> box +/-{half}")

with open(OUT, "w") as f:
    f.write("LAMMPS data file (4-ply coarse19 from make_4ply_yarn.py), units = lj\n\n")
    f.write(f"{len(new_atoms)} atoms\n1 atom types\n{len(new_bonds)} bonds\n2 bond types\n\n")
    f.write(f"{-half} {half} xlo xhi\n{-half} {half} ylo yhi\n")
    f.write(f"{box['zlo'][0]:.0f} {box['zlo'][1]:.0f} zlo zhi\n\n")
    f.write("Atoms # bpm/sphere\n\n")
    for a in new_atoms:
        img = " ".join(a[8]) if a[8] else "0 0 0"
        f.write(f"{a[0]} {a[1]} {a[2]} {a[3]} {a[4]} {a[5]:.6g} {a[6]:.6g} {a[7]:.6g} {img}\n")
    f.write("\nBonds\n\n")
    for b in new_bonds:
        f.write(f"{b[0]} {b[1]} {b[2]} {b[3]}\n")
print(f"wrote {OUT}")
