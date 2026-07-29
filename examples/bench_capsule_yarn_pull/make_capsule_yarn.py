"""Convert the coarse-19 sphere yarn into a chain-of-capsules yarn.

Each fibre (molecule) is a chain of ~283 overlapping spheres. We group every k
consecutive beads into one sphero-cylinder: cylinder axis along the run, radius =
bead radius, half-length = half the first->last bead distance in the run. Capsules
within a fibre are bonded COM-to-COM (Guo BPM), reproducing the coarse-grain beam.

Emits a JSON spec the Rust example reads:
  { "radius": r, "capsules":[{com,quat,half_length,molecule,grip}], "bonds":[[i,j,type]] }
quat rotates +z (capsule body axis) to the world segment direction.

Usage: python3 make_capsule_yarn.py [k]   (default k=4)  -> data/capsule_yarn_k{k}.json
"""
import sys, json, math, os

K = int(sys.argv[1]) if len(sys.argv) > 1 else 4
SRC = os.path.join(os.path.dirname(__file__), "..", "bench_yarn_pull_coarse19",
                   "data", "yarn_gl10000_1ply.data")
LEFT_GRIP_ZMAX, RIGHT_GRIP_ZMIN = 500.0, 9500.0

# ── parse Atoms: tag -> (mol, x, y, z), radius ────────────────────────────────
atoms, radius = {}, None
sec = None
for line in open(SRC):
    s = line.strip()
    head = s.split("#")[0].strip()  # section headers may carry a "# style" comment
    if head in ("Atoms", "Velocities", "Bonds", "Masses"):
        sec = head; continue
    if sec == "Atoms" and s and not s.startswith("#"):
        f = s.split()
        if len(f) >= 8:
            tag, mol = int(f[0]), int(f[1])
            diam, x, y, z = float(f[3]), float(f[5]), float(f[6]), float(f[7])
            atoms[tag] = (mol, x, y, z)
            radius = diam / 2.0
    elif sec == "Bonds":
        break

# group tags by fibre, preserving file order (sequential along the fibre)
fibres = {}
for tag in sorted(atoms):
    mol = atoms[tag][0]
    fibres.setdefault(mol, []).append(tag)

def quat_from_axis(ax):
    """Unit quaternion (w,x,y,z) rotating +z to unit vector ax."""
    z = (0.0, 0.0, 1.0)
    d = ax[2]                                   # z . ax
    if d > 0.999999:
        return (1.0, 0.0, 0.0, 0.0)
    if d < -0.999999:
        return (0.0, 1.0, 0.0, 0.0)             # 180° about x
    c = (z[1]*ax[2]-z[2]*ax[1], z[2]*ax[0]-z[0]*ax[2], z[0]*ax[1]-z[1]*ax[0])  # z × ax
    w = 1.0 + d
    n = math.sqrt(w*w + c[0]*c[0]+c[1]*c[1]+c[2]*c[2])
    return (w/n, c[0]/n, c[1]/n, c[2]/n)

capsules, bonds = [], []
aspect = []
for mol, tags in fibres.items():
    first_cap = len(capsules)
    for c0 in range(0, len(tags), K):
        run = tags[c0:c0 + K]
        if len(run) < 2:
            # trailing singleton: attach as a tiny capsule at the bead
            _, x, y, z = atoms[run[0]]
            capsules.append(dict(com=[x, y, z], quat=[1, 0, 0, 0], half_length=radius * 0.5,
                                 molecule=mol))
            continue
        p0 = atoms[run[0]][1:]
        p1 = atoms[run[-1]][1:]
        com = [(p0[i] + p1[i]) / 2 for i in range(3)]
        d = [p1[i] - p0[i] for i in range(3)]
        L = math.sqrt(sum(di*di for di in d))
        ax = [di / L for di in d]
        half = L / 2.0
        capsules.append(dict(com=com, quat=list(quat_from_axis(ax)),
                             half_length=half, molecule=mol))
        aspect.append((half + radius) / radius)
    # bond consecutive capsules in this fibre (type 1 = breakable interior)
    for i in range(first_cap, len(capsules) - 1):
        bonds.append([i, i + 1, 1])

# grip flag: capsule COM in a z-end region (rigid clamp, like the sphere yarn)
for c in capsules:
    z = c["com"][2]
    c["grip"] = 0 if LEFT_GRIP_ZMAX <= z <= RIGHT_GRIP_ZMIN else (1 if z < LEFT_GRIP_ZMAX else 2)

# ── emit a plain-text spec (no serde_json dep in the binary crate) ────────────
#   line 1:            radius  k  n_capsules  n_bonds
#   n_capsules lines:  cx cy cz  qw qx qy qz  half_length  molecule  grip
#   n_bonds lines:     i j type
outdir = os.path.join(os.path.dirname(__file__), "data")
os.makedirs(outdir, exist_ok=True)
path = os.path.join(outdir, f"capsule_yarn_k{K}.txt")
with open(path, "w") as fh:
    fh.write(f"{radius:.10g} {K} {len(capsules)} {len(bonds)}\n")
    for c in capsules:
        co, q = c["com"], c["quat"]
        fh.write("%.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %d %d\n" % (
            co[0], co[1], co[2], q[0], q[1], q[2], q[3], c["half_length"],
            c["molecule"], c["grip"]))
    for b in bonds:
        fh.write("%d %d %d\n" % (b[0], b[1], b[2]))
ng = sum(1 for c in capsules if c["grip"])
print(f"k={K}: {len(atoms)} spheres -> {len(capsules)} capsules ({len(bonds)} bonds), "
      f"{len(fibres)} fibres, radius={radius:.2f}")
print(f"  mean aspect ratio (tip-to-tip / diameter): {sum(aspect)/len(aspect):.2f}")
print(f"  grip capsules: {ng}   interior: {len(capsules)-ng}")
print(f"  wrote {path}")
