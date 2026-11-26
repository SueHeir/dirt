"""Build a 4-ply sphero-cylinder yarn from the 1-ply capsule spec, mirroring
make_4ply_yarn.py (spheres) but rotating each capsule's ORIENTATION with the
ply-twist as well as its position.

Four copies of the 19-fibre capsule ply at corners (+/-R, +/-R), R = 181.8 um,
each spiralling about the yarn axis at ply_rotation = 0.00056 rad/um. For a
capsule at (cx,cy,cz) with quaternion q, ply-twist theta(z) = -ply_rot*cz:
  position: rotate (cx+dx, cy+dy) about z by theta
  orientation: q_new = q_z(theta) (x) q     (q_z = rotation about world z)

Usage: python3 make_capsule_4ply.py [k]   (default 8) -> data/capsule_yarn_4ply_k{k}.txt
"""
import sys, math, os

K = int(sys.argv[1]) if len(sys.argv) > 1 else 8
SRC = os.path.join(os.path.dirname(__file__), "data", f"capsule_yarn_k{K}.txt")
OUT = os.path.join(os.path.dirname(__file__), "data", f"capsule_yarn_4ply_k{K}.txt")

YARN_RADIUS = 180.0
R = YARN_RADIUS * 1.01                                   # 181.8 um corner offset
PLY_ROTATION = 2.0 * math.pi / (2594 + 2976 + 2824 * 2.0)  # 0.000560 rad/um
PLY_OFFSETS = [(+R, +R), (-R, +R), (+R, -R), (-R, -R)]

def axis_from_quat(q):        # rotate +z by q=[w,x,y,z]
    w, x, y, z = q
    return [2*(x*z + w*y), 2*(y*z - w*x), 1 - 2*(x*x + y*y)]

def quat_from_axis(ax):       # unit quat [w,x,y,z] rotating +z to unit ax
    d = ax[2]
    if d > 0.999999:
        return [1.0, 0.0, 0.0, 0.0]
    if d < -0.999999:
        return [0.0, 1.0, 0.0, 0.0]
    c = (-ax[1], ax[0], 0.0)  # z x ax = (0,0,1) x ax = (-ax1, ax0, 0)
    w = 1.0 + d
    n = math.sqrt(w*w + c[0]*c[0] + c[1]*c[1] + c[2]*c[2])
    return [w/n, c[0]/n, c[1]/n, c[2]/n]

def ply_transform(pt, dx, dy):  # offset to ply corner, then ply-twist about yarn z by theta(z)
    px, py = pt[0] + dx, pt[1] + dy
    th = -PLY_ROTATION * pt[2]
    c, s = math.cos(th), math.sin(th)
    return [px * c - py * s, px * s + py * c, pt[2]]

# ── parse the 1-ply capsule spec ──────────────────────────────────────────────
lines = open(SRC).read().splitlines()
radius, k, ncap, nbond = lines[0].split()
radius, ncap, nbond = float(radius), int(ncap), int(nbond)
caps = []   # (com[3], quat[4], half_length, mol, grip)
for i in range(ncap):
    f = lines[1 + i].split()
    caps.append(([float(f[0]), float(f[1]), float(f[2])],
                 [float(f[3]), float(f[4]), float(f[5]), float(f[6])],
                 float(f[7]), int(f[8]), int(f[9])))
bonds = []
for i in range(nbond):
    f = lines[1 + ncap + i].split()
    bonds.append((int(f[0]), int(f[1]), int(f[2])))
NFIB = max(c[3] for c in caps)   # 19
print(f"1-ply capsules: {ncap}, bonds {nbond}, fibres {NFIB}, radius {radius:.2f}")

# ── replicate x4 with corner offset + ply-level helical twist ─────────────────
new_caps, new_bonds = [], []
for p, (dx, dy) in enumerate(PLY_OFFSETS):
    for com, q, hl, mol, grip in caps:
        # transform the capsule's two ENDPOINTS through the ply-twist (each at its
        # own z), so the rebuilt axis follows the actual twisted fibre tangent.
        ax = axis_from_quat(q)
        e1 = ply_transform([com[j] - hl * ax[j] for j in range(3)], dx, dy)
        e2 = ply_transform([com[j] + hl * ax[j] for j in range(3)], dx, dy)
        ncom = [(e1[j] + e2[j]) / 2 for j in range(3)]
        d = [e2[j] - e1[j] for j in range(3)]
        nlen = math.sqrt(sum(v * v for v in d))
        nax = [v / nlen for v in d] if nlen else ax
        new_caps.append((ncom, quat_from_axis(nax), nlen / 2.0, p * NFIB + mol, grip))
    for i, j, bt in bonds:
        new_bonds.append((p * ncap + i, p * ncap + j, bt))

rmax = max(math.hypot(c[0][0], c[0][1]) for c in new_caps)
print(f"4-ply capsules: {len(new_caps)}, bonds {len(new_bonds)}, max radius {rmax:.0f}")

with open(OUT, "w") as f:
    f.write(f"{radius:.10g} {K} {len(new_caps)} {len(new_bonds)}\n")
    for com, q, hl, mol, grip in new_caps:
        f.write("%.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %d %d\n" % (
            com[0], com[1], com[2], q[0], q[1], q[2], q[3], hl, mol, grip))
    for i, j, bt in new_bonds:
        f.write(f"{i} {j} {bt}\n")
print(f"wrote {OUT}")
