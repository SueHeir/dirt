"""Merge the per-rank MPI dumps (traj_r0..rN.lammpstrj) into one OVITO-ready
trajectory. Each rank writes only its OWNED atoms per frame; atoms migrate
between ranks over time, so we group by timestep and concatenate all ranks'
atoms (deduped by id) into a single frame.

Usage: python3 merge_traj.py <dir_with_data/>   (default: calib/4ply_mpi)
Output: <dir>/data/traj.lammpstrj
"""
import sys, glob, os, collections

D = sys.argv[1] if len(sys.argv) > 1 else "examples/bench_yarn_pull_coarse19/calib/4ply_mpi"
files = sorted(glob.glob(os.path.join(D, "data", "traj_r*.lammpstrj")))
if not files:
    sys.exit(f"no per-rank traj files in {D}/data")

# timestep -> {atom_id: line}
frames = collections.OrderedDict()
box = None
header = "ITEM: ATOMS id mol type x y z radius"
for fn in files:
    L = open(fn).read().splitlines()
    i = 0
    while i < len(L):
        if L[i].startswith("ITEM: TIMESTEP"):
            step = int(L[i + 1]); i += 2; continue
        if L[i].startswith("ITEM: BOX"):
            box = L[i:i + 4]; i += 4; continue
        if L[i].startswith("ITEM: ATOMS"):
            i += 1
            d = frames.setdefault(step, {})
            while i < len(L) and not L[i].startswith("ITEM"):
                aid = int(L[i].split(maxsplit=1)[0])
                d[aid] = L[i]; i += 1
            continue
        i += 1

out = os.path.join(D, "data", "traj.lammpstrj")
nframes = 0
with open(out, "w") as f:
    for step, d in frames.items():
        f.write(f"ITEM: TIMESTEP\n{step}\n")
        f.write(f"ITEM: NUMBER OF ATOMS\n{len(d)}\n")
        f.write("\n".join(box) + "\n" if box else
                "ITEM: BOX BOUNDS ff ff ff\n-500 500\n-500 500\n-1500 11500\n")
        f.write(header + "\n")
        for aid in sorted(d):
            f.write(d[aid] + "\n")
        nframes += 1
counts = {s: len(d) for s, d in frames.items()}
print(f"merged {len(files)} ranks -> {out}  ({nframes} frames, "
      f"{min(counts.values())}-{max(counts.values())} atoms/frame)")
