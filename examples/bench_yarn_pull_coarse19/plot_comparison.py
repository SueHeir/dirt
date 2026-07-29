"""Overlay the DIRT coarse-19 yarn pull against the LAMMPS bpm/guo reference.

Reads only committed data under examples/bench_yarn_pull_coarse19/data/:
  - stress_strain.csv          DIRT run output (step,strain,stress_zz,...,bond_count)
  - breaks.log                 DIRT true bond-break events (YARN_BREAK_LOG=1)
  - lammps_stress_strain.dat   LAMMPS reference (strain P sxx syy szz ... nbond)

Run from the repo root:  python3 examples/bench_yarn_pull_coarse19/plot_comparison.py
"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

EX = "examples/bench_yarn_pull_coarse19"
N0 = 5311                                   # initial bonds
STEP2STRAIN = 0.000216 / 9000 * 100.0       # step -> strain %
VBOX = 720 * 720 * 13000                    # LAMMPS box volume (stress normalisation)

# Wall-clock (s), native builds, 560k steps, DIRT best case (neighbour every=100 +
# molecule/intra exclusion). Full runs; DIRT via /usr/bin/time, LAMMPS loop time.
#            serial  2-rank
TIMING_DIRT   = [392.0, 190.0]
TIMING_LAMMPS = [368.0, 198.0]

# ── DIRT stress-strain + bond_count metric ────────────────────────────────────
dx, ds, dmetric = [], [], []
for l in open(f"{EX}/data/stress_strain.csv").read().splitlines()[1:]:
    f = l.split(",")
    dx.append(float(f[1]) * 100); ds.append(-float(f[2])); dmetric.append(int(f[4]))

# ── DIRT true bond survival from the break log ────────────────────────────────
breaks = sorted(int(l.split("step=")[1].split()[0]) * STEP2STRAIN
                for l in open(f"{EX}/data/breaks.log") if "[break]" in l)
bt_x, bt_y = [0.0], [N0]
for s in breaks:
    bt_x += [s, s]; bt_y += [bt_y[-1], bt_y[-1] - 1]
bt_x.append(max(dx)); bt_y.append(bt_y[-1])

# ── LAMMPS reference ──────────────────────────────────────────────────────────
lx, ls, lb = [], [], []
for l in open(f"{EX}/data/lammps_stress_strain.dat").read().splitlines():
    if l.startswith("#") or not l.strip():
        continue
    f = l.split()
    lx.append(float(f[0]) * 100); ls.append(float(f[4]) / VBOX); lb.append(int(float(f[8])) // 2)

fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.8))

# Panel 1 — stress-strain
ax[0].plot(lx, ls, "-", color="#1f77b4", lw=2.2, label="LAMMPS")
ax[0].plot(dx, ds, "--", color="#2ca02c", lw=2, label="DIRT")
ax[0].set_xlabel("axial strain (%)")
ax[0].set_ylabel(r"$\sigma_{zz} = -\mathrm{virial}_{zz}/V_{box}$ (lj)")
ax[0].set_title("Stress–strain")
ax[0].legend(); ax[0].grid(alpha=0.3); ax[0].set_xlim(0, 1.35)
ax[0].annotate("peak  LAMMPS %.4f\n          DIRT %.4f" % (max(ls), max(ds)),
               (0.60, 0.050), fontsize=8)

# Panel 2 — bond survival (true breaks; metric shown dotted for context)
ax[1].plot(lx, [100 * b / N0 for b in lb], "-", color="#1f77b4", lw=2.2,
           label=f"LAMMPS true ({N0 - lb[-1]} broke)")
ax[1].plot(bt_x, [100 * b / N0 for b in bt_y], "-", color="#2ca02c", lw=2.4,
           label=f"DIRT true ({len(breaks)} broke)")
ax[1].plot(dx, [100 * b / N0 for b in dmetric], ":", color="#888888", lw=1.8,
           label="DIRT bond_count metric*")
ax[1].set_xlabel("axial strain (%)"); ax[1].set_ylabel("bonds intact (%)")
ax[1].set_title("Bond survival"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
ax[1].set_xlim(0, 1.35); ax[1].set_ylim(74, 100.6)
ax[1].text(0.02, 75, "*metric drop after fracture = intact bonds whose partner\n"
                     "scattered past the ghost cutoff, not breakage",
           fontsize=7, color="#888888")

# Panel 3 — wall-clock: DIRT vs LAMMPS at serial and 2-rank (grouped)
import numpy as np
x = np.arange(2); w = 0.38
bd = ax[2].bar(x - w / 2, TIMING_DIRT, w, color="#2ca02c", label="DIRT")
bl = ax[2].bar(x + w / 2, TIMING_LAMMPS, w, color="#1f77b4", label="LAMMPS")
for bars in (bd, bl):
    for b in bars:
        ax[2].text(b.get_x() + b.get_width() / 2, b.get_height() + 8,
                   f"{b.get_height():.0f}s", ha="center", fontsize=9)
for k in range(2):
    ax[2].text(k, 40, "DIRT/LAMMPS\n%.2f×" % (TIMING_DIRT[k] / TIMING_LAMMPS[k]),
               ha="center", fontsize=8, color="#333333")
ax[2].set_xticks(x); ax[2].set_xticklabels(["serial (1 rank)", "MPI (2 ranks)"])
ax[2].set_ylabel("wall-clock (s), 560k steps")
ax[2].set_title("Run time (native, best case)")
ax[2].legend(loc="upper right"); ax[2].grid(alpha=0.3, axis="y")
ax[2].set_ylim(0, 470)

fig.suptitle("DIRT vs LAMMPS — coarse-19 1-ply yarn pull "
             "(matched physics; DIRT at parity on run time, serial and 2-rank)", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.93])
out = f"{EX}/plots/stress_strain_comparison.png"
fig.savefig(out, dpi=140)
print("DIRT   peak %.4f at %.3f%%  (%d bonds broke, %.3f-%.3f%%)"
      % (max(ds), dx[ds.index(max(ds))], len(breaks), breaks[0], breaks[-1]))
print("LAMMPS peak %.4f  (%d bonds broke)" % (max(ls), N0 - lb[-1]))
print("wrote", out)
