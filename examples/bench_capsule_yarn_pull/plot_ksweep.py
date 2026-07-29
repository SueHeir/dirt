"""k-sweep of the sphero-cylinder yarn (k = 2,4,8 beads/capsule) vs the sphere yarn.

Left: stress-strain (results). Right: run time + particle count (performance).
Reads the capsule runs from the scratchpad sweep dir (edit SWEEP if elsewhere) and
the committed sphere-yarn curve. Full 560k-step native serial runs.
"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SWEEP = ("/private/tmp/claude-501/-Users-suehr-Documents-GitHub-dirt/"
         "84d33cd1-a33e-4438-972e-653add533b92/scratchpad/ksweepR")
SPHERE = "examples/bench_yarn_pull_coarse19/data/stress_strain.csv"

def curve(path, scol=2):
    xs, ys = [], []
    for l in open(path).read().splitlines()[1:]:
        f = l.split(","); xs.append(float(f[1]) * 100); ys.append(-float(f[scol]))
    return xs, ys

# (label, path, colour, time_s, particles)
RUNS = [
    ("sphere (5330)", SPHERE, "#000000", 392.0, 5330),
    ("k=2 (2674)", f"{SWEEP}/k2/data/stress_strain.csv", "#1f77b4", 434.1, 2674),
    ("k=4 (1340)", f"{SWEEP}/k4/data/stress_strain.csv", "#2ca02c", 303.6, 1340),
    ("k=8 (676)",  f"{SWEEP}/k8/data/stress_strain.csv", "#d62728", 183.3, 676),
]

fig, ax = plt.subplots(1, 2, figsize=(13, 4.9))

# ── Panel 1: stress-strain ────────────────────────────────────────────────────
for label, path, col, _, _ in RUNS:
    x, y = curve(path)
    ax[0].plot(x, y, "-" if "sphere" in label else "--", color=col, lw=2, label=label)
ax[0].set_xlabel("axial strain (%)")
ax[0].set_ylabel(r"$\sigma_{zz} = -\mathrm{virial}/V_{box}$ (lj)")
ax[0].set_title("Results — stress-strain (sphere vs sphero-cylinder k-sweep)")
ax[0].legend(fontsize=9); ax[0].grid(alpha=0.3); ax[0].set_xlim(0, 1.35)
ax[0].text(0.03, 0.070, "relaxed: elastic modulus matches the sphere yarn;\n"
                        "capsules stronger + fail later (rigid-segment bending)",
           fontsize=7.5, color="#555555")

# ── Panel 2: performance (time bars + particle count) ─────────────────────────
labels = [r.split(" ")[0] for r in [x[0] for x in RUNS]]
times = [r[3] for r in RUNS]
parts = [r[4] for r in RUNS]
cols = [r[2] for r in RUNS]
xpos = np.arange(len(RUNS))
bars = ax[1].bar(xpos, times, color=cols, width=0.6)
for b, t, p in zip(bars, times, parts):
    ax[1].text(b.get_x() + b.get_width() / 2, t + 5, f"{t:.0f}s", ha="center", fontsize=9)
ax[1].set_xticks(xpos); ax[1].set_xticklabels(labels)
ax[1].set_ylabel("wall-clock (s), 660k steps (100k relax + pull)")
ax[1].set_title("Performance — run time; particle count on right axis")
ax[1].grid(alpha=0.3, axis="y"); ax[1].set_ylim(0, 440)
# particle count overlay
ax2 = ax[1].twinx()
ax2.plot(xpos, parts, "o--", color="#888888", lw=1.5, ms=7, label="particles")
for xx, p in zip(xpos, parts):
    ax2.text(xx, p + 150, f"{p}", ha="center", fontsize=8, color="#555555")
ax2.set_ylabel("particle count", color="#555555")
ax2.set_ylim(0, 6200); ax2.tick_params(axis="y", colors="#555555")

fig.suptitle("Sphero-cylinder yarn k-sweep vs sphere yarn — coarse-19 1-ply pull to failure (free-grip damped relax)", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])
out = "examples/bench_capsule_yarn_pull/plots/ksweep.png"
import os; os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=140)
for label, path, _, t, p in RUNS:
    x, y = curve(path); pk = max(y)
    print(f"{label:16s} peak {pk:.4f} @ {x[y.index(pk)]:.3f}%   {t:.0f}s   {p} particles")
print("wrote", out)
