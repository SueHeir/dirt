#!/usr/bin/env python3
"""Strong + weak scaling, DIRT-vs-LAMMPS, on controlled dense/dilute gas configs.

Self-contained companion to the perf_mpi_scaling example. Runs the *exact* dense
(phi=0.50, sc lattice) and dilute (phi=0.05, random) granular-gas configs from our
LAMMPS benchmark — not the example's gas/bed framework — in both scaling modes,
capturing a per-phase time breakdown from DIRT's scheduler report, and emits two
3-panel figures:

    plots/strong_scaling.png   # fixed total N, more cores
    plots/weak_scaling.png     # fixed particles/core, total grows with cores

Each figure: (1) DIRT/LAMMPS throughput ratio vs cores, (2) parallel efficiency
(DIRT & LAMMPS, both regimes), (3) phase panel — breakdown (weak) / per-phase
strong-scaling efficiency (strong).

    python3 examples/perf_mpi_scaling/scaling_gas.py            # build + run + graph, both modes
    python3 examples/perf_mpi_scaling/scaling_gas.py strong     # one mode
    python3 examples/perf_mpi_scaling/scaling_gas.py weak
    python3 examples/perf_mpi_scaling/scaling_gas.py graph      # replot from saved json, no runs

LAMMPS is optional but, for a fair comparison, should be built native and on PATH
(or PERF_LAMMPS=/path/to/lmp). DIRT is built `-C target-cpu=native`, so match it:
    cmake ../cmake -D BUILD_MPI=on -D PKG_GRANULAR=on -D CMAKE_BUILD_TYPE=Release \
          -D CMAKE_CXX_FLAGS="-O3 -mcpu=native"    # -march=native on x86
Knobs: PERF_NCORE (particles/core), PERF_RANKS, PERF_STEPS, PERF_NEIGH_EVERY, PERF_LAMMPS.
"""
import math, os, sys, re, json, shutil, random, statistics, subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "perf_mpi_scaling"
BIN = os.path.join(REPO_ROOT, "target", "release", "examples", EXAMPLE)
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
WORK = os.path.join(SCRIPT_DIR, "sweep", "scaling")
RESULTS = {m: os.path.join(DATA_DIR, f"{m}_scaling_results.json") for m in ("strong", "weak")}

# ── Physics (identical to the LAMMPS benchmark) ──────────────────────────────
R = 0.0011; D = 2 * R; RHO = 2500.0; SIG = 0.5
E = 5.0e7; NU = 0.3; ER = 0.9; MU = 0.3
DT = 6.3857607976e-06
STEPS = int(os.environ.get("PERF_STEPS", "5000"))
EVERY = int(os.environ.get("PERF_NEIGH_EVERY", "20"))
NPC = int(os.environ.get("PERF_NCORE", "10000"))         # particles per core
CORES = [int(x) for x in os.environ.get("PERF_RANKS", "1,2,4,8,16").split(",")]
DEC = {1: (1, 1, 1), 2: (2, 1, 1), 4: (2, 2, 1), 8: (2, 2, 2),
       16: (4, 2, 2), 32: (4, 4, 2), 64: (4, 4, 4)}
KB = 1.380649e-23; MASS = RHO * (4 / 3) * math.pi * R**3
TLMP = MASS * SIG**2 / KB; VP = (4 / 3) * math.pi * R**3
MPI = ["mpiexec", "--bind-to", "none"] + os.environ.get("PERF_MPI_EXTRA", "").split()

PHASE_CATS = {
    "force":     ["hertz_mindlin", "hooke_contact", "wall_contact_force"],
    "neighbor":  ["bin_neighbor_list", "sort_atoms_by_bin", "decide_rebuild"],
    "comm":      ["forward_comm_borders", "reverse_send_force", "comm::borders", "comm::exchange"],
    "rotation":  ["initial_rotation", "final_rotation"],
    "integrate": ["initial_integration", "final_integration"],
}
PHASES = ["force", "neighbor", "comm", "rotation", "integrate", "other"]
PCOL = {"force": "tab:red", "neighbor": "tab:green", "comm": "tab:blue",
        "rotation": "tab:orange", "integrate": "tab:purple", "other": "tab:gray"}


def phase_of(name):
    for p, keys in PHASE_CATS.items():
        if any(k in name for k in keys):
            return p
    return "other"


def parse_phases(out):
    """Median per-rank seconds per phase from the scheduler's per-system blocks.

    Under MPI each rank prints its own block (possibly truncated at exit); the
    median over blocks is the representative per-rank cost under a balanced
    decomposition."""
    per = {}
    for b in out.split("--- Per-system timing")[1:]:
        for ln in b.splitlines():
            if "TOTAL" in ln:
                break
            m = re.match(r"\s*(\S+)\s+([\d.]+)\s+[\d.]+%", ln)
            if m:
                per.setdefault(m.group(1), []).append(float(m.group(2)))
    t = {f"t_{p}": 0.0 for p in PHASES}
    for n, v in per.items():
        t[f"t_{phase_of(n)}"] += statistics.median(v)
    return t


def find_lammps():
    p = os.environ.get("PERF_LAMMPS")
    if p and os.path.exists(p):
        return p
    return shutil.which("lmp_mpi") or shutil.which("lmp")


def build():
    print("Building DIRT (release, native) ...")
    env = dict(os.environ, RUSTFLAGS="-C target-cpu=native")
    r = subprocess.run(["cargo", "build", "--release", "--example", EXAMPLE,
                        "--features", "mpi_backend"], cwd=REPO_ROOT, env=env)
    if r.returncode != 0 or not os.path.exists(BIN):
        sys.exit("DIRT build failed.")


# ── geometry: total N -> (dirt insert, lammps geom, box L, actual N) ─────────
def dense_geom(N, c):
    s = D * (0.5235987756 / 0.50) ** (1 / 3)
    ns = round(N ** (1 / 3)); L = ns * s
    csv = os.path.join(WORK, f"pos_{c}.csv"); random.seed(42)
    with open(csv, "w") as f:
        f.write("x,y,z,radius,vx,vy,vz\n")
        for i in range(ns):
            for j in range(ns):
                for k in range(ns):
                    v = [random.gauss(0, SIG) for _ in range(3)]
                    f.write(f"{(i+.5)*s:.8e},{(j+.5)*s:.8e},{(k+.5)*s:.8e},{R},{v[0]:.6e},{v[1]:.6e},{v[2]:.6e}\n")
    ins = (f'[[particles.insert]]\nsource="file"\nformat="csv"\nfile="{csv}"\nmaterial="glass"\nradius={R}\n'
           f'columns = {{ x=0, y=1, z=2, radius=3, vx=4, vy=5, vz=6 }}\ndensity={RHO}')
    geom = (f"region box block 0 {L:.8e} 0 {L:.8e} 0 {L:.8e} units box\ncreate_box 1 box\n"
            f"lattice sc {s:.8e} origin 0.5 0.5 0.5\ncreate_atoms 1 region box\n")
    return ins, geom, L, ns**3


def dilute_geom(N, c):
    L = (N * VP / 0.05) ** (1 / 3); pad = 0.05 * L
    ins = (f'[[particles.insert]]\nmaterial="glass"\ncount={N}\nradius={R}\ndensity={RHO}\nvelocity={SIG}\n'
           f'region={{ type="block", min=[{pad:.6e},{pad:.6e},{pad:.6e}], max=[{L-pad:.6e},{L-pad:.6e},{L-pad:.6e}] }}')
    geom = (f"region box block 0 {L:.8e} 0 {L:.8e} 0 {L:.8e} units box\ncreate_box 1 box\n"
            f"create_atoms 1 random {N} 12345 box overlap {D:.8e} maxtry 200000\n")
    return ins, geom, L, N


GEOM = {"dense": dense_geom, "dilute": dilute_geom}


def n_total(mode, c):
    # weak: fixed work per core; strong: fixed total (= top of the core ladder).
    return NPC * c if mode == "weak" else NPC * max(CORES)


def dirt_cfg(px, py, pz, ins, L, outd):
    return (f'[comm]\nprocessors_x={px}\nprocessors_y={py}\nprocessors_z={pz}\n[domain]\n'
            f'x_low=0.0\nx_high={L:.8e}\ny_low=0.0\ny_high={L:.8e}\nz_low=0.0\nz_high={L:.8e}\n'
            f'boundary_x="periodic"\nboundary_y="periodic"\nboundary_z="periodic"\n[neighbor]\n'
            f'skin_fraction=1.1\nevery={EVERY}\ncheck=false\nbin_size={2*D:.6e}\n[gravity]\n'
            f'gx=0.0\ngy=0.0\ngz=0.0\n[dem]\ncontact_model="hertz"\n[[dem.materials]]\nname="glass"\n'
            f'youngs_mod={E:.6e}\npoisson_ratio={NU}\nrestitution={ER}\nfriction={MU}\n{ins}\n'
            f'[output]\ndir="{outd}"\n[run]\nsteps={STEPS}\nthermo={STEPS}\ndt={DT:.10e}\n')


def run_dirt(px, py, pz, c, ins, L):
    outd = os.path.join(WORK, "dirt"); os.makedirs(outd, exist_ok=True)
    cfg = os.path.join(WORK, "c.toml"); open(cfg, "w").write(dirt_cfg(px, py, pz, ins, L, outd))
    cf = os.path.join(outd, "data", "perf_results.csv")
    if os.path.exists(cf):
        os.remove(cf)
    out = subprocess.run(MPI + ["-n", str(c), BIN, cfg], stdin=subprocess.DEVNULL,
                         capture_output=True, text=True).stdout
    if not os.path.exists(cf):
        return None
    row = open(cf).read().strip().splitlines()[-1].split(",")
    return {"tp": float(row[9]), **parse_phases(out)}


LMP_HEAD = "units si\natom_style sphere\natom_modify map array\nboundary p p p\nnewton off\ncomm_modify vel yes\n"


def run_lmp(lmp, px, py, pz, c, geom):
    if not lmp:
        return float("nan")
    inp = os.path.join(WORK, "in.lmp"); log = os.path.join(WORK, "l.log")
    tail = (f"set type 1 diameter {D:.8e}\nset type 1 density {RHO}\npair_style granular\n"
            f"pair_coeff 1 1 hertz/material {E:.6e} {ER} {NU} tangential mindlin NULL 1.0 {MU} "
            f"damping tsuji rolling none twisting none\nneighbor {0.1*D:.8e} bin\n"
            f"neigh_modify every {EVERY} delay 0 check no\n"
            f"velocity all create {TLMP:.6e} 98765 dist gaussian mom yes rot no\n"
            f"fix integrate all nve/sphere\ntimestep {DT:.10e}\nthermo {STEPS}\nrun {STEPS}\n")
    open(inp, "w").write(LMP_HEAD + f"processors {px} {py} {pz}\n" + geom + tail)
    if os.path.exists(log):
        os.remove(log)
    subprocess.run(MPI + ["-n", str(c), lmp, "-in", inp, "-log", log],
                   stdin=subprocess.DEVNULL, capture_output=True, text=True)
    try:
        txt = open(log).read()
        m = re.search(r"Loop time of ([\d.]+)", txt)
        a = re.search(r"with (\d+) atoms", txt)
        return (int(a.group(1)) * STEPS / float(m.group(1))) if m else float("nan")
    except OSError:
        return float("nan")


def run(mode):
    os.makedirs(WORK, exist_ok=True); os.makedirs(DATA_DIR, exist_ok=True)
    if not os.path.exists(BIN):
        build()
    lmp = find_lammps()
    print(f"\n=== {mode.upper()} scaling ===")
    print(f"DIRT: {BIN}\nLAMMPS: {lmp or 'not found (DIRT-only)'}")
    print(f"{NPC} particles/core, every={EVERY}, {STEPS} steps, "
          f"{'fixed total N' if mode == 'strong' else 'fixed N/core'}")
    print(f"{'regime':>7} {'cores':>5} {'N':>8} {'DIRT M/s':>9} {'LMP M/s':>8} {'D/L':>6}")
    results = {}
    for regime in ("dense", "dilute"):
        results[regime] = []
        for c in CORES:
            os.makedirs(WORK, exist_ok=True)   # robust if anything wiped it mid-run
            px, py, pz = DEC[c]
            ins, geom, L, N = GEOM[regime](n_total(mode, c), c)
            d = run_dirt(px, py, pz, c, ins, L)
            if not d:
                print(f"  ! DIRT {regime} c={c} failed"); continue
            l = run_lmp(lmp, px, py, pz, c, geom)
            rec = {"cores": c, "N": N, "dirt": d["tp"], "lmp": l,
                   **{k: d[k] for k in d if k.startswith("t_")}}
            results[regime].append(rec)
            dl = f"{d['tp']/l:.2f}" if l == l else "n/a"
            print(f"{regime:>7} {c:>5} {N:>8} {d['tp']/1e6:>9.1f} {l/1e6:>8.1f} {dl:>6}")
    json.dump(results, open(RESULTS[mode], "w"), indent=1)
    print(f"Wrote {RESULTS[mode]}")


def graph(mode):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter, NullFormatter
    if not os.path.exists(RESULTS[mode]):
        print(f"  (skip {mode}: no {RESULTS[mode]})"); return
    Rd = json.load(open(RESULTS[mode]))
    os.makedirs(PLOT_DIR, exist_ok=True)

    def cores(rg): return [r["cores"] for r in Rd[rg]]

    def logx(ax, cs):
        ax.set_xscale("log", base=2); ax.set_xticks(cs)
        f = ScalarFormatter(); f.set_scientific(False)
        ax.xaxis.set_major_formatter(f); ax.xaxis.set_minor_formatter(NullFormatter())
        ax.grid(True, which="both", ls=":", alpha=0.4)

    cs0 = cores("dense") if Rd.get("dense") else cores("dilute")
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))

    # 1. DIRT/LAMMPS throughput ratio
    ax = axes[0]; ax.axhline(1.0, color="k", ls="--", lw=1, label="parity")
    for rg, col in (("dense", "tab:blue"), ("dilute", "tab:orange")):
        if not Rd.get(rg):
            continue
        cs = cores(rg)
        ax.plot(cs, [r["dirt"] / r["lmp"] for r in Rd[rg]], "o-", color=col, lw=2, ms=7, label=f"{rg} gas")
    logx(ax, cs0)
    if Rd.get("dense"):
        for r in Rd["dense"]:
            ax.annotate(f"{r['dirt']/r['lmp']:.2f}x", (r["cores"], r["dirt"] / r["lmp"]),
                        textcoords="offset points", xytext=(0, 8), fontsize=8, ha="center")
    ax.set_xlabel(f"cores ({NPC//1000}k particles/core)" if mode == "weak" else "cores (fixed total N)")
    ax.set_ylabel("DIRT throughput / LAMMPS")
    ax.set_title("DIRT vs LAMMPS\n(>1 = DIRT faster)"); ax.legend()

    # 2. parallel efficiency (mode-appropriate)
    ax = axes[1]; ax.axhline(1.0, color="k", ls="--", lw=1, label="ideal")
    for rg, ls in (("dense", "-"), ("dilute", ":")):
        if not Rd.get(rg):
            continue
        cs = cores(rg); d = [r["dirt"] for r in Rd[rg]]; l = [r["lmp"] for r in Rd[rg]]
        if mode == "weak":   # throughput-per-core retention
            de = [(d[i] / cs[i]) / (d[0] / cs[0]) for i in range(len(cs))]
            le = [(l[i] / cs[i]) / (l[0] / cs[0]) for i in range(len(cs))]
        else:                # speedup / cores
            de = [(d[i] / d[0]) / (cs[i] / cs[0]) for i in range(len(cs))]
            le = [(l[i] / l[0]) / (cs[i] / cs[0]) for i in range(len(cs))]
        ax.plot(cs, de, "o" + ls, color="tab:red", lw=2, ms=6, label=f"DIRT {rg}")
        ax.plot(cs, le, "s" + ls, color="tab:gray", mfc="none", lw=1.6, ms=6, label=f"LAMMPS {rg}")
    logx(ax, cs0); ax.set_ylim(0, 1.08); ax.set_xlabel("cores")
    ax.set_ylabel("parallel efficiency")
    ax.set_title("Efficiency — DIRT vs LAMMPS"); ax.legend(fontsize=8)

    # 3. dense phase panel: breakdown (weak) | per-phase efficiency (strong)
    ax = axes[2]
    if Rd.get("dense"):
        cs = cores("dense")
        if mode == "weak":
            tot = [sum(r.get(f"t_{p}", 0.0) for p in PHASES) or 1.0 for r in Rd["dense"]]
            bottom = [0.0] * len(cs)
            for p in PHASES:
                vals = [100 * Rd["dense"][i].get(f"t_{p}", 0.0) / tot[i] for i in range(len(cs))]
                ax.bar([str(c) for c in cs], vals, bottom=bottom, color=PCOL[p], label=p)
                bottom = [b + v for b, v in zip(bottom, vals)]
            ax.set_ylabel("% of step time"); ax.set_title("Dense phase breakdown\n(comm grows with cores)")
            ax.legend(fontsize=8, ncol=2)
        else:
            ax.axhline(1.0, color="k", ls="--", lw=1, label="ideal")
            t = {p: [r.get(f"t_{p}", 0.0) for r in Rd["dense"]] for p in PHASES}
            for p in ("force", "neighbor", "comm", "rotation"):
                if t[p][0] <= 0:
                    continue
                eff = [t[p][0] / (cs[i] * t[p][i]) if t[p][i] > 0 else float("nan") for i in range(len(cs))]
                ax.plot(cs, eff, "o-", color=PCOL[p], lw=1.8, ms=6, label=p)
            logx(ax, cs); ax.set_ylim(0, 1.15)
            ax.set_ylabel("strong-scaling efficiency  (t₁/(cores·tₙ))")
            ax.set_title("Dense per-phase efficiency\n(force scales, comm collapses)")
            ax.legend(fontsize=9)
        ax.set_xlabel("cores")
    basis = (f"fixed total N≈{n_total('strong', max(CORES)):,}" if mode == "strong"
             else f"{NPC//1000}k particles/core")
    fig.suptitle(f"{mode.capitalize()} scaling — our gas configs, {basis}, "
                 f"every={EVERY}, DIRT(opt) vs LAMMPS(native)", fontsize=13)
    fig.tight_layout()
    out = os.path.join(PLOT_DIR, f"{mode}_scaling.png")
    fig.savefig(out, dpi=140); plt.close(fig)
    print(f"Wrote {out}")


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    if cmd in ("strong", "weak"):
        run(cmd); graph(cmd)
    elif cmd == "graph":
        for m in ("strong", "weak"):
            graph(m)
    elif cmd == "all":
        for m in ("strong", "weak"):
            run(m); graph(m)
    else:
        sys.exit("usage: scaling_gas.py [all|strong|weak|graph]")


if __name__ == "__main__":
    main()
