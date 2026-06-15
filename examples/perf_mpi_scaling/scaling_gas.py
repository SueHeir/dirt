#!/usr/bin/env python3
"""Strong + weak scaling, DIRT-vs-LAMMPS, on controlled dense/dilute gas configs.

Directory-per-simulation workflow, designed to run interactively on a laptop OR as
PBS jobs on a supercomputer. Each (code, mode, regime, cores) point gets its own
self-contained directory under `runs/` holding its config, its matched initial
state, a `pbs_submit` script, and — after it runs — its `result.json`. The graph is
built by reading those directories, so partial results plot fine.

    runs/
      _init/                         # shared, matched init (one per (regime, N))
        dense_N157464.csv  .lmp
      dirt_strong_dense_c4/          # one simulation
        meta.json  config.toml  pbs_submit  ->  perf_1.csv run_1.log ...
      lammps_strong_dense_c4/
        meta.json  in.lammps  pbs_submit  ->  lammps_1.log run_1.log ...
      ...

Each pbs_submit is:  #PBS ... ; cd $PBS_O_WORKDIR ; <module load> ; mpirun -np N <bin> ...
(it mpiruns the binary directly -- no python at run time). The graph parses the raw
perf_*.csv / lammps_*.log / run_*.log left in each dir, so partial results plot fine.

Workflow:
    scaling_gas.py generate                 # build dirs, configs, pbs_submit, start_all.sh
    bash start_all.sh                        # PBS: qsub every job   (HPC)
    scaling_gas.py run-local                 # OR run every job here  (laptop)
    scaling_gas.py run-case runs/<dir>       # run a single dir locally (same mpirun block)
    scaling_gas.py graph                     # parse runs/*/ outputs -> plots/

Two 3-panel figures land in plots/: strong_scaling.png and weak_scaling.png.

Binaries: pass a directory OR a full path with PERF_DIRT / PERF_LAMMPS. LAMMPS, for a
fair comparison, must be built native and `newton on` (matched to DIRT's half list).
    cmake ../cmake -D BUILD_MPI=on -D PKG_GRANULAR=on -D CMAKE_BUILD_TYPE=Release \
          -D CMAKE_CXX_FLAGS="-O3 -mcpu=native"   # -march=native on x86

Knobs (env): PERF_RANKS, PERF_NCORE, PERF_STEPS, PERF_NEIGH_EVERY, PERF_REPS,
             PERF_DIRT, PERF_LAMMPS, PERF_MPI (launcher), PERF_MPI_NP, and PBS_* (below).
"""
import math, os, sys, re, json, glob, random, statistics, subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "perf_mpi_scaling"
RUNS_DIR = os.path.join(SCRIPT_DIR, "runs")
INIT_DIR = os.path.join(RUNS_DIR, "_init")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
START_ALL = os.path.join(SCRIPT_DIR, "start_all.sh")
MPI = os.environ.get("PERF_MPI", "mpirun")           # launcher; used as: <MPI> -np <cores> <bin> …
MPI_NP = os.environ.get("PERF_MPI_NP", "-np")        # rank flag (mpirun: -np, some MPTs: -n)


def _resolve_bin(value, names, default):
    """Accept a directory (find one of `names` in it) or a full path. Empty → default."""
    p = value or default
    if p and os.path.isdir(p):
        for n in names:
            cand = os.path.join(p, n)
            if os.path.exists(cand):
                return os.path.abspath(cand)
        return os.path.abspath(os.path.join(p, names[0]))   # record path even if not built yet
    return os.path.abspath(p) if p and os.path.sep in p else p


# Binary paths — give a directory OR a full path for each (env PERF_DIRT / PERF_LAMMPS).
BIN = _resolve_bin(os.environ.get("PERF_DIRT"), [EXAMPLE],
                   os.path.join(REPO_ROOT, "target", "release", "examples", EXAMPLE))
LMP = _resolve_bin(os.environ.get("PERF_LAMMPS"), ["lmp", "lmp_mpi", "lammps"], "lmp_mpi")

# ── Sweep definition ─────────────────────────────────────────────────────────
MODES = os.environ.get("PERF_MODES", "strong,weak").split(",")
REGIMES = os.environ.get("PERF_REGIMES", "dense,dilute").split(",")
CORES = [int(x) for x in os.environ.get("PERF_RANKS", "1,4,8,16,32,40").split(",")]
REPS = int(os.environ.get("PERF_REPS", "3"))
NPC = int(os.environ.get("PERF_NCORE", "10000"))     # particles per core
STEPS = int(os.environ.get("PERF_STEPS", "5000"))
EVERY = int(os.environ.get("PERF_NEIGH_EVERY", "20"))
DEC = {1: (1, 1, 1), 2: (2, 1, 1), 4: (2, 2, 1), 8: (2, 2, 2), 16: (4, 2, 2),
       32: (4, 4, 2), 40: (5, 4, 2), 64: (4, 4, 4), 128: (8, 4, 4)}


def decomp(c):
    """px*py*pz == c, balanced. Table for the common counts, factorize otherwise."""
    if c in DEC:
        return DEC[c]
    grid, m, f = [1, 1, 1], c, 2
    while m > 1:
        while m % f == 0:
            grid[grid.index(min(grid))] *= f; m //= f
        f += 1
    return tuple(sorted(grid, reverse=True))

# ── Physics (identical to the LAMMPS benchmark) ──────────────────────────────
R = 0.0011; D = 2 * R; RHO = 2500.0; SIG = 0.5
E = 5.0e7; NU = 0.3; ER = 0.9; MU = 0.3
DT = 6.3857607976e-06
KB = 1.380649e-23; MASS = RHO * (4 / 3) * math.pi * R**3
TLMP = MASS * SIG**2 / KB; VP = (4 / 3) * math.pi * R**3

# ── PBS / scheduler config (env-overridable; site-specific) ──────────────────
PBS_MODEL = os.environ.get("PBS_MODEL", "rom14")     # node model (e.g. NAS: bro/sky/rom14)
PBS_NCPUS = int(os.environ.get("PBS_NCPUS", "128"))  # physical cores per node of that model
PBS_GROUP = os.environ.get("PBS_GROUP", "GROUP_ID")  # -W group_list
PBS_QUEUE = os.environ.get("PBS_QUEUE", "normal")
PBS_WALLTIME = os.environ.get("PBS_WALLTIME", "00:30:00")
# Module-load block pasted into each job's pbs_submit. DIRT/Rust just needs MPI;
# LAMMPS often needs its own module(s) too. Both default to the shared PBS_MODULES.
PBS_MODULES = os.environ.get("PBS_MODULES", "YOUR module load HERE")          # e.g. "module load mpi-hpe/mpt"
PBS_MODULES_DIRT = os.environ.get("PBS_MODULES_DIRT", PBS_MODULES)
PBS_MODULES_LAMMPS = os.environ.get("PBS_MODULES_LAMMPS", PBS_MODULES)

# ── Per-phase timing buckets ─────────────────────────────────────────────────
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
    """Median per-rank seconds per phase from the scheduler's per-system blocks."""
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


# ── Matched initial state (shared per (regime, N)) ───────────────────────────
def n_total(mode, cores):
    # weak: fixed work per core; strong: fixed total (top of the core ladder).
    return NPC * cores if mode == "weak" else NPC * max(CORES)


def build_atoms(regime, N, seed=42):
    """Deterministic atom list (x,y,z,vx,vy,vz) + box L + actual N, identical for both
    codes. dense = sc lattice (φ=0.50); dilute = uniform random (φ=0.05)."""
    rng = random.Random(seed)
    if regime == "dense":
        s = D * (0.5235987756 / 0.50) ** (1 / 3)
        ns = round(N ** (1 / 3)); L = ns * s
        atoms = [((i + .5) * s, (j + .5) * s, (k + .5) * s, *[rng.gauss(0, SIG) for _ in range(3)])
                 for i in range(ns) for j in range(ns) for k in range(ns)]
    else:  # dilute
        L = (N * VP / 0.05) ** (1 / 3); pad = 0.05 * L
        atoms = [(rng.uniform(pad, L - pad), rng.uniform(pad, L - pad), rng.uniform(pad, L - pad),
                  *[rng.gauss(0, SIG) for _ in range(3)]) for _ in range(N)]
    return atoms, L, len(atoms)


def write_init(regime, N):
    """Write (once) the shared matched init for (regime, N): DIRT csv + LAMMPS data.
    Returns (csv_path, data_path, L, actual_N)."""
    os.makedirs(INIT_DIR, exist_ok=True)
    atoms, L, aN = build_atoms(regime, N)
    csv_p = os.path.join(INIT_DIR, f"{regime}_N{aN}.csv")
    lmp_p = os.path.join(INIT_DIR, f"{regime}_N{aN}.lmp")
    meta_p = os.path.join(INIT_DIR, f"{regime}_N{aN}.json")
    if os.path.exists(csv_p) and os.path.exists(lmp_p):
        return csv_p, lmp_p, L, aN
    with open(csv_p, "w") as f:
        f.write("x,y,z,radius,vx,vy,vz\n")
        for x, y, z, vx, vy, vz in atoms:
            f.write(f"{x:.8e},{y:.8e},{z:.8e},{R},{vx:.6e},{vy:.6e},{vz:.6e}\n")
    with open(lmp_p, "w") as f:
        f.write(f"# {regime} gas, matched init to DIRT\n\n{aN} atoms\n1 atom types\n\n")
        f.write(f"0.0 {L:.8e} xlo xhi\n0.0 {L:.8e} ylo yhi\n0.0 {L:.8e} zlo zhi\n\n")
        f.write("Atoms # sphere\n\n")               # id type diameter density x y z
        for n, (x, y, z, vx, vy, vz) in enumerate(atoms, 1):
            f.write(f"{n} 1 {D:.8e} {RHO} {x:.8e} {y:.8e} {z:.8e}\n")
        f.write("\nVelocities\n\n")                  # id vx vy vz wx wy wz
        for n, (x, y, z, vx, vy, vz) in enumerate(atoms, 1):
            f.write(f"{n} {vx:.6e} {vy:.6e} {vz:.6e} 0 0 0\n")
    json.dump({"L": L, "N": aN}, open(meta_p, "w"))
    return csv_p, lmp_p, L, aN


# ── Config / input / PBS generation ──────────────────────────────────────────
def dirt_config(px, py, pz, csv_path, L):
    ins = (f'[[particles.insert]]\nsource="file"\nformat="csv"\nfile="{csv_path}"\nmaterial="glass"\nradius={R}\n'
           f'columns = {{ x=0, y=1, z=2, radius=3, vx=4, vy=5, vz=6 }}\ndensity={RHO}')
    return (f'[comm]\nprocessors_x={px}\nprocessors_y={py}\nprocessors_z={pz}\n[domain]\n'
            f'x_low=0.0\nx_high={L:.8e}\ny_low=0.0\ny_high={L:.8e}\nz_low=0.0\nz_high={L:.8e}\n'
            f'boundary_x="periodic"\nboundary_y="periodic"\nboundary_z="periodic"\n[neighbor]\n'
            f'skin_fraction=1.1\nevery={EVERY}\ncheck=false\nbin_size={2*D:.6e}\n[gravity]\n'
            f'gx=0.0\ngy=0.0\ngz=0.0\n[dem]\ncontact_model="hertz"\n[[dem.materials]]\nname="glass"\n'
            f'youngs_mod={E:.6e}\npoisson_ratio={NU}\nrestitution={ER}\nfriction={MU}\n{ins}\n'
            f'[output]\ndir="."\n[run]\nsteps={STEPS}\nthermo={STEPS}\ndt={DT:.10e}\n')


def lammps_input(px, py, pz, data_path):
    return (f"units si\natom_style sphere\natom_modify map array\nboundary p p p\nnewton on\n"
            f"comm_modify vel yes\nprocessors {px} {py} {pz}\nread_data {data_path}\n"
            f"pair_style granular\n"
            f"pair_coeff 1 1 hertz/material {E:.6e} {ER} {NU} tangential mindlin NULL 1.0 {MU} "
            f"damping tsuji rolling none twisting none\n"
            f"neighbor {0.1*D:.8e} bin\nneigh_modify every {EVERY} delay 0 check no\n"
            f"fix integrate all nve/sphere\ntimestep {DT:.10e}\nthermo {STEPS}\nrun {STEPS}\n")


def run_block(code, cores):
    """Bash that mpiruns the binary REPS times into per-rep outputs. Pasted verbatim
    into pbs_submit (no python at run time) and also used by run-local."""
    launch = f'{MPI} {MPI_NP} {cores}'
    if code == "dirt":
        return (f'for r in $(seq 1 {REPS}); do\n'
                f'  rm -rf data\n'
                f'  {launch} "{BIN}" config.toml > "run_$r.log" 2>&1\n'
                f'  [ -f data/perf_results.csv ] && cp data/perf_results.csv "perf_$r.csv"\n'
                f'done\n')
    return (f'for r in $(seq 1 {REPS}); do\n'
            f'  {launch} "{LMP}" -in in.lammps -log "lammps_$r.log" > "run_$r.log" 2>&1\n'
            f'done\n')


def pbs_submit(code, name, cores):
    nodes = max(1, math.ceil(cores / PBS_NCPUS))
    per = min(cores, PBS_NCPUS)
    modules = PBS_MODULES_DIRT if code == "dirt" else PBS_MODULES_LAMMPS
    return (f"#!/bin/bash\n"
            f"#PBS -N {name}\n"
            f"#PBS -l select={nodes}:ncpus={per}:mpiprocs={per}:model={PBS_MODEL}\n"
            f"#PBS -l walltime={PBS_WALLTIME}\n"
            f"#PBS -q {PBS_QUEUE}\n"
            f"#PBS -W group_list={PBS_GROUP}\n"
            f"#PBS -j oe\n#PBS -o run.log\n\n"
            f"cd $PBS_O_WORKDIR\n"
            f"{modules}\n\n"
            f"{run_block(code, cores)}")


def case_name(code, mode, regime, cores):
    return f"{code}_{mode}_{regime}_c{cores}"


def generate():
    os.makedirs(RUNS_DIR, exist_ok=True)
    cases = []
    for mode in MODES:
        for regime in REGIMES:
            for c in CORES:
                csv_p, lmp_p, L, aN = write_init(regime, n_total(mode, c))
                px, py, pz = decomp(c)
                for code in ("dirt", "lammps"):
                    d = os.path.join(RUNS_DIR, case_name(code, mode, regime, c))
                    os.makedirs(d, exist_ok=True)
                    meta = {"code": code, "mode": mode, "regime": regime, "cores": c,
                            "N": aN, "px": px, "py": py, "pz": pz, "reps": REPS,
                            "bin": BIN, "lmp": LMP, "steps": STEPS}
                    json.dump(meta, open(os.path.join(d, "meta.json"), "w"), indent=1)
                    if code == "dirt":
                        open(os.path.join(d, "config.toml"), "w").write(dirt_config(px, py, pz, csv_p, L))
                    else:
                        open(os.path.join(d, "in.lammps"), "w").write(lammps_input(px, py, pz, lmp_p))
                    open(os.path.join(d, "pbs_submit"), "w").write(
                        pbs_submit(code, case_name(code, mode, regime, c), c))
                    cases.append(os.path.relpath(d, SCRIPT_DIR))
    with open(START_ALL, "w") as f:
        f.write("#!/bin/bash\n# Submit every benchmark job to PBS. Run after `generate`,\n"
                "# then `scaling_gas.py graph` once they finish.\n"
                'cd "$(dirname "$0")"\n')
        for c in cases:
            f.write(f'(cd "{c}" && qsub pbs_submit)\n')
    os.chmod(START_ALL, 0o755)
    print(f"Generated {len(cases)} job dirs under {RUNS_DIR}")
    print(f"  HPC:    bash {os.path.relpath(START_ALL, os.getcwd())}   # qsub all, then 'graph'")
    print(f"  local:  python3 {os.path.relpath(__file__, os.getcwd())} run-local")


# ── Run one case locally (same mpirun block the pbs script runs) ─────────────
def run_case(d):
    meta = json.load(open(os.path.join(d, "meta.json")))
    subprocess.run(["bash", "-c", run_block(meta["code"], meta["cores"])],
                   cwd=d, stdin=subprocess.DEVNULL)
    tp, _ = parse_dir(d, meta)
    msg = f"{tp/1e6:8.1f} M part-steps/s" if tp else "(no output)"
    print(f"  {os.path.basename(d):<34} {msg}")


def run_local():
    if not os.path.exists(BIN):
        print(f"DIRT binary not found at {BIN} — build it first:\n"
              f"  RUSTFLAGS='-C target-cpu=native' cargo build --release --example {EXAMPLE} --features mpi_backend")
        sys.exit(1)
    dirs = sorted(d for d in glob.glob(os.path.join(RUNS_DIR, "*", "")) if "_init" not in d)
    if not dirs:
        print("No job dirs — run 'generate' first."); sys.exit(1)
    print(f"Running {len(dirs)} cases locally ...")
    for d in dirs:
        run_case(d.rstrip(os.sep))


# ── Graph (parses raw mpirun outputs straight from runs/<dir>/) ──────────────
def parse_dir(d, meta):
    """-> (median throughput part-steps/s, phase dict). DIRT: perf_*.csv + run_*.log;
    LAMMPS: lammps_*.log. Returns (None, {}) if the job hasn't produced output yet."""
    tps, phases = [], []
    if meta["code"] == "dirt":
        for cf in sorted(glob.glob(os.path.join(d, "perf_*.csv"))):
            try:
                tps.append(float(open(cf).read().strip().splitlines()[-1].split(",")[9]))
            except (OSError, ValueError, IndexError):
                pass
        for lg in sorted(glob.glob(os.path.join(d, "run_*.log"))):
            ph = parse_phases(open(lg, errors="ignore").read())
            if sum(ph.values()) > 0:
                phases.append(ph)
    else:
        for lg in sorted(glob.glob(os.path.join(d, "lammps_*.log"))):
            try:
                txt = open(lg, errors="ignore").read()
                m = re.search(r"Loop time of ([\d.]+)", txt)
                a = re.search(r"with (\d+) atoms", txt)
                if m and a:
                    tps.append(int(a.group(1)) * meta["steps"] / float(m.group(1)))
            except OSError:
                pass
    if not tps:
        return None, {}
    pm = {}
    if phases:
        pm = {k: statistics.median(p[k] for p in phases) for k in phases[0]}
    return statistics.median(tps), pm


def load_results():
    """-> {mode: {regime: [records sorted by cores]}}, records have dirt+lmp+phases."""
    by = {}  # (mode, regime, cores) -> {}
    for mj in glob.glob(os.path.join(RUNS_DIR, "*", "meta.json")):
        d = os.path.dirname(mj)
        meta = json.load(open(mj))
        tp, ph = parse_dir(d, meta)
        if tp is None:
            continue
        key = (meta["mode"], meta["regime"], meta["cores"])
        rec = by.setdefault(key, {"cores": meta["cores"], "N": meta["N"]})
        if meta["code"] == "dirt":
            rec["dirt"] = tp
            rec.update(ph)
        else:
            rec["lmp"] = tp
    out = {}
    for (mode, regime, _), rec in by.items():
        out.setdefault(mode, {}).setdefault(regime, []).append(rec)
    for mode in out:
        for regime in out[mode]:
            out[mode][regime].sort(key=lambda r: r["cores"])
    return out


def graph():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter, NullFormatter
    data = load_results()
    if not data:
        sys.exit("No runs/*/result.json — run the jobs first (run-local or start_all.sh).")
    os.makedirs(PLOT_DIR, exist_ok=True)

    def logx(ax, cs):
        ax.set_xscale("log", base=2); ax.set_xticks(cs)
        f = ScalarFormatter(); f.set_scientific(False)
        ax.xaxis.set_major_formatter(f); ax.xaxis.set_minor_formatter(NullFormatter())
        ax.grid(True, which="both", ls=":", alpha=0.4)

    for mode in MODES:
        d = data.get(mode)
        if not d:
            continue
        any_regime = next(iter(d.values()))
        cs0 = [r["cores"] for r in any_regime]
        fig, axes = plt.subplots(1, 3, figsize=(17, 5))

        # 1. DIRT/LAMMPS ratio
        ax = axes[0]; ax.axhline(1.0, color="k", ls="--", lw=1, label="parity")
        for regime, col in (("dense", "tab:blue"), ("dilute", "tab:orange")):
            recs = [r for r in d.get(regime, []) if "dirt" in r and "lmp" in r]
            if not recs:
                continue
            cs = [r["cores"] for r in recs]
            ax.plot(cs, [r["dirt"] / r["lmp"] for r in recs], "o-", color=col, lw=2, ms=7, label=f"{regime} gas")
        logx(ax, cs0)
        for r in d.get("dense", []):
            if "dirt" in r and "lmp" in r:
                ax.annotate(f"{r['dirt']/r['lmp']:.2f}x", (r["cores"], r["dirt"] / r["lmp"]),
                            textcoords="offset points", xytext=(0, 8), fontsize=8, ha="center")
        ax.set_xlabel("cores"); ax.set_ylabel("DIRT throughput / LAMMPS")
        ax.set_title("DIRT vs LAMMPS\n(>1 = DIRT faster)"); ax.legend()

        # 2. parallel efficiency (mode-appropriate)
        ax = axes[1]; ax.axhline(1.0, color="k", ls="--", lw=1, label="ideal")
        for regime, ls in (("dense", "-"), ("dilute", ":")):
            recs = d.get(regime, [])
            if len(recs) < 2:
                continue
            cs = [r["cores"] for r in recs]
            for who, col in (("dirt", "tab:red"), ("lmp", "tab:gray")):
                if not all(who in r for r in recs):
                    continue
                t = [r[who] for r in recs]
                if mode == "weak":
                    eff = [(t[i] / cs[i]) / (t[0] / cs[0]) for i in range(len(cs))]
                else:
                    eff = [(t[i] / t[0]) / (cs[i] / cs[0]) for i in range(len(cs))]
                mk = "o" if who == "dirt" else "s"
                ax.plot(cs, eff, mk + ls, color=col, lw=2 if who == "dirt" else 1.6, ms=6,
                        mfc="none" if who == "lmp" else col, label=f"{who.upper()} {regime}")
        logx(ax, cs0); ax.set_ylim(0, 1.08); ax.set_xlabel("cores")
        ax.set_ylabel("parallel efficiency"); ax.set_title("Efficiency — DIRT vs LAMMPS")
        ax.legend(fontsize=8)

        # 3. dense phase panel: breakdown (weak) | per-phase efficiency (strong)
        ax = axes[2]; recs = d.get("dense", [])
        if recs and any("t_force" in r for r in recs):
            recs = [r for r in recs if "t_force" in r]
            cs = [r["cores"] for r in recs]
            if mode == "weak":
                tot = [sum(r.get(f"t_{p}", 0.0) for p in PHASES) or 1.0 for r in recs]
                bottom = [0.0] * len(cs)
                for p in PHASES:
                    vals = [100 * recs[i].get(f"t_{p}", 0.0) / tot[i] for i in range(len(cs))]
                    ax.bar([str(c) for c in cs], vals, bottom=bottom, color=PCOL[p], label=p)
                    bottom = [b + v for b, v in zip(bottom, vals)]
                ax.set_ylabel("% of step time"); ax.set_title("Dense phase breakdown")
                ax.legend(fontsize=8, ncol=2)
            else:
                ax.axhline(1.0, color="k", ls="--", lw=1, label="ideal")
                t = {p: [r.get(f"t_{p}", 0.0) for r in recs] for p in PHASES}
                for p in ("force", "neighbor", "comm", "rotation"):
                    if t[p][0] <= 0:
                        continue
                    eff = [t[p][0] / (cs[i] * t[p][i]) if t[p][i] > 0 else float("nan") for i in range(len(cs))]
                    ax.plot(cs, eff, "o-", color=PCOL[p], lw=1.8, ms=6, label=p)
                logx(ax, cs); ax.set_ylim(0, 1.15)
                ax.set_ylabel("strong-scaling efficiency  (t₁/(cores·tₙ))")
                ax.set_title("Dense per-phase efficiency")
                ax.legend(fontsize=9)
            ax.set_xlabel("cores")
        fig.suptitle(f"{mode.capitalize()} scaling — dense/dilute gas, {NPC//1000}k/core, "
                     f"every={EVERY}, DIRT(opt) vs LAMMPS(native, newton on)", fontsize=13)
        fig.tight_layout()
        out = os.path.join(PLOT_DIR, f"{mode}_scaling.png")
        fig.savefig(out, dpi=140); plt.close(fig)
        print(f"Wrote {out}")


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "help"
    if cmd == "generate":
        generate()
    elif cmd == "run-case":
        run_case(os.path.abspath(sys.argv[2]))
    elif cmd == "run-local":
        run_local(); graph()
    elif cmd == "graph":
        graph()
    elif cmd == "all":            # laptop convenience: generate + run-local + graph
        generate(); run_local(); graph()
    else:
        print(__doc__)
        sys.exit(0 if cmd in ("help", "-h", "--help") else 1)


if __name__ == "__main__":
    main()
