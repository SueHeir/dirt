#!/bin/bash
cd /Users/suehr/Documents/GitHub/dirt
D=examples/bench_yarn_pull_coarse19/calib/4ply_mpi
export YARN_RELAX_STEPS=100000 YARN_RELAX_DAMPING=30 BOND_BREAK_MIN_STEP=100000 YARN_BREAK_LOG=1 DUMP_EVERY=15000
/usr/bin/time -p mpirun -np 6 target/release/examples/bench_yarn_pull_coarse19 "$D/cfg.toml" > "$D/run_dump.log" 2>&1
python3 examples/bench_yarn_pull_coarse19/merge_traj.py "$D" >> "$D/run_dump.log" 2>&1
echo DONE >> "$D/run_dump.log"
