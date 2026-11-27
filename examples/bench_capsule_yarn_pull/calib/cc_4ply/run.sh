#!/bin/bash
cd /Users/suehr/Documents/GitHub/dirt
D=examples/bench_capsule_yarn_pull/calib/cc_4ply
FIBRE_LIMIT=1000 CAPSULE_YARN_SPEC=examples/bench_capsule_yarn_pull/data/capsule_yarn_4ply_k8.txt \
CAPSULE_RELAX_STEPS=100000 CAPSULE_RELAX_DAMPING=30 BOND_BREAK_MIN_STEP=100000 DUMP_EVERY=10000 \
/usr/bin/time -p target/release/examples/bench_capsule_yarn_pull "$D/cfg.toml" > "$D/run.log" 2>&1
echo DONE >> "$D/run.log"
