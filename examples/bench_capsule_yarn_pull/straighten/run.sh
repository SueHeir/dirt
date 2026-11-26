#!/bin/bash
cd /Users/suehr/Documents/GitHub/dirt
D=examples/bench_capsule_yarn_pull/straighten
FIBRE_LIMIT=1 CAPSULE_RELAX_STEPS=15000000 CAPSULE_RELAX_DAMPING=6 BOND_BREAK_MIN_STEP=999999999 STRAIGHT_REST=1 DUMP_EVERY=150000 \
target/release/examples/bench_capsule_yarn_pull "$D/cfg.toml" > "$D/run.log" 2>&1
echo DONE >> "$D/run.log"
