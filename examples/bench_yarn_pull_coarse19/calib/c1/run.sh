#!/bin/bash
cd /Users/suehr/Documents/GitHub/dirt
D=examples/bench_yarn_pull_coarse19/calib/c1
/usr/bin/time -p target/release/examples/bench_yarn_pull_coarse19 "$D/cfg.toml" > "$D/run.log" 2>&1
echo DONE >> "$D/run.log"
