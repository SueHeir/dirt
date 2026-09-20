#!/usr/bin/env bash
set -eu

# Fifty percent in two million steps preserves the 10%-per-400k-step rate of
# the original paired experiment. DIRT derives the grip speed independently
# from each specimen's initial grip spacing.
export TENSION_SETTLE_STEPS=1000
export TENSION_RAMP_STEPS=5000
export TENSION_FINAL_STRAIN_PERCENT=50.0
export TENSION_RECORD_EVERY=2000

control_output="examples/fake_felted_tension/output_50pct"
needled_output="examples/fake_needled_tension/output_50pct"
mkdir -p "$control_output" "$needled_output"

mpirun -np 16 target/release/examples/fake_felted_tension \
  examples/fake_felted_tension/config_50pct.toml \
  >"$control_output/run.log" 2>&1 &
control_pid=$!

mpirun -np 16 target/release/examples/fake_felted_tension \
  examples/fake_needled_tension/config_50pct.toml \
  >"$needled_output/run.log" 2>&1 &
needled_pid=$!

trap 'kill "$control_pid" "$needled_pid" 2>/dev/null || true' INT TERM EXIT
wait "$control_pid"
wait "$needled_pid"
trap - INT TERM EXIT

python3 examples/fake_needled_tension/plot.py --fifty
