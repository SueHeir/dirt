#!/usr/bin/env bash
set -eu

# Matched duration and target strain give the two specimens the same imposed
# engineering strain rate despite their slightly different initial thicknesses.
export TENSION_SETTLE_STEPS=1000
export TENSION_RAMP_STEPS=5000
export TENSION_FINAL_STRAIN_PERCENT=1.0
export TENSION_RECORD_EVERY=100

mpirun -np 8 target/release/examples/fake_felted_tension \
  examples/fake_felted_tension/config_slow.toml &
control_pid=$!

mpirun -np 8 target/release/examples/fake_felted_tension \
  examples/fake_needled_tension/config_slow.toml &
needled_pid=$!

trap 'kill "$control_pid" "$needled_pid" 2>/dev/null || true' INT TERM EXIT
wait "$control_pid"
wait "$needled_pid"
trap - INT TERM EXIT

python3 examples/fake_needled_tension/plot.py --slow
