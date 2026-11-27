#!/usr/bin/env bash
set -euo pipefail

case_name="${1:?usage: hpc_run_50pct.sh <unneedled|needled> <cpu-set>}"
cpu_set="${2:?usage: hpc_run_50pct.sh <unneedled|needled> <cpu-set>}"

case "$case_name" in
  unneedled)
    config="examples/fake_felted_tension/config_50pct.toml"
    output="examples/fake_felted_tension/output_50pct"
    ;;
  needled)
    config="examples/fake_needled_tension/config_50pct.toml"
    output="examples/fake_needled_tension/output_50pct"
    ;;
  *)
    echo "unknown case: $case_name" >&2
    exit 2
    ;;
esac

if [[ -e "$output/tension.csv" ]]; then
  echo "refusing to overwrite existing $output/tension.csv" >&2
  exit 1
fi
mkdir -p "$output"

export TENSION_SETTLE_STEPS=1000
export TENSION_RAMP_STEPS=5000
export TENSION_FINAL_STRAIN_PERCENT=50.0
export TENSION_RECORD_EVERY=2000

echo "starting $case_name 50% pull on CPUs $cpu_set at $(date --iso-8601=seconds)"
echo "strain reference: initial internal non-gripped gauge length"
/usr/bin/time -p taskset -c "$cpu_set" mpirun --mca btl self,sm -np 16 \
  target/release/examples/fake_felted_tension "$config" 2>&1 \
  | tee "$output/run.log"
