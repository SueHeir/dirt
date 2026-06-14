#!/usr/bin/env bash
# Run the hopper quiescence benchmark: baseline (no optimization) vs coherence.
#   A. performance  — wall time split by phase (filling vs emptying)
#   B. validation   — settled fill height (baseline vs coherence)
#   C. validation   — emptying speed / discharge curve (baseline vs coherence)
# Usage: ./examples/hopper_quiescence/bench.sh   (runs from any directory)
set -euo pipefail

cd "$(dirname "$0")/../.."          # repo root
DIR=examples/hopper_quiescence
FILL_STEP=120000                    # last step of the "filling" stage (gate opens next)

cargo build --release --no-default-features --example hopper_quiescence

for v in baseline coherence; do
    echo "=== running $v ==="
    ./target/release/examples/hopper_quiescence "$DIR/config_$v.toml" 2>&1 |
        grep -E "gate opened|FINAL|TOTAL WALL"
done

echo
echo "=== A. performance — wall time by phase [s] ==="
printf '%-12s %8s %8s %8s %9s\n' variant fill empty total speedup
base_total=""
for v in baseline coherence; do
    read -r fill total < <(awk -F, -v fs="$FILL_STEP" '
        NR>1 { last=$2; if ($1==fs) fill=$2 }
        END  { print fill, last }' "$DIR/config_${v}_stats.csv")
    [ "$v" = baseline ] && base_total=$total
    awk -v V="$v" -v f="$fill" -v t="$total" -v bt="$base_total" \
        'BEGIN { printf "%-12s %8.1f %8.1f %8.1f %8.2fx\n", V, f, t-f, t, bt/t }'
done

echo
echo "=== B. validation: settled fill height [m] (max particle z at step $FILL_STEP) ==="
printf '%-12s %10s\n' variant top_z
for v in baseline coherence; do
    awk -F, -v V="$v" -v fs="$FILL_STEP" \
        '$1==fs { printf "%-12s %10.5f\n", V, $11 }' "$DIR/config_${v}_stats.csv"
done

echo
echo "=== C. validation: emptying speed — n_discharged vs step (flow stage) ==="
paste -d' ' \
    <(awk -F, -v fs="$FILL_STEP" 'NR>1 && $1>=fs && $1%10000==0 {print $1, $5}' "$DIR/config_baseline_stats.csv") \
    <(awk -F, -v fs="$FILL_STEP" 'NR>1 && $1>=fs && $1%10000==0 {print $5}' "$DIR/config_coherence_stats.csv") |
    awk 'BEGIN{print "step baseline coherence"} {print}'

echo "  step to 90% discharged (n_discharged >= 18000):"
for v in baseline coherence; do
    awk -F, -v V="$v" 'NR>1 && $5>=18000 {printf "    %-12s %d\n", V, $1; exit}' "$DIR/config_${v}_stats.csv"
done
