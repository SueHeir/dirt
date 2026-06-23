#!/usr/bin/env bash
# Wait until the other production closures (#4, #6, #8) finish, then kill the
# sequential #1 (shear rheology) driver and re-run its 20 cases in PARALLEL.
set -u
ROOT="/Users/suehr/Documents/GitHub/dirt"
SCD="$ROOT/examples/SPH_glass_sphere_calibration"
BIN="$ROOT/target/release/examples/sphcal_shear_rheology"
cd "$ROOT" || exit 1

echo "=== watcher: waiting for #4/#6/#8 to finish ($(date +%H:%M:%S)) ==="
deadline=$(( $(date +%s) + 7200 ))   # 2h safety cap
while :; do
  a=$(pgrep -f "04_enduring_contact/sweep.py"     | wc -l | tr -d ' ')
  b=$(pgrep -f "06_conductivity/sweep.py"         | wc -l | tr -d ' ')
  c=$(pgrep -f "08_cooperativity_length/sweep.py" | wc -l | tr -d ' ')
  [ "$a" = 0 ] && [ "$b" = 0 ] && [ "$c" = 0 ] && break
  [ "$(date +%s)" -gt "$deadline" ] && { echo "watcher TIMEOUT"; break; }
  sleep 30
done
echo "=== others finished ($(date +%H:%M:%S)) — parallelizing #1 ==="

# Kill the sequential #1 driver + its current in-flight case.
pkill -f "01_shear_rheology/sweep.py" 2>/dev/null
pkill -f "examples/sphcal_shear_rheology" 2>/dev/null
sleep 2

CFGS=$(ls "$SCD"/01_shear_rheology/sweep/*/config.toml 2>/dev/null)
N=$(echo "$CFGS" | grep -c config.toml)
echo "=== running $N shear cases, 14 at a time ($(date +%H:%M:%S)) ==="
echo "$CFGS" | xargs -P 14 -I{} sh -c '"$1" "$2" > "$2.run.log" 2>&1 && echo "  done: $(dirname $2 | xargs basename)"' _ "$BIN" {}

echo "=== graphing #1 ($(date +%H:%M:%S)) ==="
python3 "$SCD/01_shear_rheology/sweep.py" graph 2>&1 | tail -20
echo "=== #1 PARALLEL DONE ($(date +%H:%M:%S)) ==="
