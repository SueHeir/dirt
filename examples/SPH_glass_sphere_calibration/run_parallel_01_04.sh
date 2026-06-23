#!/usr/bin/env bash
# Parallelize the two remaining sequential sweeps: #1 shear (20) + #4 enduring (7).
# 27 independent cases across the cores. Then graph both.
set -u
ROOT="/Users/suehr/Documents/GitHub/dirt"
SCD="$ROOT/examples/SPH_glass_sphere_calibration"
SHEAR="$ROOT/target/release/examples/sphcal_shear_rheology"
ENDUR="$ROOT/target/release/examples/sphcal_enduring_contact"
cd "$ROOT" || exit 1

echo "=== build both (no-op if cached) ==="
cargo build --release --no-default-features --example sphcal_shear_rheology  2>&1 | tail -1
cargo build --release --no-default-features --example sphcal_enduring_contact 2>&1 | tail -1

CFGS=$(ls "$SCD"/01_shear_rheology/sweep/*/config.toml "$SCD"/04_enduring_contact/sweep/*/config.toml 2>/dev/null)
N=$(echo "$CFGS" | grep -c config.toml)
echo "=== running $N cases, 14 at a time ($(date +%H:%M:%S)) ==="
echo "$CFGS" | xargs -P 14 -I{} sh -c '
  cfg="$1"
  case "$cfg" in
    *01_shear_rheology*)  bin="'"$SHEAR"'" ;;
    *04_enduring_contact*) bin="'"$ENDUR"'" ;;
  esac
  "$bin" "$cfg" > "$cfg.run.log" 2>&1 && echo "  done: $(basename $(dirname $cfg))"
' _ {}

echo "=== graph #1 shear ($(date +%H:%M:%S)) ==="
python3 "$SCD/01_shear_rheology/sweep.py"   graph 2>&1 | tail -16
echo "=== graph #4 enduring ($(date +%H:%M:%S)) ==="
python3 "$SCD/04_enduring_contact/sweep.py" graph 2>&1 | tail -16
echo "=== #1+#4 PARALLEL DONE ($(date +%H:%M:%S)) ==="
