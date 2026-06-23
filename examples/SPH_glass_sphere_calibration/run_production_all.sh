#!/usr/bin/env bash
# Run all PRODUCTION closures with the finalized canonical material (μ_r=0.10,
# ±10% poly). Build once, then run each folder's own generate→run→graph pipeline
# concurrently (folders run in parallel; cases within a folder are sequential).
# Validation runs (KT in bench_lebc_shear, Haff in 05/config.toml) are NOT re-run.
set -u
ROOT="/Users/suehr/Documents/GitHub/dirt"
SCD="$ROOT/examples/SPH_glass_sphere_calibration"
cd "$ROOT" || exit 1

EXES="sphcal_shear_rheology sphcal_compressibility sphcal_enduring_contact \
sphcal_cooling_dissipation sphcal_conductivity sphcal_cooperativity_length"

echo "=== build all production examples once ($(date +%H:%M:%S)) ==="
for ex in $EXES; do
  echo "  building $ex"
  cargo build --release --no-default-features --example "$ex" 2>&1 | tail -1
done

echo "=== launch production closures concurrently ($(date +%H:%M:%S)) ==="
python3 "$SCD/01_shear_rheology/sweep.py"        all   > "$SCD/01_shear_rheology/prod.log"        2>&1 &
python3 "$SCD/02_compressibility/sweep.py"       all   > "$SCD/02_compressibility/prod.log"       2>&1 &
python3 "$SCD/04_enduring_contact/sweep.py"      all   > "$SCD/04_enduring_contact/prod.log"      2>&1 &
python3 "$SCD/06_conductivity/sweep.py"          all   > "$SCD/06_conductivity/prod.log"          2>&1 &
python3 "$SCD/08_cooperativity_length/sweep.py"  --run > "$SCD/08_cooperativity_length/prod.log"  2>&1 &
# 05 production variant (Haff validation left untouched).
"$ROOT/target/release/examples/sphcal_cooling_dissipation" \
    "$SCD/05_cooling_dissipation/config_production.toml" \
    > "$SCD/05_cooling_dissipation/prod_production.log" 2>&1 &

wait
echo "=== ALL PRODUCTION DONE ($(date +%H:%M:%S)) ==="
echo "--- per-folder tail ---"
for d in 01_shear_rheology 02_compressibility 04_enduring_contact 06_conductivity 08_cooperativity_length; do
  echo "### $d"; tail -3 "$SCD/$d/prod.log" 2>/dev/null
done
echo "### 05_cooling_dissipation (production)"; tail -3 "$SCD/05_cooling_dissipation/prod_production.log" 2>/dev/null
