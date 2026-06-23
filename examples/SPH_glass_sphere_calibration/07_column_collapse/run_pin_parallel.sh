#!/usr/bin/env bash
# PARALLEL polydisperse (±10%) μ_r re-pin. 18 independent runs (3 μ_r × 6 aspects)
# run concurrently (each DEM run is single-threaded; box has 6 perf + 12 eff cores).
set -u
ROOT="/Users/suehr/Documents/GitHub/dirt"
DIR="$ROOT/examples/SPH_glass_sphere_calibration/07_column_collapse"
BIN="$ROOT/target/release/examples/sphcal_column_collapse"
cd "$ROOT" || exit 1
export POLY=0.10
MURS="0.05 0.1 0.15"
NPAR=12

echo "=== generate configs ($(date +%H:%M:%S)) ==="
for MUR in $MURS; do MU_R=$MUR python3 "$DIR/sweep.py" generate >/dev/null; done

echo "=== build once ==="
cargo build --release --no-default-features --example sphcal_column_collapse 2>&1 | tail -2

# Collect configs, wipe stale deposits.
CFGS=$(for MUR in $MURS; do ls "$DIR"/sweep_mu${MUR}/*/config.toml; done)
echo "$CFGS" | sed 's#/config.toml#/data/column_collapse_results.csv#' | xargs rm -f 2>/dev/null

N=$(echo "$CFGS" | wc -l | tr -d ' ')
echo "=== running $N cases, $NPAR at a time ($(date +%H:%M:%S)) ==="
echo "$CFGS" | xargs -P "$NPAR" -I{} sh -c '"$1" "$2" > "$2.run.log" 2>&1 && echo "  done: $2"' _ "$BIN" {}

echo "=== AGGREGATE (POLY ${POLY}) ($(date +%H:%M:%S)) ==="
python3 - "$DIR" "$MURS" <<'PY'
import sys, os, importlib.util
d=sys.argv[1]; murs=sys.argv[2].split()
spec=importlib.util.spec_from_file_location("sweepmod", f"{d}/sweep.py")
sw=importlib.util.module_from_spec(spec); spec.loader.exec_module(sw)
def law(a): return 1.2*a if a<=3.0 else 1.6*a**(2.0/3.0)
ASPECTS=[0.5,1.0,2.0,3.0,4.0,5.0]
def tag(a): return ("a%g"%a).replace(".","p")
data={}
for m in murs:
    data[m]={}
    for a in ASPECTS:
        dep=f"{d}/sweep_mu{m}/{tag(a)}/data/column_collapse_results.csv"
        if os.path.exists(dep):
            _,lf=sw.measure_column(dep); data[m][a]=(lf-sw.L0)/sw.L0
    # write per-mu runout csv
    import csv
    with open(f"{d}/data/runout_mu{m}.csv","w",newline="") as f:
        w=csv.writer(f); w.writerow(["aspect","runout_norm"])
        for a in ASPECTS:
            if a in data[m]: w.writerow([a,data[m][a]])
print(f"{'a':>4} {'exp':>7} " + " ".join(f"mu={m:>5}" for m in murs))
for a in ASPECTS:
    print(f"{a:>4} {law(a):>7.2f} " + " ".join(f"{(('%.2f'%data[m][a]) if a in data[m] else 'NA'):>8}" for m in murs))
print("\nmean |rel.err| vs law (all aspects | clean a>0.5):")
best=None
for m in murs:
    ae=[abs(data[m][a]-law(a))/law(a) for a in ASPECTS if a in data[m]]
    ce=[abs(data[m][a]-law(a))/law(a) for a in ASPECTS if a in data[m] and a>0.5]
    if not ae: print(f"  mu_r={m}: NA"); continue
    ma=sum(ae)/len(ae); mc=sum(ce)/len(ce)
    print(f"  mu_r={m}: all={ma:.3f}  clean={mc:.3f}")
    if best is None or mc<best[1]: best=(m,mc)
if best: print(f"\n==> best (clean-regime) match: mu_r = {best[0]} (rel.err {best[1]:.3f})")
PY
echo "=== PARALLEL PIN DONE ($(date +%H:%M:%S)) ==="
