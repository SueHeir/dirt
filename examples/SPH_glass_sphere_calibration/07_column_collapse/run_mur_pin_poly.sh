#!/usr/bin/env bash
# Narrow POLYDISPERSE (±10%) re-pin of μ_r against the Lube/Lajeunesse runout law.
# Mono baseline bracketed μ_r≈0.1; this refines it on the production material.
set -u
DIR="/Users/suehr/Documents/GitHub/dirt/examples/SPH_glass_sphere_calibration/07_column_collapse"
cd "$DIR" || exit 1
export POLY=0.10
MURS="0.05 0.1 0.15"

for MUR in $MURS; do
  echo "############ μ_r = ${MUR}  POLY=${POLY}  ($(date +%H:%M:%S)) ############"
  MU_R=$MUR python3 sweep.py generate
  MU_R=$MUR python3 sweep.py start 2>&1 | grep -E "Building|a=|wrote|ERROR|no deposit" | tail -12
done

echo "============ AGGREGATE (POLY ${POLY}): runout vs Lube/Lajeunesse ============"
python3 - "$DIR" "$MURS" <<'PY'
import sys, csv, os
d=sys.argv[1]; murs=sys.argv[2].split()
def law(a): return 1.2*a if a<=3.0 else 1.6*a**(2.0/3.0)
data={}
for m in murs:
    p=f"{d}/data/runout_mu{m}.csv"
    data[m]={float(r['aspect']):float(r['runout_norm']) for r in csv.DictReader(open(p))} if os.path.exists(p) else {}
aspects=sorted({a for m in murs for a in data[m]})
print(f"{'a':>4} {'exp':>7} " + " ".join(f"mu={m:>4}" for m in murs))
for a in aspects:
    print(f"{a:>4} {law(a):>7.2f} " + " ".join(f"{(('%.2f'%data[m][a]) if a in data[m] else 'NA'):>7}" for m in murs))
print("\nmean |rel.err| vs law (all aspects, and excl. squat a<=0.5):")
best=None
for m in murs:
    all_e=[abs(data[m][a]-law(a))/law(a) for a in aspects if a in data[m]]
    cl_e=[abs(data[m][a]-law(a))/law(a) for a in aspects if a in data[m] and a>0.5]
    if not all_e: print(f"  mu_r={m}: NA"); continue
    ma=sum(all_e)/len(all_e); mc=sum(cl_e)/len(cl_e)
    print(f"  mu_r={m}: all={ma:.3f}  clean={mc:.3f}")
    if best is None or mc<best[1]: best=(m,mc)
if best: print(f"\n==> best (clean-regime) match: mu_r = {best[0]} (rel.err {best[1]:.3f})")
PY
echo "============ POLY MUR PIN DONE ($(date +%H:%M:%S)) ============"
