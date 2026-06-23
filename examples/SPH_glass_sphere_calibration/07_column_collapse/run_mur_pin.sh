#!/usr/bin/env bash
# Pin μ_r: run the full aspect-ratio column-collapse sweep at several μ_r, then
# compare each μ_r's dimensionless runout curve to the Lube/Lajeunesse law
# ( (L_f-L0)/L0 ≈ 1.2·a for a≲3 ; ≈ 1.6·a^(2/3) for a≳3 ). Smooth spheres
# (μ_r=0) over-run; the pinned μ_r is the one whose curve best matches the law.
set -u
DIR="/Users/suehr/Documents/GitHub/dirt/examples/SPH_glass_sphere_calibration/07_column_collapse"
cd "$DIR" || exit 1
MURS="0.0 0.1 0.2 0.3"

for MUR in $MURS; do
  echo "############ μ_r = ${MUR}  ($(date +%H:%M:%S)) ############"
  MU_R=$MUR python3 sweep.py generate
  MU_R=$MUR python3 sweep.py start 2>&1 | grep -E "Building|a=|wrote|ERROR|no deposit" | tail -12
done

echo "============ AGGREGATE: runout vs Lube/Lajeunesse ============"
python3 - "$DIR" "$MURS" <<'PY'
import sys, csv, os, math
d=sys.argv[1]; murs=sys.argv[2].split()
def law(a):  # experimental dimensionless runout
    return 1.2*a if a<=3.0 else 1.6*a**(2.0/3.0)
print(f"{'a':>4} {'exp':>7} " + " ".join(f"mu={m:>4}" for m in murs))
# load
data={}
for m in murs:
    p=f"{d}/data/runout_mu{m}.csv"
    if not os.path.exists(p): data[m]={}; continue
    data[m]={float(r['aspect']):float(r['runout_norm']) for r in csv.DictReader(open(p))}
aspects=sorted({a for m in murs for a in data[m]})
for a in aspects:
    row=f"{a:>4} {law(a):>7.2f} "
    for m in murs:
        v=data[m].get(a); row+=f"{('%.2f'%v) if v is not None else 'NA':>7} "
    print(row)
# pin metric: mean |relative error| vs law across aspects
print("\nmean |rel.err| vs experimental law (lower = better):")
best=None
for m in murs:
    errs=[abs(data[m][a]-law(a))/law(a) for a in aspects if a in data[m]]
    if not errs: print(f"  mu_r={m}: NA"); continue
    me=sum(errs)/len(errs)
    print(f"  mu_r={m}: {me:.3f}")
    if best is None or me<best[1]: best=(m,me)
if best: print(f"\n==> best match: mu_r = {best[0]} (mean rel.err {best[1]:.3f})")
PY
echo "============ MUR PIN DONE ($(date +%H:%M:%S)) ============"
