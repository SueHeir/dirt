#!/usr/bin/env bash
# μ_r pinning sweep at REAL e=0.926, polydisperse. μ_r × 2 seeded reps. Fit θ_r.
set -u
ROOT="/Users/suehr/Documents/GitHub/dirt"
DIR="$ROOT/examples/SPH_glass_sphere_calibration/03_angle_of_repose"
TEMPLATE="$DIR/pin_template.toml"
cd "$ROOT" || exit 1
MURS="0.1 0.2 0.3 0.4 0.6"
SEEDS="1 2"

echo "=== build ==="
cargo build --release --no-default-features --example sphcal_angle_of_repose 2>&1 | tail -3 || exit 1

for MUR in $MURS; do
  for SEED in $SEEDS; do
    OUT="$DIR/data/pin_mu${MUR}_s${SEED}"
    CFG="$DIR/pin_mu${MUR}_s${SEED}.toml"
    mkdir -p "$OUT/data"
    sed -e "s|__MUR__|${MUR}|" -e "s|__SEED__|${SEED}|" -e "s|__OUTDIR__|${OUT}|" "$TEMPLATE" > "$CFG"
    echo "=== mu_r=${MUR} seed=${SEED} ($(date +%H:%M:%S)) ==="
    cargo run --release --no-default-features --example sphcal_angle_of_repose -- "$CFG" 2>&1 \
      | grep -E "wrote|panic|error|Error" | tail -3
  done
done

echo "=== FIT θ_r (mean ± half-spread over reps) ==="
python3 - "$DIR" "$MURS" "$SEEDS" <<'PY'
import sys, csv, math, os
d=sys.argv[1]; murs=sys.argv[2].split(); seeds=sys.argv[3].split()
def fit(path):
    if not os.path.exists(path): return None
    xs=[];ys=[];zs=[]
    with open(path) as f:
        for r in csv.DictReader(f):
            xs.append(float(r['x']));ys.append(float(r['y']));zs.append(float(r['z']))
    n=len(xs)
    if n<30: return None
    r=[math.hypot(x,y) for x,y in zip(xs,ys)]
    rmax=sorted(r)[int(0.97*n)]
    nb=20; binz=[[] for _ in range(nb)]
    for ri,zi in zip(r,zs):
        b=min(nb-1,int(ri/rmax*nb)); binz[b].append(zi)
    surf=[]
    for b in range(nb):
        if len(binz[b])>=3:
            zz=sorted(binz[b]); surf.append(((b+0.5)*rmax/nb, zz[int(0.85*len(zz))]))
    # clean flank: radius in [0.2,0.85]*rmax (skip rounded apex + toe)
    flank=[(rr,ss) for rr,ss in surf if 0.2*rmax<=rr<=0.85*rmax]
    if len(flank)<4: return None
    rr=[p[0] for p in flank]; ss=[p[1] for p in flank]
    m=len(rr);sx=sum(rr);sy=sum(ss);sxx=sum(a*a for a in rr);sxy=sum(a*b for a,b in zip(rr,ss))
    den=m*sxx-sx*sx
    if den==0: return None
    return math.degrees(math.atan(abs((m*sxy-sx*sy)/den)))
print(f"{'mu_r':>6} {'theta_r(deg) per rep':>26} {'mean':>7}")
for mur in murs:
    vals=[fit(f"{d}/data/pin_mu{mur}_s{s}/data/repose_results.csv") for s in seeds]
    vv=[v for v in vals if v is not None]
    mean=sum(vv)/len(vv) if vv else float('nan')
    s=' '.join(('%.1f'%v) if v is not None else 'NA' for v in vals)
    flag=' <-- IN BAND' if 22.0<=mean<=26.0 else ''
    print(f"{mur:>6} {s:>26} {mean:>7.1f}{flag}")
PY
echo "=== PIN DONE ($(date +%H:%M:%S)) ==="
