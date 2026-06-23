#!/usr/bin/env bash
# μ_r-sensitivity smoke: gentle settle at REAL e=0.926, μ_r ∈ {0.0, 0.5, 1.5}.
# Run from anywhere; paths are absolute. Build once, run 3 cases, fit θ_r.
set -u
ROOT="/Users/suehr/Documents/GitHub/dirt"
DIR="$ROOT/examples/SPH_glass_sphere_calibration/03_angle_of_repose"
TEMPLATE="$DIR/smoke_template.toml"
cd "$ROOT" || exit 1

echo "=== building sphcal_angle_of_repose ==="
cargo build --release --no-default-features --example sphcal_angle_of_repose 2>&1 | tail -5 || exit 1

for MUR in 0.0 0.5 1.5; do
  OUT="$DIR/data/smoke_mu${MUR}"
  CFG="$DIR/smoke_mu${MUR}.toml"
  mkdir -p "$OUT/data"
  sed -e "s|__MUR__|${MUR}|" -e "s|__OUTDIR__|${OUT}|" "$TEMPLATE" > "$CFG"
  echo "=== RUN mu_r=${MUR} ($(date +%H:%M:%S)) ==="
  cargo run --release --no-default-features --example sphcal_angle_of_repose -- "$CFG" 2>&1 \
    | grep -E "settled|rest|wrote|panic|Error|error" | tail -8
done

echo "=== FIT θ_r ==="
python3 - "$DIR" <<'PY'
import sys, csv, math, os
d = sys.argv[1]
def fit(path):
    if not os.path.exists(path): return None, 0
    xs=[]; ys=[]; zs=[]
    with open(path) as f:
        for row in csv.DictReader(f):
            xs.append(float(row['x'])); ys.append(float(row['y'])); zs.append(float(row['z']))
    n=len(xs)
    if n<20: return None, n
    r=[math.hypot(x,y) for x,y in zip(xs,ys)]
    rmax=sorted(r)[int(0.97*n)]
    nb=16; binz=[[] for _ in range(nb)]
    for ri,zi in zip(r,zs):
        b=min(nb-1,int(ri/rmax*nb))
        binz[b].append(zi)
    # surface = 90th pct height per bin (bin center radius)
    surf=[]
    for b in range(nb):
        if len(binz[b])>=3:
            zz=sorted(binz[b]); surf.append(((b+0.5)*rmax/nb, zz[int(0.9*len(zz))]))
    if len(surf)<4: return None, n
    # flank: from the bin of peak surface outward
    pk=max(range(len(surf)), key=lambda i:surf[i][1])
    flank=surf[pk:]
    if len(flank)<3: return None, n
    rr=[p[0] for p in flank]; ss=[p[1] for p in flank]
    m=len(rr); sx=sum(rr); sy=sum(ss); sxx=sum(a*a for a in rr); sxy=sum(a*b for a,b in zip(rr,ss))
    denom=(m*sxx-sx*sx)
    if denom==0: return None, n
    slope=(m*sxy-sx*sy)/denom
    return math.degrees(math.atan(abs(slope))), n
print(f"{'mu_r':>6} {'N':>5} {'theta_r(deg)':>12}")
for mur in ["0.0","0.5","1.5"]:
    th,n=fit(f"{d}/data/smoke_mu{mur}/data/repose_results.csv")
    print(f"{mur:>6} {n:>5} {('%.1f'%th) if th is not None else 'NA':>12}")
PY
echo "=== SMOKE DONE ($(date +%H:%M:%S)) ==="
