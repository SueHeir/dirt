#!/usr/bin/env python3
"""Fixed-free bonded-fiber timestep experiment after Guo et al. (2013).

The independently cited reference is Eq. 41, dt=sqrt(2/3)Lsqrt(rho/E).
Unlike the superseded one-bond oscillator check, this launches an actual DIRT
wave/bending pulse along a many-sphere fixed-free fiber.  The axial high mode
has the mesh-scale wavelength that produces the cited wave-transit bound.
"""
from __future__ import annotations
import csv, math, os, shutil, subprocess, sys, tomllib
import numpy as np
from pathlib import Path

ROOT=Path(__file__).resolve().parent; REPO=ROOT.parents[1]
DATA,PLOTS,RUNS=ROOT/'data',ROOT/'plots',ROOT/'runs'
EXE=REPO/'target/release/examples/bench_fiber_timestep'

def load():
    with (ROOT/'config.toml').open('rb') as f: return tomllib.load(f)
def dt_guo(m): return math.sqrt(2/3)*m['bond_length']*math.sqrt(m['density']/m['youngs_modulus'])

def lattice_limits(m):
    """Independent small-amplitude fixed-free lattice limits for this probe.

    The axial element energy is k_n/2 (u_j-u_i)^2.  The transverse element
    energy is k_s/2 (z_j-z_i+L(theta_i+theta_j)/2)^2 +
    k_b/2 (theta_j-theta_i)^2.  With node 0 fixed, the largest eigenvalue of
    M^-1/2 K M^-1/2 gives omega_max and the explicit central-difference bound
    dt=2/omega_max.  This includes the translational DOF coupled through the
    shear lever arm; it is not a fitted DIRT threshold.
    """
    n, E, G, rho, r, rb, L = (m[k] for k in ('nodes','youngs_modulus','shear_modulus','density','sphere_radius','bond_radius','bond_length'))
    A, I = math.pi*rb**2, math.pi*rb**4/4
    mass, moi = rho*4*math.pi*r**3/3, 0.4*rho*4*math.pi*r**5/3
    ka, ks, kb = E*A/L, G*A/L, E*I/L
    axial = np.zeros((n-1,n-1)); transverse = np.zeros((2*(n-1),2*(n-1)))
    def add(K, entries, k):
        for i, a in entries:
            for j, b in entries: K[i,j] += k*a*b
    for e in range(n-1):
        axial_entries=[]; shear_entries=[]; bend_entries=[]
        for node, sign in ((e,-1.0),(e+1,1.0)):
            if node:
                axial_entries.append((node-1,sign))
                shear_entries += [(2*(node-1),sign),(2*(node-1)+1,L/2)]
                bend_entries.append((2*(node-1)+1,sign))
        add(axial, axial_entries, ka); add(transverse, shear_entries, ks); add(transverse, bend_entries, kb)
    def limit(K, diagonal_mass):
        invsqrt=np.diag([1/math.sqrt(x) for x in diagonal_mass])
        omega=math.sqrt(float(np.linalg.eigvalsh(invsqrt@K@invsqrt)[-1]))
        return 2/omega
    return {'axial':limit(axial,[mass]*(n-1)), 'bending':limit(transverse,[mass,moi]*(n-1))}
def geometry(m):
    RUNS.mkdir(exist_ok=True); p=RUNS/'fiber.csv'
    p.write_text(''.join(f'{i*m["bond_length"]:.12e},0,0\n' for i in range(m['nodes'])))
    return p
def case(m,mode,frac,steps):
    L=m['bond_length']; out=RUNS/f'{mode}_{frac:.4f}'; p=RUNS/f'{mode}_{frac:.4f}.toml'; geo=geometry(m)
    p.write_text(f'''[comm]
processors_x=1
processors_y=1
processors_z=1
[domain]
x_low={-2*L:.12e}
x_high={(m['nodes']+2)*L:.12e}
y_low={-4*L:.12e}
y_high={4*L:.12e}
z_low={-4*L:.12e}
z_high={4*L:.12e}
boundary_x="fixed"
boundary_y="fixed"
boundary_z="fixed"
[neighbor]
skin_fraction=1.1
bin_size={2.5*L:.12e}
every=1
[dem]
contact_model="hertz"
[[dem.materials]]
name="bpm"
youngs_mod={m['youngs_modulus']:.12e}
poisson_ratio=0.25
restitution=1.0
friction=0.0
[[particles.insert]]
source="file"
file="{geo.relative_to(REPO)}"
format="csv"
material="bpm"
radius={m['sphere_radius']:.12e}
density={m['density']:.12e}
columns={{x=0,y=1,z=2}}
[bonds]
auto_bond=true
bond_tolerance=1.001
bond_radius_ratio={m['bond_radius']/m['sphere_radius']:.12e}
youngs_modulus={m['youngs_modulus']:.12e}
shear_modulus={m['shear_modulus']:.12e}
beta_normal=0.0
beta_shear=0.0
beta_twist=0.0
beta_bending=0.0
[[group]]
name="anchor"
region={{type="block",min=[{-0.25*L:.12e},{-L:.12e},{-L:.12e}],max=[{0.25*L:.12e},{L:.12e},{L:.12e}]}}
dynamic=false
[[freeze]]
group="anchor"
[output]
dir="{out.relative_to(REPO)}"
[run]
steps={steps}
thermo={steps}
dt={frac*dt_guo(m):.12e}
''')
    return p,out
def probe(path):
    with path.open(newline='') as f: return [{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
def run(c,mode,frac):
    m=c['material']; limits=lattice_limits(m); p,out=case(m,mode,frac,c['probe']['steps']); shutil.rmtree(out,ignore_errors=True)
    env=os.environ|{'FIBER_TIMESTEP_MODE':mode,'FIBER_TIMESTEP_DRIVE':str(c['probe'][f'{mode}_velocity'])}
    cp=subprocess.run([str(EXE),str(p.relative_to(REPO))],cwd=REPO,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    rows=probe(out/'data/timestep_probe.csv') if cp.returncode==0 and (out/'data/timestep_probe.csv').exists() else []
    finite=bool(rows) and cp.returncode==0 and all(r['finite']==1 and math.isfinite(r['q']) for r in rows)
    # A stable linear pulse remains O(initial amplitude); this generous bound
    # only classifies explosive growth after it has crossed six decades.
    amp=max((abs(r['q']) for r in rows if math.isfinite(r['q'])),default=math.inf)
    stable=finite and amp <= c['acceptance']['max_tip_displacement']
    return {'mode':mode,'dt_over_guo':frac,'dt_over_lattice':frac*dt_guo(m)/limits[mode], 'dt':frac*dt_guo(m),'dt_guo':dt_guo(m),'dt_lattice':limits[mode],'exit_code':cp.returncode,'finite':finite,'max_tip_displacement':amp,'stable':stable}
def bracket(rows,mode,c):
    q=sorted((r for r in rows if r['mode']==mode),key=lambda r:r['dt_over_lattice']); bound=c['acceptance']['transition_bracket']
    for a,b in zip(q,q[1:]):
        if a['stable'] and not b['stable']:
            straddles=a['dt_over_lattice'] < 1.0 < b['dt_over_lattice']
            return straddles and b['dt_over_lattice']-a['dt_over_lattice']<=bound,(a,b)
    return False,None
def output(rows,c):
    DATA.mkdir(exist_ok=True)
    with (DATA/'stability_curve.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    import matplotlib.pyplot as plt
    PLOTS.mkdir(exist_ok=True); fig,ax=plt.subplots(figsize=(8,4.8))
    for mode,color in [('axial','tab:blue'),('bending','tab:orange')]:
        q=[r for r in rows if r['mode']==mode]
        diagnostic=c['acceptance']['max_tip_displacement']
        visible=[max(r['max_tip_displacement'],1e-18) if r['stable'] else diagnostic*10 for r in q]
        ax.semilogy([r['dt_over_lattice'] for r in q],visible,'o-',color=color,label=f'DIRT {mode} tip amplitude')
        for r, y in zip(q, visible):
            if not r['stable']:
                ax.plot(r['dt_over_lattice'],y,'x',ms=9,mew=2,color=color,label=f'{mode} failed trace' if r is q[-1] else None)
        ok,pair=bracket(rows,mode,c)
        if pair: ax.axvspan(pair[0]['dt_over_lattice'],pair[1]['dt_over_lattice'],color=color,alpha=.13,label=f'{mode} measured bracket')
    ax.axvline(1,color='k',ls=':',label='theoretical lattice limit')
    ax.axhline(c['acceptance']['max_tip_displacement'],color='crimson',ls='--',label='explosive-growth diagnostic')
    ax.set(xlabel='dt / dt_lattice',ylabel='maximum free-tip displacement (m)',title='Fixed-free DIRT fiber: empirical stability versus lattice spectrum')
    ax.grid(alpha=.25,which='both');ax.legend(fontsize=7,ncol=2);fig.tight_layout();fig.savefig(PLOTS/'fiber_timestep_stability.png',dpi=160)
def main():
    c=load(); subprocess.run(['cargo','build','--release','--example','bench_fiber_timestep','--no-default-features','--features','precision-double'],cwd=REPO,check=True)
    rows=[run(c,mode,f) for mode in ('axial','bending') for f in c['sweep']['dt_over_guo']]
    output(rows,c); checks=[]
    for mode in ('axial','bending'):
        at=[r for r in rows if r['mode']==mode and r['dt_over_guo']==1]
        br,pair=bracket(rows,mode,c)
        criterion=len(at)==1 and at[0]['stable'] if mode=='axial' else True
        checks += [bool(criterion),br]
        print(f'{mode}: lattice limit={next(r["dt_lattice"] for r in rows if r["mode"]==mode):.6e} s; bounded spectral bracket={br}' + (f' ({pair[0]["dt_over_lattice"]:.4f}, {pair[1]["dt_over_lattice"]:.4f})' if pair else ''))
    axial_limit=next(r['dt_lattice'] for r in rows if r['mode']=='axial')
    guo_agreement=abs(axial_limit/dt_guo(c['material'])-1)<=c['acceptance']['guo_lattice_relative_tolerance']
    checks.append(guo_agreement)
    print(f'Guo Eq. 41 versus axial infinite-lattice limit (tolerance={c["acceptance"]["guo_lattice_relative_tolerance"]:.1%})={guo_agreement}')
    print(f'VALIDATION: {"PASS" if all(checks) else "FAIL"}')
    return 0 if all(checks) else 1
if __name__=='__main__': sys.exit(main())
