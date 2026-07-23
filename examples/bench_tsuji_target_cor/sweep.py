#!/usr/bin/env python3
"""Calibrate physical normal COR to DIRT's legacy Hertz--Tsuji input and validate it.

The inversion deliberately integrates the contact ODE here, independently of the
DIRT contact kernel.  The generated TOML keeps `restitution` as the legacy raw
Tsuji value; no existing input semantics are changed.
"""
from __future__ import annotations
import csv, math, os, shutil, subprocess, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA, PLOTS, CASES = HERE / "data", HERE / "plots", HERE / "sweep"
TARGETS, VELOCITIES, RADII, DENSITIES, DT_FRACS = [0.50, 0.70, 0.90], [0.25, 1.0], [0.0025, 0.005], [1000.0, 2500.0], [0.05, 0.15]
E, NU = 70e9, 0.22
TARGET_TOL, ODE_TOL, PARITY_TOL = 0.015, 5e-4, 0.015
LMP = next((shutil.which(x) for x in ("lmp_serial", "lmp", "lmp_mpi", "lammps") if shutil.which(x)), None)
BENCH_BIN = os.environ.get("DIRT_REBOUND_BIN", str(ROOT / "target/release/examples/bench_hertz_rebound"))
CALIBRATOR_BIN = os.environ.get("DIRT_TSUJI_CALIBRATOR", str(ROOT / "target/release/examples/bench_tsuji_target_cor"))

TOML = '''[comm]\nprocessors_x=1\nprocessors_y=1\nprocessors_z=1\n[domain]\nx_low=-0.01\nx_high=0.01\ny_low=-0.01\ny_high=0.01\nz_low=0.0\nz_high=0.1\nboundary_x="fixed"\nboundary_y="fixed"\nboundary_z="fixed"\n[neighbor]\nskin_fraction=1.1\nbin_size=0.015\nevery=1\n[dem]\ncontact_model="hertz"\nlimit_damping=false\n[[dem.materials]]\nname="glass"\nyoungs_mod={E}\npoisson_ratio={nu}\nrestitution={raw}\nfriction=0.0\n[[particles.insert]]\nmaterial="glass"\ncount=1\nradius={r}\ndensity={rho}\nvelocity_z=-{v}\nregion={{ type="block", min=[-0.001,-0.001,{z}], max=[0.001,0.001,{z2}] }}\n[[wall]]\npoint_x=0.0\npoint_y=0.0\npoint_z=0.0\nnormal_x=0.0\nnormal_y=0.0\nnormal_z=1.0\nmaterial="glass"\n[output]\ndir="{out}"\n[run]\nsteps={steps}\ndt={dt}\nthermo=1000\n'''
LMP_IN = '''units si\natom_style sphere\natom_modify map array\ndimension 3\nboundary f f f\nnewton off\ncomm_modify vel yes\nregion box block -0.01 0.01 -0.01 0.01 0 0.1 units box\ncreate_box 1 box\ncreate_atoms 1 single 0 0 {z} units box\nset group all diameter {diam}\nset group all density {rho}\npair_style granular\npair_coeff 1 1 hertz/material {E} {raw} {nu} tangential linear_nohistory 0 0 damping tsuji rolling none twisting none\nfix wall all wall/gran granular hertz/material {E} {raw} {nu} tangential linear_nohistory 0 0 damping tsuji rolling none twisting none zplane 0 NULL\nfix nve all nve/sphere\nvelocity all set 0 0 -{v} units box\ntimestep {dt}\nvariable zpos equal z[1]\nvariable zvel equal vz[1]\nfix rec all print 1 "${{zpos}} ${{zvel}}" file {trace} screen no title "z vz"\nrun {steps}\n'''

def alpha(e): return 1.2728-4.2783*e+11.087*e**2-22.348*e**3+27.467*e**4-18.022*e**5+4.8218*e**6
def beta(raw): return 0.0 if raw >= .9999 else alpha(max(.001, min(.9999, raw)))/math.sqrt(5)
def ode_cor(raw):
    """Dimensionless RK4 Hertz--Tsuji collision, independent reference."""
    c = 2*beta(raw)*math.sqrt(5/6)*math.sqrt(2); d,v,dt = 0.,1.,1e-4
    def a(x,y): return 0. if x <= 0 else -(4/3)*x**1.5-c*x**.25*y
    for _ in range(2000000):
        k1d,k1v=v,a(d,v); k2d,k2v=v+.5*dt*k1v,a(d+.5*dt*k1d,v+.5*dt*k1v)
        k3d,k3v=v+.5*dt*k2v,a(d+.5*dt*k2d,v+.5*dt*k2v); k4d,k4v=v+dt*k3v,a(d+dt*k3d,v+dt*k3v)
        dn=d+dt*(k1d+2*k2d+2*k3d+k4d)/6; vn=v+dt*(k1v+2*k2v+2*k3v+k4v)/6
        if d > 0 and dn <= 0: return abs(v+d/(d-dn)*(vn-v))
        d,v=dn,vn
    raise RuntimeError("ODE did not separate")
def raw_for_target(target):
    lo,hi=.001,.9999
    for _ in range(48):
        mid=(lo+hi)/2
        if ode_cor(mid) > target: hi=mid
        else: lo=mid
    return (lo+hi)/2
def dirt_raw_for_target(target):
    """Invoke DIRT's public conversion; the Python ODE remains the oracle."""
    text = subprocess.check_output([CALIBRATOR_BIN, "calibrate", f"{target:.17g}"], text=True).strip()
    got_target, raw, _beta = map(float, text.split(','))
    if abs(got_target-target) > 1e-12: raise ValueError("calibrator returned wrong target")
    return raw
def dt_for(r,rho,frac):
    g=E/(2*(1+NU)); return frac*math.pi*r/(.1631*NU+.876605)*(rho/g)**.5
def key(t,v,r,rho,f): return f"t{t:.2f}_v{v:g}_r{r*1e3:g}mm_rho{rho:g}_dt{f:g}"
def parse(path, radius):
    row=next(csv.DictReader(open(path))); return float(row["cor_measured"])
def lmp_cor(trace,radius):
    hit=False; prev=None
    for line in open(trace):
        p=line.split()
        # `fix print` emits exactly the two declared observables, z and vz.
        # Treat any other non-header record as malformed rather than quietly
        # dropping the independent-code comparison.
        if p == ["z", "vz"]: continue
        if len(p) != 2: raise ValueError(f"malformed LAMMPS rebound trace: {line.rstrip()}")
        z,v=map(float,p); contact=radius-z>0
        if not hit and not contact: prev=v
        elif not hit and contact: hit=True; vin=prev
        elif hit and not contact: return abs(v/vin)
    return None
def cases():
    for t in TARGETS:
      for v in VELOCITIES:
       for r in RADII:
        for rho in DENSITIES:
         for frac in DT_FRACS: yield t,v,r,rho,frac
def run():
    DATA.mkdir(exist_ok=True); CASES.mkdir(exist_ok=True)
    if not Path(BENCH_BIN).is_file() or not Path(CALIBRATOR_BIN).is_file():
      subprocess.run(["cargo","build","--release","--example","bench_hertz_rebound","--example","bench_tsuji_target_cor","--no-default-features","--features","precision-double"],cwd=ROOT,check=True)
    if not Path(BENCH_BIN).is_file() or not Path(CALIBRATOR_BIN).is_file():
      raise RuntimeError("rebound or target-COR calibration executable was not produced")
    rows=[]
    for t,v,r,rho,frac in cases():
      raw=dirt_raw_for_target(t)
      name=key(t,v,r,rho,frac); d=CASES/name; d.mkdir(parents=True,exist_ok=True); dt=dt_for(r,rho,frac); steps=int((5e-5/v+5e-5)/dt)+1000
      # A 10 µm gap avoids a long free-flight loop while preserving a clean impact.
      z=r+1e-5
      (d/"config.toml").write_text(TOML.format(E=E,nu=NU,raw=raw,r=r,rho=rho,v=v,z=z,z2=z+.001,out=d,steps=steps,dt=dt))
      subprocess.run([BENCH_BIN,str(d/"config.toml")],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
      dirt=parse(d/"data/rebound_results.csv",r); lcor=""
      if LMP:
       trace=d/"trace.lmp"; (d/"in.lammps").write_text(LMP_IN.format(E=E,nu=NU,raw=raw,rho=rho,v=v,z=z,diam=2*r,dt=dt,trace=trace,steps=steps))
       subprocess.run([LMP,"-in",str(d/"in.lammps"),"-log",str(d/"lammps.log")],cwd=ROOT,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
       lcor=lmp_cor(trace,r)
      rows.append(dict(target=t,raw=raw,ode=ode_cor(raw),velocity=v,radius=r,density=rho,dt_fraction=frac,dirt=dirt,lammps=lcor if lcor is not None else ""))
      print(name, f"raw={raw:.6f}", f"DIRT={dirt:.5f}")
    with open(DATA/"results.csv","w",newline="") as f:
      writer=csv.DictWriter(f,fieldnames=rows[0],lineterminator="\n")
      writer.writeheader(); writer.writerows(rows)
    return rows
def graph(rows):
    PLOTS.mkdir(exist_ok=True)
    targets=[float(r['target']) for r in rows]; dirt=[float(r['dirt']) for r in rows]; lammps=[float(r['lammps']) for r in rows]
    fig, ax = plt.subplots(figsize=(7.4, 5.2), constrained_layout=True)
    x=[0.45, 0.95]
    ax.fill_between(x, [v-TARGET_TOL for v in x], [v+TARGET_TOL for v in x], color='#dcefd8', label='DIRT target gate ±0.015')
    ax.plot(x, x, 'k--', label='declared target'); ax.scatter(targets, dirt, c='#195ac8', label='DIRT', zorder=3)
    ax.scatter(targets, lammps, facecolors='none', edgecolors='black', marker='s', label='LAMMPS', zorder=3)
    ax.set(xlim=(.45,.95), ylim=(.45,.95), xlabel='Declared physical target COR', ylabel='Realized rebound COR', title='Hertz–Tsuji target-COR calibration (48 cases)')
    ax.set_xticks([.5,.6,.7,.8,.9]); ax.set_yticks([.5,.6,.7,.8,.9]); ax.grid(alpha=.25); ax.legend(loc='upper left', fontsize=8)
    fig.savefig(PLOTS/'target_cor.png', dpi=160); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7.4, 5.2), constrained_layout=True); index=list(range(1,len(rows)+1))
    ax.axhspan(-TARGET_TOL,TARGET_TOL,color='#dcefd8',label='target gate ±0.015'); ax.axhline(0,color='black',linewidth=.8)
    ax.scatter(index,[d-t for d,t in zip(dirt,targets)],c='#195ac8',label='DIRT',zorder=3)
    ax.scatter(index,[l-t for l,t in zip(lammps,targets)],facecolors='none',edgecolors='black',marker='s',label='LAMMPS',zorder=3)
    ax.set(xlabel='Sweep case (target × velocity × radius × density × timestep)', ylabel='Realized − target COR', title='Target-COR error across the 48-case campaign', xlim=(0,49), ylim=(-.02,.02))
    ax.set_xticks([1,16,17,32,33,48], ['0.50', '', '0.70', '', '0.90', '']); ax.grid(alpha=.25); ax.legend(loc='upper left', fontsize=8)
    fig.savefig(PLOTS/'calibration_error.png', dpi=160); plt.close(fig)
def validate(rows):
    failures=[]
    for r in rows:
      if abs(float(r['ode'])-float(r['target']))>ODE_TOL: failures.append('independent ODE inversion')
      if abs(float(r['dirt'])-float(r['target']))>TARGET_TOL: failures.append('DIRT target')
      if LMP and r['lammps']=='': failures.append('missing LAMMPS rebound')
      if r['lammps']!='' and abs(float(r['dirt'])-float(r['lammps']))>PARITY_TOL: failures.append('LAMMPS parity')
    print(f"Calibration: {len(rows)} cases; LAMMPS={'yes' if LMP else 'not available'}")
    print("ALL CHECKS PASSED" if not failures else f"CHECKS FAILED: {len(failures)}")
    return not failures
def main():
    if len(sys.argv)>1 and sys.argv[1]=='calibrate':
      t=float(sys.argv[2]); raw=raw_for_target(t); print(f"target={t:.6f} raw_tsuji_input={raw:.9f} beta={beta(raw):.9f} ode_cor={ode_cor(raw):.9f}"); return
    rows=run(); graph(rows); sys.exit(0 if validate(rows) else 1)
if __name__=='__main__': main()
