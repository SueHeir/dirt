import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# DIRT: step,strain,stress_zz,virial_zz,bond_count,gauge,mean_bond_eps,max_bond_eps
dirt=[]
for l in open("examples/bench_yarn_pull_coarse19/data/stress_strain.csv").read().splitlines()[1:]:
    f=l.split(","); dirt.append((float(f[1])*100, -float(f[2]), int(f[4])))
# LAMMPS: strain P sxx syy szz ... nbond
lam=[]
for l in open("examples/bench_yarn_pull_coarse19/data/lammps_stress_strain.dat").read().splitlines():
    if l.startswith("#") or not l.strip(): continue
    f=l.split(); lam.append((float(f[0])*100, float(f[4]), int(float(f[8]))//2))

N=5311; Vbox=720*720*13000
d_s=[p[0] for p in dirt]; d_stress=[p[1] for p in dirt]; d_bf=[100*p[2]/N for p in dirt]
l_s=[p[0] for p in lam]; l_stress=[p[1]/Vbox for p in lam]; l_bf=[100*p[2]/N for p in lam]

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,4.6))
ax1.plot(l_s,l_stress,'-',color="#1f77b4",lw=2,label="LAMMPS")
ax1.plot(d_s,d_stress,'--',color="#d62728",lw=2,label="DIRT")
ax1.set_xlabel("axial strain (%)"); ax1.set_ylabel(r"axial stress $\sigma_{zz}=-\mathrm{virial}_{zz}/V_{box}$ (lj)")
ax1.set_title("Stress–strain"); ax1.legend(); ax1.grid(alpha=0.3); ax1.set_xlim(0,1.75)

ax2.plot(l_s,l_bf,'-',color="#1f77b4",lw=2,label="LAMMPS")
ax2.plot(d_s,d_bf,'--',color="#d62728",lw=2,label="DIRT")
ax2.set_xlabel("axial strain (%)"); ax2.set_ylabel("bonds intact (%)")
ax2.set_title("Bond survival"); ax2.legend(); ax2.grid(alpha=0.3); ax2.set_xlim(0,1.75)

fig.suptitle("DIRT vs LAMMPS — coarse-19 yarn tensile pull (relaxed, rigid grips)",fontsize=12)
fig.tight_layout(rect=[0,0,1,0.95])
fig.savefig("examples/bench_yarn_pull_coarse19/plots/stress_strain_comparison.png",dpi=140)
# stats
dpk=max(d_stress); dpk_s=d_s[d_stress.index(dpk)]
lpk=max(l_stress); lpk_s=l_s[l_stress.index(lpk)]
print(f"DIRT   peak σ={dpk:.4f} at {dpk_s:.2f}% strain")
print(f"LAMMPS peak σ={lpk:.4f} at {lpk_s:.2f}% strain")
print("wrote plots/stress_strain_comparison.png")
