#!/usr/bin/env python3
"""
Granular-temperature conductivity from a gas-free vibro-fluidized bed.

A bed of glass beads on an oscillating base wall (periodic x,z, gravity in −y)
reaches a steady fluidized state with profiles Φ(y), T(y). The base injects
fluctuation energy that conducts upward and is dissipated inelastically, so in
steady state (no mean shear) the energy balance is pure conduction vs dissipation:

    d q_y/dy = −Γ(Φ,T),   q_y = −κ dT/dy,   Γ = (12/d√π)(1−e²) ρ_s Φ² g₀ T^{3/2}.

Integrating dissipation from the top down gives the upward energy flux q_y(y), and
κ(y) = q_y(y) / (−dT/dy). Because Φ varies with height, one run sweeps **κ(Φ)**,
which we compare to the KT (Lun/Gidaspow) conductivity. The measured kinetic heat
flux q_y(y) (recorded directly) is overlaid as an independent cross-check.

Commands:
    python3 examples/bench_granular_conductivity/sweep.py start   # build + run
    python3 examples/bench_granular_conductivity/sweep.py graph   # profiles + κ(Φ)
    python3 examples/bench_granular_conductivity/sweep.py         # both

Outputs:
    data/conductivity_profiles.csv     per-y-bin time series      (gitignored)
    plots/profiles.png, kappa_of_phi.png                          (tracked)

To run the **de-fluidization** transient instead, set the base `oscillate.amplitude`
to 0 in config.toml (restart from a fluidized state) and watch T(y,t) decay.
"""

import os
import sys
import csv
import math
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_granular_conductivity"
CONFIG = os.path.join(SCRIPT_DIR, "config.toml")
DATA = os.path.join(SCRIPT_DIR, "data", "conductivity_profiles.csv")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")

# Must match config.toml.
RHO_S = 2500.0
D_MEAN = 0.000225 + 0.000275   # mean diameter [m]
E = 0.926                       # canonical glass-bead COR (matches config restitution)


def g0_cs(phi):
    return (2.0 - phi) / (2.0 * (1.0 - phi) ** 3) if phi < 0.999 else 1e6


def kt_kappa_star(phi, e):
    """Dimensionless KT conductivity κ* = κ/(ρ_s d √T), Lun/Gidaspow. Dilute limit
    κ*→150√π/768 ≈ 0.346 (e→1), giving κ*/η* = 15/4 as kinetic theory requires."""
    g0 = g0_cs(phi)
    kin = (150.0 * math.sqrt(math.pi)) / (384.0 * (1.0 + e) * g0) \
        * (1.0 + (6.0 / 5.0) * (1.0 + e) * g0 * phi) ** 2
    col = 2.0 * phi * phi * g0 * (1.0 + e) / math.sqrt(math.pi)
    return kin + col


def dissipation(phi, T, e):
    """Collisional dissipation rate per unit volume Γ(Φ,T) [W/m³]."""
    return (12.0 / (D_MEAN * math.sqrt(math.pi))) * (1.0 - e * e) * RHO_S \
        * phi * phi * g0_cs(phi) * T ** 1.5


def start():
    subprocess.run(["cargo", "build", "--release", "--no-default-features", "--example", EXAMPLE],
                   cwd=REPO_ROOT, check=True)
    exe = os.path.join(REPO_ROOT, "target", "release", "examples", EXAMPLE)
    subprocess.run([exe, CONFIG], cwd=REPO_ROOT, check=True)


def steady_profiles(frac=0.5):
    """Time-average each y-bin over the steady window (last `frac` of snapshots)."""
    rows = list(csv.DictReader(open(DATA)))
    steps = sorted({int(r["step"]) for r in rows})
    keep = set(steps[int(len(steps) * frac):])
    acc = {}  # y -> [phi_sum, T_sum, qy_sum, n]
    for r in rows:
        if int(r["step"]) not in keep:
            continue
        y = float(r["y"])
        a = acc.setdefault(y, [0.0, 0.0, 0.0, 0])
        a[0] += float(r["phi"]); a[1] += float(r["T"]); a[2] += float(r["qy"]); a[3] += 1
    ys = sorted(acc)
    phi = [acc[y][0] / acc[y][3] for y in ys]
    T = [acc[y][1] / acc[y][3] for y in ys]
    qy = [acc[y][2] / acc[y][3] for y in ys]
    return ys, phi, T, qy


def graph():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(PLOT_DIR, exist_ok=True)
    if not os.path.isfile(DATA):
        print("no data — run `start` first")
        return
    ys, phi, T, qy = steady_profiles()
    dy = ys[1] - ys[0]

    # --- Profiles ---
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(13, 4.2))
    a1.plot(phi, ys, "o-"); a1.set_xlabel("Φ"); a1.set_ylabel("y [m]"); a1.set_title("Solid fraction")
    a2.plot(T, ys, "o-", color="tab:red"); a2.set_xlabel(r"$T$ [m²/s²]"); a2.set_title("Granular temperature"); a2.set_xscale("log")
    a3.plot(qy, ys, "o-", color="tab:green"); a3.set_xlabel(r"$q_y$ (kinetic) [W/m²]"); a3.set_title("Heat flux"); a3.axvline(0, c="k", lw=0.5)
    for a in (a1, a2, a3):
        a.grid(True, alpha=0.3)
    fig.suptitle("Vibro-fluidized bed: steady profiles")
    fig.tight_layout(); fig.savefig(os.path.join(PLOT_DIR, "profiles.png"), dpi=130); plt.close(fig)

    # --- κ(Φ): energy-balance (integrated dissipation) + direct kinetic flux ---
    n = len(ys)
    Gamma = [dissipation(phi[i], T[i], E) for i in range(n)]
    # Upward flux at each bin = total dissipation at/above it (per unit area).
    q_up = [0.0] * n
    run = 0.0
    for i in range(n - 1, -1, -1):
        run += Gamma[i] * dy
        q_up[i] = run
    phis, k_eb, k_kin = [], [], []
    for i in range(1, n - 1):
        dTdy = (T[i + 1] - T[i - 1]) / (2 * dy)
        if dTdy >= 0 or T[i] <= 0 or phi[i] < 0.02 or phi[i] > 0.62:
            continue
        norm = RHO_S * D_MEAN * math.sqrt(T[i])      # κ → κ* = κ/(ρ_s d √T)
        phis.append(phi[i])
        k_eb.append(q_up[i] / (-dTdy) / norm)        # total κ* (energy balance)
        k_kin.append(qy[i] / (-dTdy) / norm)         # kinetic-only κ* (direct, lower bound)

    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    phi_line = [0.02 + 0.005 * k for k in range(0, 120)]
    ax.plot(phi_line, [kt_kappa_star(x, E) for x in phi_line], "k-", lw=2, label=f"KT (Lun/Gidaspow, e={E})")
    if phis:
        ax.scatter(phis, k_eb, c="tab:blue", zorder=3, label="DEM κ* (energy balance, total)")
        ax.scatter(phis, k_kin, facecolors="none", edgecolors="tab:green", zorder=3, label="DEM κ* (kinetic flux only)")
    ax.set_xlabel("Φ"); ax.set_ylabel(r"$\kappa^* = \kappa / (\rho_s d \sqrt{T})$")
    ax.set_yscale("log"); ax.set_title("Granular-temperature conductivity vs Φ"); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(PLOT_DIR, "kappa_of_phi.png"), dpi=130); plt.close(fig)
    print(f"wrote {PLOT_DIR}/profiles.png and kappa_of_phi.png")
    if phis:
        print(f"κ* range (energy-balance): {min(k_eb):.2f}–{max(k_eb):.2f} over Φ {min(phis):.3f}–{max(phis):.3f}")


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    if cmd == "start":
        start()
    elif cmd == "graph":
        graph()
    else:
        start(); graph()


if __name__ == "__main__":
    main()
