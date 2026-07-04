#!/usr/bin/env python3
"""gen_series — generate the monkey-barrel LEBC volume-fraction shear campaign.

Emits a 3 (particle type) × 8 (volume fraction Φ) grid of DECLARATIVE TOML
configs for the `monkey_shear` binary, all at a COMMON equivalent-volume
diameter D_eq = 0.1 m so the three shapes displace the *identical* solid volume
per particle (V_eq = π·D_eq³/6 = 5.236e-4 m³) and share one N(Φ) table.

Three particle types (all D_eq = 0.1):
  * sphere : a single sphere, radius 0.05 m.
  * rigid  : the 44-sub-sphere "monkey" clump (multisphere rigid body).
  * bpm    : the same monkey as a bonded-particle model — 44 free sub-spheres
             per monkey, welded INTRA-monkey by `auto_bond` at setup. Sub-sphere
             centres are emitted as a CSV; monkeys are placed with a gap so no
             two monkeys' sub-spheres are ever close enough to auto-bond.

Protocol per (type, Φ), 3 stages (cf. examples/bench_lebc_shear/config.toml):
  1. settle  — insert a LOOSE pack in an ENLARGED, aspect-preserving (2:2:1) box,
               relax with no deform (mild prep-damping cools it — see main.rs).
  2. compress — quasi-static isotropic `erate` box shrink to the final 2×2×1 box
               (Φ target), atoms + clump COMs remapped. SKIPPED when s == 1.
  3. shear   — Lees–Edwards `xy erate` at γ̇ (constant volume), strain ≈ 2.

Subcommands:
  build-monkey   scale the raw monkey to D_eq=0.1, MC-verify, write monkey_Deq0.1.toml
  series         emit the full 3×8 config grid (+ bpm CSVs)   [--gdot, --smoke]
  all (default)  build-monkey then series

Run (from the dirt repo root):
  python3 examples/monkey_shear/tools/gen_series.py all
  python3 examples/monkey_shear/tools/gen_series.py series --smoke   # phi=0.1 only, short
"""
from __future__ import annotations
import argparse
import math
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EX_DIR = os.path.dirname(HERE)                 # examples/monkey_shear
MONKEY_TOML = os.path.join(EX_DIR, "monkey_Deq0.1.toml")

# ── Campaign constants ───────────────────────────────────────────────────────
D_EQ = 0.1                                     # common equivalent-volume diameter [m]
V_EQ = math.pi * D_EQ**3 / 6.0                 # 5.235988e-4 m³
DOMAIN_VOL = 4.0                               # final box 2×2×1
HALF_FINAL = (1.0, 1.0, 0.5)                   # half-extents of the final box
CENTER = (1.0, 1.0, 0.5)                       # box centre (invariant under erate/vel deform)

# Reference glass (softened rigid-grain limit), identical to bench_lebc_shear.
GLASS = dict(youngs_mod=7.0e7, poisson_ratio=0.245, restitution=0.926,
             friction=0.16, density=2500.0)
E_BOND = 7.0e7
G_BOND = E_BOND / (2.0 * (1.0 + GLASS["poisson_ratio"]))   # 2.811e7 Pa

# N per Φ — IDENTICAL across types, N = round(Φ · DOMAIN_VOL / V_EQ).
PHIS = [0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55]
N_TABLE = {0.025: 191, 0.05: 382, 0.1: 764, 0.2: 1528,
           0.3: 2292, 0.4: 3056, 0.5: 3820, 0.55: 4202}

CAP_SPHERE = 0.25          # loose-insertion solid fraction cap for spheres
# Monkeys insert looser than spheres (the task allows <= 0.15). The bpm gap rule
# (reject a monkey if any sub-sphere is within GAP*(Ri+Rj) of another monkey's)
# excludes far more than the union volume for these elongated, randomly-oriented
# bodies, so union fraction 0.15 is unreachable; 0.10 places reliably. The gap
# only has to hold at SETUP (auto_bond runs once); compress then closes the gaps.
CAP_MONKEY = 0.10          # loose-insertion solid fraction cap for monkeys

# ── Stable-timestep (CFL / Rayleigh) model ───────────────────────────────────
# The SHEAR stage runs fully UNDAMPED (see main.rs: damping acts only during the
# prep stages), so any spurious energy an over-large dt injects at a stiff
# collision has nowhere to go — it accumulates, heats the pack, and finally
# trips the excessive-overlap guard (crates/dirt_granular/src/contact.rs). The
# cure is to make dt small enough to integrate the STIFFEST interaction present
# in each run, NOT to relax the guard. We size dt as a safe fraction of the
# critical timestep of that stiffest element, per particle type, per Φ:
#
#   * sphere : the stiffest element is the Hertz contact between the (large,
#              r = 0.05 m) grains. Its Rayleigh time is long, so the reference
#              dt below is already deep in the stable regime — LEFT UNCHANGED
#              (the whole sphere series validated at this dt).
#   * rigid / bpm : the 44-sub-sphere monkey. Two stiff scales compete — the
#              Hertz contact of the smallest sub-sphere (r_min ≈ 8.5 mm) and,
#              for bpm, the E_bond WELD network, which is ~8× stiffer and hence
#              sets the smaller (governing) critical step. We compute the
#              bond-network critical step explicitly from the engine's own bond
#              stiffnesses so E_bond enters the timestep, and use that same
#              conservative scale for the rigid monkey (identical geometry).
#
# Denser packings add a high-Φ margin: a grain sharing n_c ≈ 1 + Z·Φ/Φ_rcp
# simultaneous contacts feels ~n_c× the stiffness, so its stable dt falls
# ~1/√n_c. This is why the flat legacy dt survived the dilute cases but blew up
# the dense high-Φ compressions.
DT_SPHERE = 2.0e-5         # r=0.05 m contact-Rayleigh margin — validated, UNCHANGED.

G_GLASS = GLASS["youngs_mod"] / (2.0 * (1.0 + GLASS["poisson_ratio"]))  # shear modulus [Pa]
PHI_RCP = 0.64             # random-close-pack fraction (coordination ~saturates here)
COORD_Z = 6.0              # contacts added per grain by an isostatic pack at Φ_rcp
DT_SHEAR_SAFETY = 0.02     # fraction of the stiffest-element critical step (settle & shear)
DT_COMPRESS_QS = 2.5       # quasi-static compression tolerates a modestly larger dt
EDOT_COMPRESS = 3.0        # quasi-static isotropic compression strain rate [1/s]

BOND_TOL = 1.1             # auto_bond eligibility multiplier on (R_i+R_j) (welds intra-monkey)
GAP = 1.25                 # inter-monkey rejection factor (> BOND_TOL ⇒ no inter-monkey bonds)

# Raw 44-sphere monkey (examples/monkey_barrel/settle_cad.toml, longest axis 6 mm).
RAW_SPHERES = [
    (+0.000312, -0.000062, +0.001003, 0.000640),
    (+0.000312, -0.000062, -0.001961, 0.000591),
    (+0.000312, -0.000062, -0.001302, 0.000591),
    (+0.000312, -0.000062, -0.000644, 0.000591),
    (+0.000312, -0.000062, +0.000015, 0.000591),
    (-0.001444, +0.000047, +0.001113, 0.000439),
    (-0.001335, -0.000062, +0.000344, 0.000439),
    (-0.000457, +0.000047, +0.001552, 0.000439),
    (-0.000237, +0.000047, +0.000344, 0.000439),
    (+0.000861, +0.000047, +0.000344, 0.000439),
    (+0.000970, +0.000047, +0.001552, 0.000439),
    (+0.001958, +0.000047, +0.001222, 0.000439),
    (+0.002068, +0.000047, +0.000674, 0.000439),
    (-0.001115, -0.000062, +0.001442, 0.000411),
    (-0.000457, +0.000047, -0.002510, 0.000411),
    (+0.000861, -0.000062, -0.002180, 0.000411),
    (+0.001080, +0.000047, -0.002619, 0.000411),
    (+0.001519, -0.000062, +0.001552, 0.000411),
    (-0.000896, -0.000062, +0.000234, 0.000396),
    (-0.000347, +0.000047, -0.002071, 0.000396),
    (+0.001300, -0.000062, +0.000234, 0.000396),
    (+0.001739, -0.000062, +0.000234, 0.000396),
    (-0.000676, -0.000062, -0.001741, 0.000329),
    (+0.000861, -0.000062, -0.001741, 0.000329),
    (-0.001664, -0.000062, +0.000674, 0.000310),
    (-0.000896, +0.000047, -0.001302, 0.000310),
    (-0.000237, +0.000047, -0.001632, 0.000310),
    (-0.000017, +0.000157, +0.001552, 0.000310),
    (+0.000312, -0.000062, +0.001661, 0.000310),
    (+0.000312, +0.000377, +0.000454, 0.000310),
    (-0.000017, -0.000172, +0.001552, 0.000269),
    (+0.000422, +0.000267, +0.001552, 0.000269),
    (-0.000896, +0.000157, +0.001771, 0.000245),
    (-0.000786, -0.000172, +0.001771, 0.000245),
    (+0.002068, +0.000157, +0.000234, 0.000245),
    (-0.001554, +0.000267, +0.000564, 0.000220),
    (-0.001444, -0.000282, +0.000783, 0.000220),
    (-0.000676, -0.000062, +0.000564, 0.000220),
    (-0.000676, +0.000267, +0.000234, 0.000220),
    (-0.000566, -0.000282, +0.000234, 0.000220),
    (-0.000566, -0.000062, -0.002949, 0.000220),
    (-0.000566, -0.000062, +0.000015, 0.000220),
    (-0.000237, -0.000282, -0.001851, 0.000220),
    (-0.000127, +0.000157, -0.000973, 0.000220),
]


# ── Geometry helpers ─────────────────────────────────────────────────────────
def mc_union_deq(centers, radii, n_samples=4_000_000, seed=1):
    """Monte-Carlo union volume of a sphere set → equivalent-volume diameter.
    Returns (D_eq, V_union, sum_of_sphere_volumes)."""
    centers = np.asarray(centers, float)
    radii = np.asarray(radii, float)
    lo = (centers - radii[:, None]).min(axis=0)
    hi = (centers + radii[:, None]).max(axis=0)
    box_vol = float(np.prod(hi - lo))
    rng = np.random.default_rng(seed)
    inside = 0
    batch = 500_000
    done = 0
    r2 = radii**2
    while done < n_samples:
        m = min(batch, n_samples - done)
        pts = lo + (hi - lo) * rng.random((m, 3))
        hit = np.zeros(m, dtype=bool)
        for k in range(0, len(centers), 8):
            c = centers[k:k+8]
            rr = r2[k:k+8]
            d2 = ((pts[:, None, :] - c[None, :, :])**2).sum(axis=2)
            hit |= (d2 <= rr[None, :]).any(axis=1)
        inside += int(hit.sum())
        done += m
    v_union = box_vol * inside / n_samples
    d_eq = (6.0 * v_union / math.pi)**(1.0 / 3.0)
    v_sum = float((4.0 / 3.0 * math.pi * radii**3).sum())
    return d_eq, v_union, v_sum


def scaled_monkey():
    """Recentre + scale the raw monkey so its union D_eq = 0.1 m. Returns
    (centers[N,3], radii[N], report_dict)."""
    arr = np.array(RAW_SPHERES, float)
    centers = arr[:, :3].copy()
    radii = arr[:, 3].copy()
    centers -= centers.mean(axis=0)                 # recentre on sub-sphere centroid
    d_eq_raw, v_u_raw, _ = mc_union_deq(centers, radii, seed=7)
    scale = D_EQ / d_eq_raw
    centers *= scale
    radii *= scale
    d_eq, v_u, v_sum = mc_union_deq(centers, radii, seed=11)
    r_bound = float((np.linalg.norm(centers, axis=1) + radii).max())
    n_bonds = 0
    for i in range(len(centers)):
        for j in range(i + 1, len(centers)):
            d = float(np.linalg.norm(centers[i] - centers[j]))
            if d <= BOND_TOL * (radii[i] + radii[j]):
                n_bonds += 1
    report = dict(scale=scale, d_eq_raw=d_eq_raw, d_eq=d_eq, v_union=v_u,
                  v_sum=v_sum, v_eq=V_EQ, r_min=float(radii.min()),
                  r_max=float(radii.max()), r_bound=r_bound, n_bonds=n_bonds,
                  overlap=1.0 - v_u / v_sum)
    return centers, radii, report


def random_rotation(rng):
    """Uniform random 3×3 rotation matrix (QR of a Gaussian matrix, sign-fixed)."""
    a = rng.standard_normal((3, 3))
    q, r = np.linalg.qr(a)
    q *= np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    return q


# ── Stable-timestep helpers ──────────────────────────────────────────────────
def t_rayleigh(radius):
    """Rayleigh surface-wave transit time of a grain — the standard DEM contact
    critical timestep: T_R = π R √(ρ/G) / (0.1631 ν + 0.8766)."""
    return (math.pi * radius * math.sqrt(GLASS["density"] / G_GLASS)
            / (0.1631 * GLASS["poisson_ratio"] + 0.8766))


def bond_network_critical_dt(centers, radii):
    """Highest-frequency stable timestep of the welded intra-monkey bond network.

    Mirrors the engine's per-bond stiffnesses (crates/dirt_bond/src/lib.rs):
        r_b = bond_radius_ratio·min(R_i,R_j)   (bond_radius_ratio = 1.0 here)
        A = π r_b²,  J = ½π r_b⁴,  I = ½J,  L = r₀ (rest length)
        k_n = E_b·A/L,  k_t = G_b·A/L,  k_tor = G_b·J/L,  k_bend = E_b·I/L
    A node bonded to several neighbours sees the SUM of their stiffnesses, so the
    governing frequency is a network property, not a single-bond one. We take a
    Gershgorin bound on the mass- and inertia-normalised stiffness assembly,
        ω_max² ≤ 2·max_i( Σ_j k_ij / m_i ),
    with bending/shear coupled into translation through the L/2 lever arm, and
    return the explicit-integration stability limit dt_crit = 2/ω_max. This is
    where E_bond enters the timestep."""
    n = len(centers)
    rho = GLASS["density"]
    mass = 4.0 / 3.0 * math.pi * rho * radii**3
    inertia = 0.4 * mass * radii**2                 # solid sphere I = 2/5 m R²
    kt_sum = np.zeros(n)                            # translational stiffness per node
    kr_sum = np.zeros(n)                            # rotational stiffness per node
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            d = float(np.linalg.norm(centers[i] - centers[j]))
            if d <= BOND_TOL * (radii[i] + radii[j]):
                r_b = min(radii[i], radii[j])       # bond_radius_ratio = 1.0
                area = math.pi * r_b * r_b
                jpol = 0.5 * math.pi * r_b**4
                iben = 0.5 * jpol
                length = d
                k_n = E_BOND * area / length
                k_t = G_BOND * area / length
                k_tor = G_BOND * jpol / length
                k_bend = E_BOND * iben / length
                lever = 0.5 * length
                kt_sum[i] += k_n + k_t + k_bend / (lever * lever)
                kr_sum[i] += k_tor + k_bend + k_t * lever * lever
    wmax2 = max(float((2.0 * kt_sum / mass).max()),
                float((2.0 * kr_sum / inertia).max()))
    return 2.0 / math.sqrt(wmax2)


def phi_margin(phi):
    """High-Φ multi-contact stiffening margin (≤1): a grain sharing
    n_c ≈ 1 + Z·Φ/Φ_rcp simultaneous contacts oscillates ~√n_c faster, so the
    stable dt shrinks ~1/√n_c."""
    return 1.0 / math.sqrt(1.0 + COORD_Z * phi / PHI_RCP)


def monkey_shear_dt(phi, t_bond, r_min):
    """Stable settle/shear dt for a monkey (rigid or bpm) run at solid fraction Φ.
    The stiffest element is the smaller of the smallest sub-sphere's contact
    Rayleigh time and the bond-network critical step (the latter governs for the
    glass E_bond used here); the rigid monkey reuses this conservative scale."""
    t_crit = min(t_bond, t_rayleigh(r_min))
    return DT_SHEAR_SAFETY * t_crit * phi_margin(phi)


# ── TOML emission ────────────────────────────────────────────────────────────
def clump_def_block(centers, radii, name="monkey"):
    lines = ["[[clump.definitions]]", f'name = "{name}"', "spheres = ["]
    for c, r in zip(centers, radii):
        lines.append(f"    {{ offset = [{c[0]:+.6f}, {c[1]:+.6f}, {c[2]:+.6f}], radius = {r:.6f} }},")
    lines.append("]")
    return "\n".join(lines)


def header_common(kind, phi, n, s, dt, gdot, extra=""):
    return (f"# monkey_shear — {kind} | Φ = {phi} | N = {n} | D_eq = {D_EQ}\n"
            f"# Constant-volume Lees-Edwards simple shear. Final box 2x2x1 (vol {DOMAIN_VOL}),\n"
            f"# fully periodic: flow=x, gradient=y, vorticity=z. Insertion box scaled x{s:.4f}\n"
            f"# (aspect 2:2:1 preserved), then erate-compressed to the final box (skipped if s=1).\n"
            f"# dt={dt:.2e}s, gdot={gdot} s^-1, shear strain target ~2. {extra}\n"
            f"# GENERATED by tools/gen_series.py — declarative config, do not hand-edit.\n")


def domain_block(s):
    hx, hy, hz = HALF_FINAL[0]*s, HALF_FINAL[1]*s, HALF_FINAL[2]*s
    cx, cy, cz = CENTER
    return (hx, hy, hz,
            "[domain]\n"
            f"x_low = {cx-hx:.6f}\nx_high = {cx+hx:.6f}\n"
            f"y_low = {cy-hy:.6f}\ny_high = {cy+hy:.6f}\n"
            f"z_low = {cz-hz:.6f}\nz_high = {cz+hz:.6f}\n"
            'boundary_x = "periodic"\nboundary_y = "periodic"\nboundary_z = "periodic"\n')


def material_block():
    return ("[dem]\ncontact_model = \"hertz\"\n\n"
            "[[dem.materials]]\n"
            'name = "glass"\n'
            f"youngs_mod = {GLASS['youngs_mod']:.1e}\n"
            f"poisson_ratio = {GLASS['poisson_ratio']}\n"
            f"restitution = {GLASS['restitution']}\n"
            f"friction = {GLASS['friction']}\n")


def stage_dts(kind, phi, t_bond, r_min):
    """(settle_dt, compress_dt, shear_dt) for a particle type at solid fraction Φ.

    Sphere keeps its validated flat dt. Monkeys (rigid/bpm) get a Φ-dependent,
    bond-network-limited dt (see monkey_shear_dt); settle and shear share it, and
    the quasi-static compression may run at DT_COMPRESS_QS× that value."""
    if kind == "sphere":
        return DT_SPHERE, DT_SPHERE, DT_SPHERE
    dt = monkey_shear_dt(phi, t_bond, r_min)
    return dt, DT_COMPRESS_QS * dt, dt


def runs_block(s, kind, gdot, smoke, phi, t_bond, r_min):
    """settle → (quasi-static compress if s>1) → constant-volume shear.

    Compression runs at a fixed strain rate EDOT_COMPRESS (quasi-static, so dense
    tangling monkeys can rearrange instead of overlapping); it may use a modestly
    larger dt than shear (DT_COMPRESS_QS×), but still scaled by the same Φ margin
    — the dense end of compression is exactly where the flat legacy dt blew up.
    Shear targets total strain ≈ 2 (full) or a fixed short window (smoke)."""
    dt_settle, dt_compress, dt_shear = stage_dts(kind, phi, t_bond, r_min)
    thermo = 2000
    settle_n = 5000 if smoke else 20000
    shear_n = 20000 if smoke else int(round(2.0 / (gdot * dt_shear)))
    out = []
    out.append("[[run]]\n"
               f'name = "settle"\nsteps = {settle_n}\nthermo = {thermo}\ndt = {dt_settle:.3e}\n')
    if s > 1.0 + 1e-9:
        compress_n = int(math.ceil((1.0 - 1.0 / s) / (dt_compress * EDOT_COMPRESS)))
        rate = (1.0 / s - 1.0) / (dt_compress * compress_n)   # ≈ -EDOT_COMPRESS
        out.append("[[run]]\n"
                   f'name = "compress"\nsteps = {compress_n}\nthermo = {thermo}\ndt = {dt_compress:.3e}\n'
                   f"deform = {{ x = {{ style = \"erate\", rate = {rate:.6e} }}, "
                   f"y = {{ style = \"erate\", rate = {rate:.6e} }}, "
                   f"z = {{ style = \"erate\", rate = {rate:.6e} }} }}\n")
    out.append("[[run]]\n"
               f'name = "shear"\nsteps = {shear_n}\nthermo = {thermo}\ndt = {dt_shear:.3e}\n'
               f"deform = {{ xy = {{ style = \"erate\", rate = {gdot} }} }}\n")
    return "\n".join(out)


def common_tail(out_dir):
    return ("[comm]\nprocessors_x = 1\nprocessors_y = 1\nprocessors_z = 1\n\n"
            "[gravity]\ngx = 0.0\ngy = 0.0\ngz = 0.0\n\n"
            f"[output]\ndir = \"{out_dir}\"\n\n"
            "[vtp]\ninterval = 100000000\n")


def bin_size_for(kind, radii):
    if kind == "sphere":
        return 1.2 * 2 * 0.05
    return 1.2 * 2 * float(np.max(radii))           # >= largest sub-sphere diameter


def neighbor_block(kind, radii):
    return (f"[neighbor]\nskin_fraction = 1.2\nbin_size = {bin_size_for(kind, radii):.5f}\nevery = 1\n")


# ── Per-type config writers ──────────────────────────────────────────────────
def write_sphere(phi, n, gdot, smoke, root):
    s = (phi / min(phi, CAP_SPHERE))**(1.0 / 3.0)
    hx, hy, hz, dom = domain_block(s)
    out_dir = f"examples/monkey_shear/data/sphere/phi_{phi}"
    body = "\n".join([
        header_common("sphere (radius 0.05)", phi, n, s, DT_SPHERE, gdot),
        dom, neighbor_block("sphere", np.array([0.05])), material_block(),
        "[[particles.insert]]\n"
        'material = "glass"\n'
        f"count = {n}\nradius = 0.05\ndensity = {GLASS['density']}\n"
        "velocity = 0.0\n"
        f"region = {{ type = \"block\", min = [{CENTER[0]-hx:.6f}, {CENTER[1]-hy:.6f}, {CENTER[2]-hz:.6f}], "
        f"max = [{CENTER[0]+hx:.6f}, {CENTER[1]+hy:.6f}, {CENTER[2]+hz:.6f}] }}\n"
        f"seed = {12345 + n}\n",
        runs_block(s, "sphere", gdot, smoke, phi, None, 0.05),
        common_tail(out_dir),
    ])
    _write(root, "sphere", phi, body)


def write_rigid(phi, n, gdot, smoke, root, centers, radii, rbound, t_bond):
    dt = monkey_shear_dt(phi, t_bond, float(radii.min()))
    # The engine's `[[clump.insert]]` rejects overlaps by each monkey's BOUNDING
    # sphere (radius rbound = 0.126, i.e. the elongated arm reach), NOT the union
    # D_eq. So the insertion box must be sized for loose bounding-sphere packing
    # (F_BOUND) — sizing it by union fraction leaves it far too small and the
    # placement retry loop stalls. The pack is then compressed to the union Φ.
    F_BOUND = 0.22
    v_bound = (4.0 / 3.0) * math.pi * rbound**3
    v_ins = n * v_bound / F_BOUND
    s = max((v_ins / DOMAIN_VOL)**(1.0 / 3.0), (phi / min(phi, CAP_MONKEY))**(1.0 / 3.0))
    hx, hy, hz, dom = domain_block(s)
    # inset the COM insertion region by the monkey bounding radius
    ix = max(hx - rbound, 1e-3); iy = max(hy - rbound, 1e-3); iz = max(hz - rbound, 1e-4)
    out_dir = f"examples/monkey_shear/data/rigid/phi_{phi}"
    body = "\n".join([
        header_common("rigid clump (44-sphere monkey)", phi, n, s, dt, gdot),
        dom, neighbor_block("rigid", radii), material_block(),
        clump_def_block(centers, radii),
        "\n[[clump.insert]]\n"
        'definition = "monkey"\n'
        f"count = {n}\ndensity = {GLASS['density']}\n"
        'material = "glass"\n'
        "random_orientation = true   # each monkey independently oriented (no coherent interlock)\n"
        "# (no velocity key: clump insert samples U(-v,v); v=0 is an empty range)\n"
        f"region = {{ type = \"block\", min = [{CENTER[0]-ix:.6f}, {CENTER[1]-iy:.6f}, {CENTER[2]-iz:.6f}], "
        f"max = [{CENTER[0]+ix:.6f}, {CENTER[1]+iy:.6f}, {CENTER[2]+iz:.6f}] }}\n",
        runs_block(s, "rigid", gdot, smoke, phi, t_bond, float(radii.min())),
        common_tail(out_dir),
    ])
    _write(root, "rigid", phi, body)


def _place_gapped_monkeys(n, frac, centers, radii, rbound, seed):
    """Place N monkeys (random pos+orientation) into an aspect-2:2:1 box whose
    union solid fraction is `frac`, rejecting any monkey with a sub-sphere within
    GAP*(Ri+Rj) of an already-placed monkey's sub-sphere. Returns (rows|None, box_s)."""
    v_box = n * V_EQ / frac
    box_s = (v_box / DOMAIN_VOL)**(1.0 / 3.0)          # box = final × box_s
    hx, hy, hz = HALF_FINAL[0]*box_s, HALF_FINAL[1]*box_s, HALF_FINAL[2]*box_s
    ix = max(hx - rbound, 1e-3); iy = max(hy - rbound, 1e-3); iz = max(hz - rbound, 1e-4)
    rng = np.random.default_rng(seed)
    cell = GAP * 2.0 * float(radii.max())
    grid = {}
    lo = np.array([CENTER[0]-ix, CENTER[1]-iy, CENTER[2]-iz])
    hi = np.array([CENTER[0]+ix, CENTER[1]+iy, CENTER[2]+iz])

    def cidx(p):
        return (int(math.floor(p[0]/cell)), int(math.floor(p[1]/cell)), int(math.floor(p[2]/cell)))

    def conflicts(pos, rad):
        ci = cidx(pos)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for (q, rq) in grid.get((ci[0]+dx, ci[1]+dy, ci[2]+dz), ()):
                        d = math.sqrt((pos[0]-q[0])**2 + (pos[1]-q[1])**2 + (pos[2]-q[2])**2)
                        if d < GAP * (rad + rq):
                            return True
        return False

    rows = []
    placed = 0
    attempts = 0
    max_attempts = 400 * n + 50000
    while placed < n and attempts < max_attempts:
        attempts += 1
        com = lo + (hi - lo) * rng.random(3)
        R = random_rotation(rng)
        sub = com[None, :] + centers @ R.T
        if any(conflicts(sub[k], radii[k]) for k in range(len(sub))):
            continue
        for k in range(len(sub)):
            grid.setdefault(cidx(sub[k]), []).append((sub[k], radii[k]))
            rows.append((sub[k, 0], sub[k, 1], sub[k, 2], radii[k]))
        placed += 1
    return (rows if placed == n else None), box_s


def write_bpm(phi, n, gdot, smoke, root, centers, radii, rbound, t_bond):
    # Try the insertion cap; if the gap rule can't fit N monkeys, enlarge the box
    # (lower the fraction) and retry. box_s (>1 ⇒ compress) is set by whatever
    # fraction actually succeeded.
    frac0 = min(phi, CAP_MONKEY)
    placed_rows, box_s = None, 1.0
    for k in range(6):
        placed_rows, box_s = _place_gapped_monkeys(
            n, frac0 * (0.85 ** k), centers, radii, rbound,
            seed=20260704 + int(round(phi * 1000)) + k)
        if placed_rows is not None:
            break
    if placed_rows is None:
        raise RuntimeError(f"bpm Φ={phi}: gap placement failed even after enlarging the box")
    s = box_s
    hx, hy, hz, dom = domain_block(s)

    csv_rel = f"configs/bpm/phi_{phi}.csv"
    csv_path = os.path.join(root, csv_rel)
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    with open(csv_path, "w") as f:
        f.write("x,y,z,radius\n")
        for (x, y, z, r) in placed_rows:
            f.write(f"{x:.8e},{y:.8e},{z:.8e},{r:.8e}\n")

    out_dir = f"examples/monkey_shear/data/bpm/phi_{phi}"
    body = "\n".join([
        header_common("bpm (bonded 44-sphere monkey)", phi, n, s,
                      monkey_shear_dt(phi, t_bond, float(radii.min())), gdot,
                      extra=f"{len(placed_rows)} sub-spheres, gap {GAP}x(Ri+Rj)."),
        dom, neighbor_block("bpm", radii), material_block(),
        "[[particles.insert]]\n"
        'source = "file"\n'
        f"file = \"examples/monkey_shear/{csv_rel}\"\n"
        'format = "csv"\n'
        'material = "glass"\n'
        f"density = {GLASS['density']}\n"
        "columns = { x = 0, y = 1, z = 2, radius = 3 }\n",
        "[bonds]\n"
        "auto_bond = true\n"
        f"bond_tolerance = {BOND_TOL}\n"
        "bond_radius_ratio = 1.0\n"
        f"youngs_modulus = {E_BOND:.4e}\n"
        f"shear_modulus  = {G_BOND:.4e}\n"
        "# Critical damping (β=1) on all four channels keeps the welded monkey quasi-rigid:\n"
        "# it bleeds the sheared cantilever/arm-tip bond modes that otherwise resonate.\n"
        "beta_normal  = 1.0\nbeta_shear   = 1.0\nbeta_twist   = 1.0\nbeta_bending = 1.0\n"
        "# No [bonds.breakage] table => bonds are UNBREAKABLE (elastic) for this campaign.\n",
        runs_block(s, "bpm", gdot, smoke, phi, t_bond, float(radii.min())),
        common_tail(out_dir),
    ])
    _write(root, "bpm", phi, body)
    return len(placed_rows)


def _write(root, kind, phi, body):
    d = os.path.join(root, "configs", kind)
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, f"phi_{phi}.toml")
    with open(path, "w") as f:
        f.write(body.rstrip() + "\n")
    print(f"  wrote {os.path.relpath(path, root)}")


# ── Commands ─────────────────────────────────────────────────────────────────
def cmd_build_monkey():
    centers, radii, rep = scaled_monkey()
    hdr = (
        f"# monkey_Deq0.1.toml — the 44-sub-sphere \"monkey\" scaled to union D_eq = {D_EQ} m.\n"
        f"#\n"
        f"# Source: examples/monkey_barrel/settle_cad.toml (raw longest-axis 6 mm monkey),\n"
        f"# recentred on its sub-sphere centroid and scaled x{rep['scale']:.4f} so the\n"
        f"# Monte-Carlo union equivalent-volume diameter D_eq = (6 V_union/pi)^(1/3) = {D_EQ}.\n"
        f"#\n"
        f"# Monte-Carlo VERIFICATION (4e6 samples):\n"
        f"#   D_eq (raw, pre-scale) = {rep['d_eq_raw']:.6e} m   ->  scale = {rep['scale']:.5f}\n"
        f"#   D_eq (scaled)         = {rep['d_eq']:.6e} m   ({100*(rep['d_eq']/D_EQ-1):+.3f} % vs {D_EQ})\n"
        f"#   V_union               = {rep['v_union']:.6e} m^3  (V_eq target = {rep['v_eq']:.6e})\n"
        f"#   sub-sphere volume sum = {rep['v_sum']:.6e} m^3  -> overlap {100*rep['overlap']:.1f} %\n"
        f"#   sub-sphere radius     = min {rep['r_min']:.6e}  max {rep['r_max']:.6e} m\n"
        f"#   intra-monkey bonds @ tol {BOND_TOL} (bonds_per_monkey) = {rep['n_bonds']}\n"
        f"#   contact-Rayleigh(r_min)      = {t_rayleigh(rep['r_min']):.3e} s\n"
        f"#   bond-network critical dt     = {bond_network_critical_dt(centers, radii):.3e} s  (E_bond={E_BOND:.1e} weld)\n"
        f"#   -> monkey dt = {DT_SHEAR_SAFETY}·min(above)·Φ-margin  (see tools/gen_series.py);\n"
        f"#      e.g. Φ=0.05 → {monkey_shear_dt(0.05, bond_network_critical_dt(centers, radii), rep['r_min']):.2e}s,"
        f" Φ=0.4 → {monkey_shear_dt(0.4, bond_network_critical_dt(centers, radii), rep['r_min']):.2e}s (settle/shear).\n"
        f"#\n"
        f"# Reusable clump definition: [[clump.insert]] references it (rigid type); the bpm\n"
        f"# type expands the same sub-spheres to a per-monkey CSV via gen_series.py.\n"
    )
    with open(MONKEY_TOML, "w") as f:
        f.write(hdr + "\n" + clump_def_block(centers, radii) + "\n")
    print("== build-monkey ==")
    for k, v in rep.items():
        print(f"  {k:12s}: {v}")
    print(f"  wrote {os.path.relpath(MONKEY_TOML, EX_DIR)}")
    return centers, radii, rep


def cmd_series(gdot, smoke, only_phi=None):
    centers, radii, rep = scaled_monkey()
    root = EX_DIR
    t_bond = bond_network_critical_dt(centers, radii)
    r_min = float(radii.min())
    phis = [only_phi] if only_phi else PHIS
    print(f"== series (gdot={gdot}, {'SMOKE' if smoke else 'FULL'}) ==")
    print(f"   stiffest-element critical dt: contact-Rayleigh(r_min)={t_rayleigh(r_min):.3e}s, "
          f"bond-network={t_bond:.3e}s  → monkey dt = {DT_SHEAR_SAFETY}·min·Φ-margin")
    for phi in phis:
        n = N_TABLE[phi]
        dtm = monkey_shear_dt(phi, t_bond, r_min)
        print(f"[Φ={phi} N={n}]  monkey dt(settle/shear)={dtm:.3e}s compress={DT_COMPRESS_QS*dtm:.3e}s")
        write_sphere(phi, n, gdot, smoke, root)
        write_rigid(phi, n, gdot, smoke, root, centers, radii, rep["r_bound"], t_bond)
        nb = write_bpm(phi, n, gdot, smoke, root, centers, radii, rep["r_bound"], t_bond)
        print(f"    bpm sub-spheres = {nb}  (~{nb//n}/monkey), "
              f"expected bonds ≈ N*{rep['n_bonds']} = {n*rep['n_bonds']}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", nargs="?", default="all",
                    choices=["all", "build-monkey", "series"])
    ap.add_argument("--gdot", type=float, default=1.0,
                    help="Lees-Edwards shear rate γ̇ [1/s] (tuned: keeps monkeys in the "
                         "rigid-grain limit; mid-series Φ≈0.3-0.4 gives I≈0.05-0.1)")
    ap.add_argument("--smoke", action="store_true", help="short runs, Φ=0.1 only")
    ap.add_argument("--phi", type=float, default=None, help="restrict to a single Φ")
    args = ap.parse_args()

    only_phi = args.phi
    if args.smoke and only_phi is None:
        only_phi = 0.1

    if args.cmd in ("all", "build-monkey"):
        cmd_build_monkey()
    if args.cmd in ("all", "series"):
        cmd_series(args.gdot, args.smoke, only_phi)


if __name__ == "__main__":
    main()
