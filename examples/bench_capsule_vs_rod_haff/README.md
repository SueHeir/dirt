# Haff Cooling — Single Capsule vs. Rod-of-Spheres

Validates the single-body **sphero-cylinder (capsule)** particle (`dirt_capsule`)
by cooling a periodic gas of them and comparing against a **rod-of-spheres
multisphere clump** (`dirt_clump`) of the same envelope. Both are elongated,
non-spherical particles that should obey Haff's law
`T(t) = T0 / (1 + t/tc)^2`, and their cooling curves should broadly agree.

## Two particle models, one shape

| | Capsule (`config_capsule.toml`) | Rod (`config_rod.toml`) |
|---|---|---|
| Representation | one rigid atom: cylinder r=0.5 mm, L=3 mm + hemispherical caps | 7 overlapping spheres (r=0.5 mm) at x = −1.5…+1.5 mm |
| Tip-to-tip | 4 mm (aspect ≈ 4) | 4 mm (aspect ≈ 4) |
| Rotation | full anisotropic inertia tensor, angular-momentum/Richardson | rigid-body clump integrator |
| Contact | capsule–capsule closest-segment + Hertz-Mindlin-Tsuji | sub-sphere Hertz-Mindlin-Tsuji |
| DOF per particle | 1 atom | 7 sub-sphere atoms |

Both use the **same** softened-glass material (E=70 MPa, ν=0.245, e=0.9, μ=0.3),
400 particles in a 40 mm periodic cube, random orientations, σ=0.5 m/s initial
velocity, no gravity — and the **same** Hertz-Mindlin-Tsuji force law and material
coefficients, so the comparison isolates *shape representation*.

## Run

```bash
for m in capsule rod; do
  CAPSULE_VS_ROD_MODE=$m cargo run --release --example bench_capsule_vs_rod_haff \
    --no-default-features --features precision-double -- \
    examples/bench_capsule_vs_rod_haff/config_$m.toml
done
```

Each writes `data/cooling_<mode>.csv` (`step,time,T_trans,T_rot,T_total`).

## Result

| Check | Capsule | Rod |
|-------|---------|-----|
| Haff's law: `1/√T` linear in `t`, R² | 0.998 | 0.9998 |
| No energy growth: `max(T)/T0` | 1.0000 | 1.0000 |
| Rotational equilibration (T_rot 0 → peak → decay) | ✓ | ✓ |

Both obey Haff's law and conserve energy. The two curves start identical (same
seed/initial field) and track within ~20% through the decay; the smooth capsule
cools slightly faster than the scalloped rod — the expected effect of surface
roughness on the tangential/rotational coupling. The capsule is the smooth
idealization of the bumpy rod, at 1/7 the atom count.
