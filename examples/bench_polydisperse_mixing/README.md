# bench_polydisperse_mixing — polydisperse / multi-material pair-mixing validation

Validates that DIRT's **per-pair mixing rules** enter the Hertz–Mindlin contact
force with the correct values, using single binary collisions between spheres of
**unequal radius** and/or **different material**. The mixing rules under test
(built in `dirt_atom::MaterialTable::build_pair_tables`) are:

| quantity | rule | where it enters |
| --- | --- | --- |
| reduced radius `R*` | `r1·r2 / (r1 + r2)` | Hertz spring `F = (4/3)E*√R* δ^{3/2}` |
| effective modulus `E*` (`e_eff_ij`) | `1 / ((1−ν1²)/E1 + (1−ν2²)/E2)` | Hertz spring |
| restitution `e_ij` | `√(e1·e2)` (geometric mean) → `beta_ij` | normal damping |
| sliding friction `μ_ij` (`friction_ij`) | `√(μ1·μ2)` (geometric mean) | tangential Coulomb cap |

The normal contact model itself is validated separately by
[`bench_hertz_rebound`](../bench_hertz_rebound); the tangential model by
[`bench_oblique_impact`](../bench_oblique_impact). **This** benchmark isolates the
*mixing* — it checks that for a heterogeneous pair the code combines the two
particles' properties correctly, not that a single-material contact is right.

## Method

Two families of single-collision cases (`sweep.py` builds one config each):

### Head-on (free–free) — `R*` and `E*`

Two **free** spheres collide along the line of centres. For **elastic** contacts
(`e = 1`, `μ = 0`) the collision is conservative and Hertz theory is exact:

```
δ_max = (15 m* v0² / (16 E* √R*))^(2/5)
t_c   = 2.868 (m*² / (R* E*² v0))^(1/5)
```

with `m* = m1 m2/(m1+m2)`, `R*`, `E*` as above (K. L. Johnson, *Contact
Mechanics*, Cambridge University Press, 1985, §11.4 & §11.5). Because these depend
on `R*` and `E*` **only through the mixing rules**, matching measured peak overlap
and contact duration for unequal-radius / cross-material pairs directly pins
`r_eff = r1 r2/(r1+r2)` and `e_eff_ij`. Cases isolate each rule:

- `N_R_*` — same material, radii from 5+5 to 2+8 mm → isolates `R*`.
- `N_E_*` — equal radii, soft/stiff material pairs → isolates `E*`.
- `N_RE_*` — unequal radii **and** cross material → both together.

### Restitution mixing — cross `(e1,e2)` vs a matched reference

The **realized** COR carries a known velocity/viscoelastic offset from the nominal
input `e` (see `bench_hertz_rebound`), so it is *not* compared to `√(e1 e2)`
directly — that would fold the calibration offset into a mixing check. Instead each
cross pair `(e1,e2)` is paired with a same-material **reference** at
`e_ref = √(e1 e2)`. Equal `e_ij` ⇒ equal `beta_ij` ⇒ identical realized COR, so
matching the two (`COR_A`, `COR_B`) isolates the geometric-mean rule and would
still fail if the code mixed arithmetically (the reference sits at the geometric
`e_ij`, which differs from the arithmetic mean printed alongside).

### Oblique (frozen target) — `μ_ij`

A projectile strikes an immovable **frozen** target sphere (`[[freeze]]`) with a
large tangential velocity chosen in the **gross-sliding** regime
(`v_t/v_n = 3`, well above the `(7/2)·μ_ij·(1+e)` stick threshold). There the
tangential force sits on the Coulomb cap for the whole contact, so the ratio of
tangential to normal impulse delivered to the projectile equals the pair friction:

```
|J_t| / |J_n| = μ_ij = √(μ1·μ2).
```

Impulses are read from the projectile's velocity change decomposed in the actual
line-of-centres frame captured at first contact (the frozen target keeps that
frame stable, exactly as in `bench_oblique_impact`). Cross pairs `(0.2,0.8)`,
`(0.1,0.9)`, `(0.16,0.64)` sit far from their **arithmetic** means, so the test
distinguishes geometric- from arithmetic-mean mixing; `(0.4,0.4)` is a
same-material sanity point; `F_02_08_46` combines unequal radii with friction.

## Running

```bash
python3 examples/bench_polydisperse_mixing/sweep.py            # generate + run + validate
python3 examples/bench_polydisperse_mixing/sweep.py graph      # re-validate existing data
```

`graph` prints per-case PASS/FAIL and exits non-zero if any check fails.

## Result

All cases pass. Elastic head-on peak overlap and contact duration match Hertz
theory (mixed `R*`, `E*`) to **≤ 0.1 %** across `R*` from 1.6–2.5 mm and `E*` from
4.8e9–3.7e10 Pa; elastic COR = 1.000. Cross-restitution collisions reproduce their
matched-`e_ij` reference COR to Δ ≤ 0.001. Oblique gross-sliding impulse ratios
match `√(μ1 μ2)` to **≤ 3.5 %** and lie far closer to the geometric than the
arithmetic mean.

![Mixing validation](plots/mixing_validation.png)

*Left: measured vs Hertz-theory peak overlap for the elastic head-on cases (mixed
`R*`, `E*`) — on the 1:1 line. Right: measured tangential/normal impulse ratio for
the oblique gross-sliding cases tracks the geometric mean `√(μ1 μ2)` (green), not
the arithmetic mean (red).*

## Where the test is weak / honest caveats

- **Shared-model check, within one code.** Like the rest of the tier-1 suite, the
  head-on and restitution references are analytical **for the model DIRT
  implements** (Hertz spring, geometric/effective mixing). They confirm the mixing
  math is wired up correctly; they do not adjudicate which mixing rule is
  physically "true" for a given material pair — that is a modelling choice.
- **Friction ~3 %, not exact.** The gross-sliding impulse ratio has a small
  systematic deficit (~2–3.5 %) below `μ_ij` because the tangential spring
  micro-slips near the end of contact (low normal load), so the cap is not held for
  the entire contact. The tolerance (5 %) reflects this; the sign is consistent and
  the ratio still clearly separates from the arithmetic mean.
- **Equal density.** All particles share one density, so `m*` varies through radius
  only. `R*` and `E*` are exercised across a wide range; the mass mixing `m*` is
  the standard reduced mass and is not separately stressed.
