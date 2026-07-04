# bench_hooke_rebound — Hooke (linear-spring) normal-contact rebound

Exercises DIRT's **linear spring-dashpot** normal contact
(`contact_model = "hooke"`, per-material `kn`/`kt`) — the contact-force branch in
[`crates/dirt_granular/src/contact.rs`](../../crates/dirt_granular/src/contact.rs)
that every other benchmark leaves untouched, because they all use the nonlinear
**Hertz** model. Two identical spheres collide head-on and rebound; the benchmark
gates the measured coefficient of restitution, contact duration, and peak overlap
against the **exact** analytical collision.

## Why the linear contact has an exact reference

During contact the mutual overlap `x(t)` of two spheres obeys a linear
spring-dashpot with exactly the coefficients DIRT integrates
(`contact.rs`: `f_n = kn·δ − γ_n·v_n`, `γ_n = 2β·√(kn·m_eff)`):

```
m_eff·ẍ + γ_n·ẋ + kn·x = 0 ,   m_eff = m/2   (two identical spheres)
```

That is a **constant-coefficient damped harmonic oscillator** — the one contact
law with a closed-form collision (unlike Hertz, whose nonlinear stiffness has no
elementary damped solution). With

```
ω₀ = √(kn/m_eff) ,   ζ = γ_n/(2√(kn·m_eff)) = β ,   ω_d = ω₀·√(1−β²)
```

and DIRT deriving the damping ratio from the restitution input `e` as the **exact
linear inversion** (`dirt_atom::build_pair_tables`, Hooke branch)

```
β = −ln e / √(π² + ln²e)
```

the standard oscillator results reduce to closed forms with **no free constants**:

| quantity | measured from DIRT | linear-contact analytical reference |
|---|---|---|
| coefficient of restitution `COR` | \|v_rebound\| / \|v_impact\| (relative, normal) | `exp(−π·ζ/√(1−ζ²)) = e` (exact) |
| contact duration `t_c` | steps in contact × `dt` | `π/ω_d = √(π²+ln²e)·√(m_eff/kn)` |
| peak overlap `δ_max` | peak geometric overlap `R₁+R₂−d` | `(v/ω_d)·e^(−ζω₀t*)·sin(ω_d t*)`, `t*=atan(√(1−β²)/β)/ω_d` |

**COR and `t_c` are velocity-independent** — the defining signature of a *linear*
contact (the undamped Hertz `t_c ∝ v^{−1/5}`, by contrast). `δ_max` scales linearly
with the impact speed. All three, plus the two velocity-independence properties, are
checked. The reference is **theory only** — neither DIRT's own output nor another
code (anti-gaming: analytical, not self-consistent).

## Running

```bash
# one representative case
cargo run --release --example bench_hooke_rebound --no-default-features \
    --features precision-double -- examples/bench_hooke_rebound/config.toml

# full sweep (generate → build → run → validate + plot)
python3 examples/bench_hooke_rebound/sweep.py
```

`sweep.py graph` prints `ALL CHECKS PASSED` and exits 0 when every case is within
tolerance (COR 1 %, contact-time 2 %, overlap 2 %, plus velocity-independence).

## Result

Across the 5 × 4 sweep (restitution 0.3–1.0 × relative impact velocity 0.5–4 m/s)
DIRT matches every analytical quantity to **≤ 0.05 %**: COR to within 0.0005 of the
input `e`, contact time and peak overlap to a few 0.01 %, and both COR and `t_c` are
flat across impact speed to < 0.05 %. This makes it the strongest normal-contact
check in the suite — it validates the linear stiffness `kn`, the restitution→damping
derivation, and the integrator simultaneously against an exact damped closed form,
not just an elastic limit or a calibrated mapping.

![Measured vs input COR](plots/cor_validation.png)
![Contact duration](plots/contact_duration.png)
![Peak overlap](plots/peak_overlap.png)

## References

- Y. Tsuji, T. Tanaka, T. Ishida, "Lagrangian numerical simulation of plug flow of
  cohesionless particles in a horizontal pipe", *Powder Technology* **71**:239–250
  (1992) — linear COR ↔ damping relation `β = −ln e/√(π²+ln²e)`.
- H.-G. Schäfer, S. Dippel, D.E. Wolf, "Force schemes in simulations of granular
  materials", *J. Phys. I France* **6**:5–20 (1996) — linear spring-dashpot contact
  duration and restitution, Eqs. (2.10)–(2.16).
