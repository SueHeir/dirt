# gpu_contact

One DEM run on any of three paths, for checking the GPU kernels against the CPU and
measuring their speed.

| `DIRT_CONTACT` | path |
|---|---|
| `cpu` | `GranularDefaultPlugins` |
| `gpu` | `GranularGpuPlugins`: contact force on the GPU, everything else on the CPU |
| `resident` | `GranularResidentPlugins`: the whole step on the GPU, `DIRT_WINDOW` steps per tick |

```sh
cargo build --release --features gpu --example gpu_contact
cargo build --release --no-default-features --features gpu,precision-single \
    --example gpu_contact --target-dir target-single
python3 examples/gpu_contact/validate.py
python3 examples/gpu_contact/speed.py
```

All numbers below are `f32` on an Apple M5 Pro (Metal), on 800 frictional grains settling in
a box for `validate.py`, and a dense periodic granular gas at solid fraction 0.5 for
`speed.py`.

## Agreement with the CPU

A settled pile is chaotic, so any change in rounding eventually sends grains elsewhere.
The fair question is whether the GPU differs from double-precision CPU by about as much as
single-precision CPU does. Largest difference in position over all grains, in metres:

| steps | CPU single | GPU force | GPU whole step |
|---:|---:|---:|---:|
| 1000 | 3.6e-7 | 9.9e-9 | 4.9e-6 |
| 2000 | 2.8e-4 | 2.6e-6 | 4.4e-4 |
| 5000 | 5.6e-3 | 5.3e-3 | 6.5e-3 |

Three things to know when reading this:

- **The whole step starts differently.** The CPU integrator's first half-kick uses a force
  of zero, because nothing has been computed yet. The GPU loop evaluates the force first.
  That leaves the CPU's velocity behind by exactly `g dt / 2`: measured 1.13e-4 m/s on one
  grain in free fall, calculated 1.134e-4. Position then drifts by that times `t`, which is
  the 1000-step entry. Contact timing shifts with it, and the difference joins the chaotic
  band by 2000 steps.
- **Atoms re-sort every 200 steps in these runs**, so springs surviving a sort is part of
  what this checks.
- **The test can see lost springs.** With history deliberately dropped at each rebuild, the
  GPU force's difference at 1000 steps goes from 1.0e-8 to 2.5e-5.

Rank counts, GPU force at 2000 steps, against one rank: 4.8e-7 on 2 ranks, 7.7e-8 on 4. That
is inside single-precision reach, though unlike the CPU path (about 1e-16) it is not
identical, because `f32` sums depend on the order atoms are visited.

Whole-step window length, at 2000 steps against a window of 500: 7.5e-7 for both 10 and 100.

## Speed

Simulated steps per second, best of two runs, with 1 rank unless noted:

| grains | CPU, 1 rank | CPU, 16 ranks | GPU force | GPU whole step | whole step / 1 rank | whole step / 16 ranks |
|---:|---:|---:|---:|---:|---:|---:|
| 8,000 | 2,520 | 12,853 | 1,409 | 13,420 | 5.3× | 1.04× |
| 27,000 | 737 | 5,011 | 605 | 7,400 | 10.0× | 1.48× |
| 64,000 | 283 | 2,244 | 287 | 3,900 | 13.8× | 1.74× |
| 216,000 | 68 | 518 | 90 | 1,230 | 18.0× | 2.37× |

The force-only path only passes one CPU core at about 64,000 grains. Every step waits on the
GPU, and neighbor build, communication and integration still run on the CPU. Its value is that it composes with everything the whole-step path cannot run.

These runs use a neighbor rebuild every 10 steps. The whole-step path used a window of 100.

## The f64 path

The same kernels build in `f64` on CubeCL's CPU runtime, which checks them without an NVIDIA
card: `cargo build --release --features gpu-cpu-f64 --example gpu_contact`. Against the CPU
contact force in double precision, position differs by 1.1e-15 at 1000 steps and 3.8e-14 at
2000, which is rounding. By 5000 steps chaos has grown that to 1.4e-4, the same growth rate
the single-precision runs show. This exercises the kernels' logic in `f64`, not CUDA's compiler
or hardware.

With two materials, walls and no gravity, where the CPU's first half kick loses nothing, the
whole-step path can be compared directly too. Largest position difference from the CPU force in
double precision, with `validate.py` printing this when the f64 build is present:

| steps | GPU force | GPU whole step |
|---:|---:|---:|
| 1000 | 4.5e-15 | 4.4e-15 |
| 2000 | 8.9e-14 | 6.8e-14 |
| 4000 | 1.9e-11 | 1.9e-11 |

## Rolling and twisting friction

These are checked differently, because the constant-torque models are singular. When two grains
spin alike their relative spin is rounding noise, yet the torque keeps its full size, with a
direction set by that noise. So any change in arithmetic moves trajectories from the first such
contact on, and a trajectory comparison says nothing. The CPU force does the same under a small
perturbation.

Instead: with rolling friction on, the f64 kernels' torque matched `contact_force_core_views`
given the same state, atom by atom, to 1e-9 relative, on every one of 1000 steps of a 400-grain
run. That was a one-off check with a temporary patch. Twisting has not had that check. `dirt_gpu`'s
tests cover both through the torque formulas on one pair, and the linear spin-down of a sphere
on a floor.

An earlier hand-written WGSL version of the same kernels reached 12.3× at 216,000 grains, so
CubeCL's version is faster here. None of this has run on CUDA yet.
