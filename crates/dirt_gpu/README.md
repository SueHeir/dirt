# dirt_gpu

GPU contact kernels for DIRT, on [`grass_gpu`](https://github.com/SueHeir/grass/tree/main/crates/grass_gpu)
and `soil_gpu`. Kernels are written once in CubeCL: `f32` on Metal for debugging on a Mac, and
`f64` on CUDA with the `cuda` feature. Results and how to reproduce them are in
[`examples/gpu_contact`](../../examples/gpu_contact/README.md).

## Two ways to use it

| | `GranularResidentPlugins` | `GranularGpuPlugins` |
|---|---|---|
| On the GPU | the whole step | the contact force |
| Speed at 216,000 dense grains, on an M5 Pro | 18× one CPU core | 1.3× one CPU core |
| MPI | one rank | yes |
| CPU fixes | every `window` steps | every step |
| Walls | fixed infinite planes, on the GPU | anything, on the CPU |
| Use it when | the kernels cover the run's physics | it needs CPU physics or several ranks |

```rust,ignore
app.add_plugins(CorePlugins).add_plugins(GranularResidentPlugins { window: 100 });
// or
app.add_plugins(CorePlugins).add_plugins(GranularGpuPlugins);
```

With the resident path, a `[[run]]` stage's `steps` count windows, so simulated steps are
`steps * window`.

## What the kernels cover

Hertz-Mindlin with spring history and torque, constant-torque or `sds` rolling, constant-torque,
`sds` or `marshall` twisting, and any number of materials. Prefer the spring-based models, which
are what LAMMPS offers: the constant-torque ones cannot be reproduced grain by grain. Cohesion,
adhesion, liquid bridges, other contact or tangential models, and rigid clumps all panic at
the first step, so a run cannot silently lose physics.

Bonds: on the force-only path, add `DemBondPlugin` and everything it supports stays on the CPU, with
bonded pairs left out of GPU contact. On the whole-step path, `GranularResidentPlugins { bonds: true }`
runs elastic beams with `combined_stress` breakage or none, and refuses plasticity and other criteria.

Neither path adds to the virial, so per-contact stress output
needs the CPU force.

The resident path also leaves out three things. Springs are not written to restart files,
periodic image counters are not updated, and quaternions are not integrated. It also
refuses tilted or sheared boxes.

## Precision

`GpuReal` is `f32` on the Mac, since Apple GPUs have no `f64`, so compare with a CPU run within
the single-precision band there. The `cuda` feature selects `f64`. That path compiles but has
not yet run on an NVIDIA card, so treat it as untested.

## License

MIT OR Apache-2.0
