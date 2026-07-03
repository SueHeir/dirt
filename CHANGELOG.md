# Changelog

All notable changes to DIRT are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

While DIRT is in the `0.1.x` pre-1.0 series, minor internal APIs may change
between releases without a major-version bump.

## [0.1.3] - 2026-07-03

### Changed
- Removed the experimental GPU (`dirt_gpu`) path from `main`; all crate
  dependencies now point at the Gitea remote (primary) instead of GitHub.
- Benchmark sweep drivers (`examples/*/sweep.py`) now pass `precision-double`
  explicitly so validation runs use the intended precision.

### Added
- `bench_twisting_friction`: twisting-friction spin-down benchmark validated
  against the analytical torque.

## [0.1.2] - 2026-06-23

### Added
- Compile-time precision abstraction (`Real`/`Accum`) with `precision-double`,
  `precision-mixed`, and `precision-single` feature flags.
- GPU contact-force kernels (`dirt_gpu`, wgpu/Metal) and a GPU granular-force
  plugin, including a CPU-vs-GPU comparison harness and pile examples.
- Coefficient-of-restitution (COR) calibration helpers and a calibration /
  benchmark example suite.
- Validation runners: CPU precision baseline (double/mixed/single), GPU-vs-CPU
  trajectory checks, and GPU-resident (windowed) + 2-rank MPI validation.
- Host↔device coherence path (`gpu_coherence`), kept opt-in and default-off
  until further validation.
- Documentation: guides, CI, and README updates.

### Changed
- Contact force factored into interior/boundary passes so interior compute
  overlaps halo communication.
- GPU-resident stepping keeps locals resident, writing only the ghost slice
  each tick; per-rank GPU binding added for MPI runs.

### Removed
- Dead `granular_basic` / `granular_gas_benchmark` example entries.

## [0.1.1] - 2026-06-15

### Added
- Preliminary strong- and weak-scaling studies and a scaling-gas Python driver.

### Fixed
- Damping bug fix and `grass` dependency resolution.
- Cargo.toml dependency-pinning fixes.

## [0.1.0] - 2026-06-14

### Added
- Initial DIRT DEM stack on top of the GRASS → SOIL → DIRT architecture.
- Core granular physics: Hertz–Mindlin contact, rotational dynamics, parallel
  bonds, walls (plane/cylinder/sphere/cone/region-surface), multisphere clumps,
  contact analysis, measurement planes, and DEM group fixes.
- Particle-insertion fix and pin/freeze controls.
- Example suite and initial documentation.

[0.1.3]: https://github.com/SueHeir/dirt/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/SueHeir/dirt/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/SueHeir/dirt/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/SueHeir/dirt/releases/tag/v0.1.0
