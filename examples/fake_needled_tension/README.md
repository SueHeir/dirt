# Fake-needled felt tension

This is the paired through-thickness tensile test for the 20-ply fake felt.
It uses the same DIRT driver, materials, grip fraction, loading history, MPI
decomposition, timestep, and output fields as `fake_felted_tension`. The only
intended specimen difference is TANGLE's needling continuation.

First generate the needled specimen layer by layer. Every ply follows the
control insertion/lowering/relaxation sequence, then the newly deposited ply
is needled before the next one is inserted:

```console
cd ../tangle
cargo run --release -p fake_felted -- --needled
```

If that continuation is interrupted, resume it without discarding progress:

```console
cargo run --release -p fake_felted -- --resume-needled
```

After layer construction, run the same staged contact/curvature cleanup and
controlled compaction used to prepare the unneedled control for DIRT:

```console
cargo run --release -p fake_felted -- --needled-cleanup
```

Resume an interrupted cleanup with `--resume-needled-cleanup`.

Finally reduce the capsule overlap to the stricter DIRT startup tolerance:

```console
cargo run --release -p fake_felted -- --needled-dem-polish
```

Then run the same eight-rank DIRT pull from the DIRT repository root:

```console
cargo build --release --example fake_felted_tension
TENSION_SETTLE_STEPS=1000 \
TENSION_RAMP_STEPS=5000 \
TENSION_GRIP_SPEED=0.2 \
TENSION_RECORD_EVERY=100 \
mpirun -np 8 target/release/examples/fake_felted_tension \
  examples/fake_needled_tension/config_hour.toml
```

Plot the needled curve and, when available, overlay the unneedled control:

```console
python3 examples/fake_needled_tension/plot.py
```

## Paired lower-rate pull

The two `config_slow.toml` inputs preserve the original 400,000 steps and
timestep but write to separate `output_slow` directories. Setting the final
engineering strain to 1% makes DIRT derive the required grip speed separately
from each specimen's initial grip spacing. This gives both specimens the same
engineering strain rate and prevents either original dataset from being
overwritten.

Run the unneedled control, then the needled specimen:

```console
TENSION_SETTLE_STEPS=1000 TENSION_RAMP_STEPS=5000 \
TENSION_FINAL_STRAIN_PERCENT=1.0 TENSION_RECORD_EVERY=100 \
mpirun -np 8 target/release/examples/fake_felted_tension \
  examples/fake_felted_tension/config_slow.toml

TENSION_SETTLE_STEPS=1000 TENSION_RAMP_STEPS=5000 \
TENSION_FINAL_STRAIN_PERCENT=1.0 TENSION_RECORD_EVERY=100 \
mpirun -np 8 target/release/examples/fake_felted_tension \
  examples/fake_needled_tension/config_slow.toml
```

Plot only this matched pair with:

```console
python3 examples/fake_needled_tension/plot.py --slow
```

## Fifty-percent HPC pull

The `config_50pct.toml` pair runs two million steps to reach 50% engineering
strain measured against the initial internal, non-gripped gauge length. The
driver records every 2,000 steps, producing approximately 1,000 samples over
the complete job.

Build in a configured HPC shell, or name an optional site-local environment
file explicitly:

```console
DIRT_BUILD_ENV=/path/to/site-build-env \
  bash examples/fake_needled_tension/hpc_build_50pct.sh
```

On a 32-core node, run the matched cases concurrently on disjoint CPU sets:

```console
bash examples/fake_needled_tension/hpc_run_50pct.sh unneedled 0-15 &
bash examples/fake_needled_tension/hpc_run_50pct.sh needled 16-31 &
wait
```

Each case uses sixteen MPI ranks and refuses to overwrite an existing
`tension.csv`. Results are written to the corresponding `output_50pct`
directory, which is intentionally excluded from version control.
