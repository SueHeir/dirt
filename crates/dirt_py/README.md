# dirt_py (preview, test branch)

`import dirt`: build and run DIRT simulations from Python, and couple them
with Apps from other GRASS modules.

This branch builds DIRT against a local grass checkout on its
`python-preview` branch (see the `[patch]` section in the workspace
`Cargo.toml`), because the shared Python App layer (`grass_pyo3`) lives
there.

```sh
# grass and dirt checked out side by side, both on python-preview
cargo build --release -p grass_py            # in grass/  -> libgrass.so
cargo build --release -p dirt_py             # in dirt/   -> libdirt.so
cp grass/target/release/libgrass.so mods/grass.abi3.so
cp dirt/target/release/libdirt.so  mods/dirt.abi3.so
PYTHONPATH=mods python dirt/crates/dirt_py/python/fluidized_bed.py
```

* `dirt.App` has the same API as `grass.App` (typed resources, phases,
  Python plugins, checkpoints).
* `dirt.plugins` adds `Core`, `GranularDefault`, `Gravity` and `Wall`.
* `dirt.phases` exposes DIRT's schedule phases, so a Python system can run
  in `PostForce`.
* `dirt.Particles` gives a Python system the local particles as NumPy
  arrays, and adds its `f_ext` back into the particle forces.
* A `dirt.App` can be a child of a `grass.App` and the other way round.
  The two modules are separate shared libraries, each with its own copy of
  GRASS, so they meet through the stable C ABI in `grass_pyo3::abi`. Data
  crosses as Python objects.

`python/fluidized_bed.py` couples a DIRT bed of 200 spheres with a Python
fluid column, using grass.Map in both directions. The bed expands as the
inlet velocity rises.
