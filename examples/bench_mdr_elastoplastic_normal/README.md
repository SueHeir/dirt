# MDR Elastic-Plastic Normal Contact

This benchmark exercises DIRT's selectable `contact_model = "mdr"` normal branch with a quasi-static two-sphere loading/unloading path. It validates the normal force trace against the rigid-flat particle-pair equations used by LAMMPS `GranSubModNormalMDR::calculate_forces` (`src/GRANULAR/gran_sub_mod_normal.cpp`).

The DIRT implementation follows LAMMPS's per-side rigid-flat placement for particle pairs, first yield, plastic unloading with the `deltaR` correction, MDR contact-radius stiffness, and an adhesive tensile branch. It deliberately does not yet update apparent particle radii, include the multi-contact free-surface/bulk response from LAMMPS `fix GRANULAR/MDR`, apply the topological contact-penalty screen, or implement MDR walls; those differences are documented in the LAMMPS parity table.

![MDR normal force trace](plots/mdr_force_trace.png)

The figure shows the DIRT loading/unloading force-displacement trace overlaid on the LAMMPS-source reference with the tolerance band used by the gate. Current result: PASS.

Run:

```bash
python3 examples/bench_mdr_elastoplastic_normal/sweep.py
```
