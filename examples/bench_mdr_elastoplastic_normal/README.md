# MDR Elastic-Plastic Normal Contact

This benchmark exercises DIRT's selectable `contact_model = "mdr"` normal branch with a quasi-static two-sphere loading/unloading path. It validates the normal force trace against an independent Python implementation of the documented one-dimensional MDR subset in `sweep.py`.

The DIRT implementation follows the Zunker/Kamrin/LAMMPS MDR shape transform for elastic loading, first yield, plastic unloading, MDR contact-radius stiffness, and an adhesive tensile branch. It deliberately does not yet update apparent particle radii or the multi-contact free-surface/bulk response used by LAMMPS `fix GRANULAR/MDR`; those differences are documented in the LAMMPS parity table.

![MDR normal force trace](plots/mdr_force_trace.png)

The figure shows the DIRT loading/unloading force-displacement trace overlaid on the analytic reference with the tolerance band used by the gate. Current result: PASS.

Run:

```bash
python3 examples/bench_mdr_elastoplastic_normal/sweep.py
```
