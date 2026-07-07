# `bench_wall_activate_by_name`

This validation checks the runtime named-wall control path for a wall that is
needed again after being removed. A single sphere is held at fixed overlap with a
named plane wall. The example samples the wall force while the wall is active,
calls `Walls::deactivate_by_name("gate")`, then calls
`Walls::activate_by_name("gate")` on the same resource and samples again.

The expected response is binary because geometry and material state are unchanged:
the active and reactivated windows must produce the same nonzero normal force,
and the deactivated window must be force-free.

![Named wall force response](plots/wall_activate_by_name_force.png)

*Measured particle-wall normal force during active, deactivated, and reactivated
windows. Latest run: PASS, the inactive force is zero within `1e-14` N and the
reactivated mean force matches the initial active force within `1e-12` relative
error.*

Run it with:

```bash
python3 examples/bench_wall_activate_by_name/sweep.py
```
