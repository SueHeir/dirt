# bench_mindlin_rescale_tangential — Mindlin unloading rescale

Validates the selectable `tangential_model = "mindlin_rescale"` and
`"mindlin_rescale/force"` variants against the documented LAMMPS unloading rule
in `doc/src/pair_granular.rst`: when the normal contact unloads and the Hertz
contact radius decreases, the tangential history is scaled by `a/a_prev`.

## Setup

Two identical spheres are held at prescribed overlaps. The contact is loaded
tangentially at fixed peak overlap, then the normal overlap is reduced while
`v_t = 0`. With `restitution = 1` and a high friction cap, damping and Coulomb
sliding are inactive, so the distinguishing recurrence is exact:

```
history:               xi stays constant, so F_t scales as a
mindlin_rescale:       xi <- xi * a/a_prev, so F_t scales as a^2
mindlin_rescale/force: Fte <- Fte * a/a_prev
linear_nohistory:      no stored history, so F_t = 0 when v_t = 0
```

![Mindlin unloading rescale](plots/mindlin_rescale_unload.png)

*DIRT unloading forces against the documented LAMMPS recurrence. PASS means every
model matches its recurrence to numerical precision, `mindlin_rescale` drops
quadratically relative to displacement-history Mindlin, and `linear_nohistory`
remains zero during the zero-slip unload.*
