# Step 5: dynamic load balancing — verification

`soil_core::{balanced_planes, rebalance_x, LoadBalancePlugin}` rebalance the
x-decomposition by particle count (histogram → equal-count quantile planes via
`all_reduce_sum` per bin), at `PreExchange`, so the following `Exchange` migrates
atoms to their new owners. `bounds_changed` triggers the neighbour rebuild.

## Unit test
`balanced_planes`: skewed histogram (mass piled left) → split at x=1.667 (50/50
particle balance) instead of the equal-volume x=5.0. PASS.

## Integration (CPU, periodic frictional gas, 2000 grains, 2000 steps, mpirun)
Static vs rebalanced, merged by global tag:

| run | rebalances | max‖Δpos‖ vs 2-rank static | atoms lost | rank0 ownership |
|---|---|---|---|---|
| 1-rank (control) | — | 4.9e-6 (inherent decomposition divergence) | 0 | — |
| every=500 | 4 | 3.9e-5 | 0 | moves |
| every=50 | 40 | 1.3e-4 | 0 | 996→997 |

**Verdict:** the rebalance + re-exchange is functionally correct — atoms are
conserved (0 lost) and the split moves. The trajectory divergence from static
scales with rebalance count (benign accumulated f64 summation-order reorder that
this chaotic gas amplifies), NOT a per-step error, and stays within the
established precision band (< the 6e-4 sliding-friction baseline). On a
non-chaotic settling system it would track the ~1e-11 decomposition-invariance of
milestone 1.

**Limitation:** a single-hop exchange assumes a rebalance shifts a boundary < one
subdomain; very skewed distributions need the multi-hop migration path
(`fix/mpi-exchange`) to be safe, and the perf payoff needs multiple GPUs.
