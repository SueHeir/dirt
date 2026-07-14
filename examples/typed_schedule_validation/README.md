# Typed scheduler contract

This executable schedule-integration check runs the five declared contracts in
[`data/schedule_contract.csv`](data/schedule_contract.csv) under the non-MPI,
precision-double configuration.  It verifies that a setforce-only setup remains
valid, addforce precedes the final setforce write, a missing contact provider is
rejected at schedule preparation, a supplied provider is accepted, and both
supported Hooke and Hertz contact choices publish the same typed seam.

Run `python3 examples/typed_schedule_validation/sweep.py`. The driver fails on an
unexpected test result and regenerates the tracked graph below. This is scheduler
contract coverage, not a DEM physics calibration.

![Observed typed scheduler contract outcomes; all five must pass](plots/typed_schedule_contract.svg)

Current result: PASS, 5/5 declared contracts match.
