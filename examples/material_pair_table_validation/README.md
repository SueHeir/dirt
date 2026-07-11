# Typed Material pair-table compatibility

This regression records the mixed **soft/stiff** pair values produced by the
pre-redesign `origin/main` `MaterialTable` implementation, then compares every
generated pair property with a fresh execution of this branch's typed
`MaterialTable` API. The golden table covers Hertz damping, elastic moduli, all
friction/adhesion modes, Hooke and SDS values, MDR values, and liquid-bridge
values. The command is a standalone executable that constructs the inputs only
through the public typed API and emits its generated pair table; it is not a
unit-test printout or a static typed-output fixture.

Run `python3 examples/material_pair_table_validation/sweep.py` to compile and
execute the typed API validation executable, regenerate the figure,
and enforce the exact `1e-12` relative criterion. The sweep does not read a
static typed-output CSV.

![Typed material values versus legacy golden table](plots/typed_vs_legacy_pair_table.png)

*PASS: all 21 typed pair properties match the pre-redesign golden table within
the visible `1e-12` relative-error limit (the plotted values are exact, shown at
the `1e-18` plotting floor).*
