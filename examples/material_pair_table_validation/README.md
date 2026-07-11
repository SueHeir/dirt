# Typed Material pair-table compatibility

This regression records the mixed **soft/stiff** pair values produced by the
pre-redesign `origin/main` `MaterialTable` implementation, then compares every
generated pair property with a fresh execution of this branch's typed
`MaterialTable` API. The golden table covers Hertz damping, elastic moduli, all
friction/adhesion modes, Hooke and SDS values, MDR values, and liquid-bridge
values. The corresponding Rust unit test independently transcribes the old
mixing rules and checks these same 21 properties.

Run `python3 examples/material_pair_table_validation/sweep.py` to compile and
execute the typed API's independent mixing-rule test, regenerate the figure,
and enforce the exact `1e-12` relative criterion. The sweep does not read a
static typed-output CSV.

![Typed material values versus legacy golden table](plots/typed_vs_legacy_pair_table.png)

*PASS: all 21 typed pair properties match the pre-redesign golden table within
the visible `1e-12` relative-error limit (the plotted values are exact, shown at
the `1e-18` plotting floor).*
