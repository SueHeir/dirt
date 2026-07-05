# peri_dem_interop

`peri_dem_interop` is the same-substrate dev_soil_peri+DIRT example for the
`dev-peri-dem-interop` / `ce-peri-dem-same-substrate-example` reconciliation.
It runs peridynamic fracture and DEM contact in one `App`, on one soil `Atom`
set, one neighbor list, and one shared `BondStore`. There is no coupler,
transport layer, or `grass_multi` sub-app.

Each lattice point is both a peridynamic point and a DEM sphere:

- `PeriPoint` is an example-local soil `AtomData` column for peridynamic volume
  and damage bookkeeping.
- `DemAtom` is DIRT's regular DEM atom column.
- Peridynamic bonds live in soil's `BondStore`; while a pair is still bonded,
  DIRT contact skips it through `BondStore::are_excluded`.
- When a peridynamic bond breaks, the example removes that bond entry. The pair
  then becomes eligible for ordinary Hertz-Mindlin DEM contact through the same
  neighbor list.

The configured case launches one brittle bar into another. DEM contact at the
unbonded interface drives stress waves, peridynamic bonds spall, and the
resulting fragments continue interacting through DEM contact. The hard gate is
closed-system conservation: mass and total linear momentum must remain conserved
through the peri-to-DEM transition. The sweep also requires that fracture
actually occurred and that active DEM contacts appeared after breakage.

Run the validation with the usual precision feature:

```bash
python3 examples/peri_dem_interop/sweep.py
```

Expected result:

```text
4/4 checks passed
ALL CHECKS PASSED
```
