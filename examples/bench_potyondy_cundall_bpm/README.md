# bench_potyondy_cundall_bpm — Potyondy-Cundall BPM rock compression

This benchmark runs a small bonded-particle compression specimen and validates
its normalized stress-strain response, peak strength, failure strain, and
bond-break progression against Potyondy & Cundall's Lac du Bonnet granite BPM
calibration.

The specimen and loading reduction live entirely in this example. Each load
increment evaluates DIRT's `CombinedStress` breakage criterion, the
Potyondy-Cundall extreme-fibre stress form cited in `dirt_bond`, for every
surviving bond. Broken bonds are removed from the specimen graph and plotted as
crack progression.

## Run

```bash
python3 examples/bench_potyondy_cundall_bpm/sweep.py
```

The committed target file,
`data/potyondy_cundall_2004_fig8a_digitized.csv`, is an approximate
normalization of Potyondy & Cundall (2004) Fig. 8(a), low-confinement PFC2D
Lac du Bonnet granite, using Table 2's PFC2D `q_u = 199.1 MPa` and
`E = 70.9 GPa` (`q_u/E = 0.00281`) as the peak stress and strain scales.

**Gate:** normalized peak strength must be within 12% of Table 2 PFC2D `q_u`,
peak/failure strain must be within 18% of `q_u/E`, and the run must produce a
near-peak-to-post-peak bond-break progression with at least 20 broken bonds.

![Stress-strain and crack progression](plots/stress_strain_and_cracks.png)

*DIRT reduced BPM stress-strain curve against the digitized Potyondy-Cundall
Fig. 8(a) target, plus the spatial sequence of CombinedStress bond breaks.
Latest run: PASS.*

## Reference

Potyondy, D. O. & Cundall, P. A. (2004). "A bonded-particle model for rock."
*International Journal of Rock Mechanics and Mining Sciences*, 41(8),
1329-1364. Table 2 reports the Lac du Bonnet PFC2D model macroproperties used
here (`E = 70.9 GPa`, `q_u = 199.1 MPa`), and Fig. 8(a) shows the biaxial
stress-strain response and post-peak damage patterns used for the normalized
curve target.
