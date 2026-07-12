#!/usr/bin/env python3
"""Prepare and audit a Guo flexible-fibre shear-cell solver input.

This is intentionally not a synthetic rheology model.  It only writes the
particle and bond topology which a DIRT executable must consume, and records
the source parameters alongside it.  A campaign result still has to be made
by the solver and passed to ``validate.py``.
"""
import argparse
import csv
import json
import math
from pathlib import Path

# Guo et al., AIChE J. 65 (2019), doi:10.1002/aic.16397, Table 2.
DIAMETER = 0.0024
LENGTH = 0.0216
BEADS = 17
SPACING = 0.0012
DENSITY = 1157.5
YOUNGS_MODULUS = 6.28e6
TIMESTEP = 2.3e-7
BED_HEIGHT = 0.045
# The source gives a densely packed bed depth of approximately 45 mm.
INITIAL_LID_HEIGHT = 0.050
# The paper explicitly gives 500 rubber-cord fibres in the 64 x 36 mm control
# cell.  It reports a 96-mm-x sensitivity calculation but not its population;
# a 750-fibre population below preserves the cited 64-mm areal loading and is
# labelled as a derived size-sensitivity input, never as a Table-2 value.
FIBRES_AT_64_MM = 500
CELL_Z = 0.036


def write_case(width_mm, output):
    if width_mm not in (64, 96):
        raise ValueError("width must be 64 or 96 mm")
    width = width_mm / 1000.0
    fibres = FIBRES_AT_64_MM * width_mm // 64
    output.mkdir(parents=True, exist_ok=True)
    # Stagger long, x-aligned chains in separated lanes.  A planar x/z grid
    # would place collinear chains on top of one another because each chain is
    # 21.6 mm long; this layout keeps the source chains non-overlapping before
    # gravity settling.
    lanes, x_offset, lane_spacing = ((2, 0.010, 0.022)
                                     if width_mm == 64 else (3, 0.010, 0.028))
    rows = math.ceil(fibres / (lanes * 14))
    y_spacing = (BED_HEIGHT - DIAMETER) / (rows - 1)
    if y_spacing < DIAMETER:
        raise ValueError("fibre layers overlap")
    points = []
    for lane in range(lanes):
        for iy in range(rows):
            for iz in range(14):
                if len(points) // BEADS == fibres:
                    break
                for bead in range(BEADS):
                    points.append((x_offset + DIAMETER / 2 + lane * lane_spacing + bead * SPACING,
                                   DIAMETER / 2 + iy * y_spacing,
                                   (iz + 0.5) * CELL_Z / 14))
    if len(points) != fibres * BEADS:
        raise ValueError("incomplete fibre population")
    with (output / "pack.csv").open("w", newline="") as stream:
        csv.writer(stream).writerows(points)
    with (output / "pack.bonds").open("w") as stream:
        stream.write(f"{fibres * (BEADS - 1)} bonds\n1 bond types\n\nBonds\n\n")
        for fibre in range(fibres):
            for bead in range(BEADS - 1):
                stream.write(f"{fibre * (BEADS - 1) + bead + 1} 1 {fibre * BEADS + bead} {fibre * BEADS + bead + 1}\n")
    # Adjacent 2.4-mm spheres are separated by 1.2 mm, so summing their
    # volumes double-counts the material in their overlaps.  Fig. 7 reports
    # the physical rubber-cord volume fraction.  The Table-2 chain has the
    # same end-to-end length and diameter as a spherocylinder, whose volume
    # is the appropriate quantity for that comparison.  Keep the raw sphere
    # sum as an audit value, but never use it as experimental solid fraction.
    raw_bead_volume = len(points) * 4.0 * math.pi * (DIAMETER / 2.0) ** 3 / 3.0
    fibre_volume = math.pi * (DIAMETER / 2.0) ** 2 * (LENGTH - DIAMETER) + 4.0 * math.pi * (DIAMETER / 2.0) ** 3 / 3.0
    physical_cord_volume = fibres * fibre_volume
    planform_area = width * CELL_Z
    initial_cell_volume = planform_area * INITIAL_LID_HEIGHT
    metadata = {"paper": "Guo et al., AIChE J. 65 (2019), doi:10.1002/aic.16397, Table 2",
                "cell_x_m": width, "cell_z_m": CELL_Z, "planform_area_m2": planform_area,
                "n_fibres": fibres,
                "population_basis": "source Table 2" if width_mm == 64 else "derived: preserve 64-mm areal loading",
                "n_atoms": len(points), "n_bonds": fibres * (BEADS - 1),
                "bead_diameter_m": DIAMETER, "fibre_length_m": LENGTH,
                "beads_per_fibre": BEADS, "bond_spacing_m": SPACING,
                "density_kg_m3": DENSITY, "youngs_modulus_pa": YOUNGS_MODULUS,
                "timestep_s": TIMESTEP,
                # Retain both conventions: raw spheres audit the DEM input,
                # while the physical cord is the Fig.-7 comparison quantity.
                "raw_overlapping_bead_volume_m3": raw_bead_volume,
                "physical_cord_volume_m3": physical_cord_volume,
                "initial_lid_height_m": INITIAL_LID_HEIGHT,
                "raw_bead_solid_fraction_at_initial_lid": raw_bead_volume / initial_cell_volume,
                "physical_cord_solid_fraction_at_initial_lid": physical_cord_volume / initial_cell_volume}
    (output / "source_parameters.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def audit(output):
    metadata = json.loads((output / "source_parameters.json").read_text())
    points = list(csv.reader((output / "pack.csv").open()))
    if len(points) != metadata["n_atoms"] or metadata["n_atoms"] != metadata["n_fibres"] * BEADS:
        raise ValueError("atom count disagrees with source topology")
    for x, _, z in points:
        if not (0.0 <= float(x) < metadata["cell_x_m"] and 0.0 <= float(z) < metadata["cell_z_m"]):
            raise ValueError("a bead centre lies outside the periodic control cell")
    if not math.isclose((BEADS - 1) * metadata["bond_spacing_m"] + metadata["bead_diameter_m"], metadata["fibre_length_m"], abs_tol=1e-12):
        raise ValueError("bead chain does not reproduce the cited fibre length")
    bonds = (output / "pack.bonds").read_text().splitlines()[0]
    if bonds != f'{metadata["n_bonds"]} bonds':
        raise ValueError("bond count disagrees with source topology")
    expected_volume = metadata["n_atoms"] * 4.0 * math.pi * (metadata["bead_diameter_m"] / 2.0) ** 3 / 3.0
    if not math.isclose(metadata["raw_overlapping_bead_volume_m3"], expected_volume, rel_tol=1e-12):
        raise ValueError("represented bead volume disagrees with source topology")
    expected_physical_volume = metadata["n_fibres"] * (
        math.pi * (metadata["bead_diameter_m"] / 2.0) ** 2
        * (metadata["fibre_length_m"] - metadata["bead_diameter_m"])
        + 4.0 * math.pi * (metadata["bead_diameter_m"] / 2.0) ** 3 / 3.0)
    if not math.isclose(metadata["physical_cord_volume_m3"], expected_physical_volume, rel_tol=1e-12):
        raise ValueError("physical cord volume disagrees with Table-2 fibre geometry")
    expected_phi = expected_physical_volume / (metadata["planform_area_m2"] * metadata["initial_lid_height_m"])
    if not math.isclose(metadata["physical_cord_solid_fraction_at_initial_lid"], expected_phi, rel_tol=1e-12):
        raise ValueError("initial physical solid fraction disagrees with source topology")
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--width-mm", type=int, choices=(64, 96))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--audit", type=Path)
    args = parser.parse_args()
    if bool(args.width_mm is not None) == bool(args.audit is not None):
        parser.error("provide exactly one of --width-mm or --audit")
    if args.audit:
        metadata = audit(args.audit)
        print(f"PASS: {metadata['n_atoms']} beads, {metadata['n_bonds']} bonds, dt={metadata['timestep_s']:.2e} s")
    else:
        if args.output is None:
            parser.error("--output is required with --width-mm")
        metadata = write_case(args.width_mm, args.output)
        print(f"prepared {metadata['n_atoms']} beads and {metadata['n_bonds']} bonds")


if __name__ == "__main__":
    main()
