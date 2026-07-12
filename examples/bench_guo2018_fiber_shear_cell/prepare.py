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
DEPTH = 0.036
BED_HEIGHT = 0.045
# This is the lid gap in config.toml.  It is deliberately an input to the
# audit rather than a fitted quantity: the source topology must be capable of
# reaching the external solid-fraction range before spending a long solver run.
INITIAL_LID_HEIGHT = 0.050
LANES = {64: (500, 2), 96: (750, 3)}


def write_case(width_mm, output):
    if width_mm not in LANES:
        raise ValueError("width must be 64 or 96 mm")
    fibres, lanes = LANES[width_mm]
    width = width_mm / 1000.0
    output.mkdir(parents=True, exist_ok=True)
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
                    points.append((DIAMETER / 2 + lane * 0.032 + bead * SPACING,
                                   DIAMETER / 2 + iy * y_spacing,
                                   (iz + 0.5) * DEPTH / 14))
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
    initial_cell_volume = width * DEPTH * INITIAL_LID_HEIGHT
    metadata = {"paper": "Guo et al., AIChE J. 65 (2019), doi:10.1002/aic.16397, Table 2",
                "width_m": width, "depth_m": DEPTH, "n_fibres": fibres,
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
    expected_phi = expected_physical_volume / (metadata["width_m"] * DEPTH * metadata["initial_lid_height_m"])
    if not math.isclose(metadata["physical_cord_solid_fraction_at_initial_lid"], expected_phi, rel_tol=1e-12):
        raise ValueError("initial physical solid fraction disagrees with source topology")
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--width-mm", type=int, choices=LANES)
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
