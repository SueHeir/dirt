#!/usr/bin/env python3
"""Materialise non-calibrated sphere-built Guo wall/blade candidates.

This is geometry evidence, not a DIRT result.  Guo et al. specify sphere-built
walls, blade lengths and pitch, but not the wall-sphere diameter or lattice.
The three resolutions are therefore a predeclared *non-equivalent* sensitivity
set.  They may not be selected using Fig. 6/7, nor called the paper's mesh.
"""
import argparse, csv, json, math
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONTRACT = HERE / "data" / "source_geometry_contract.json"
RESOLUTIONS_MM = (0.6, 1.2, 2.4)

def centres(length, spacing):
    n = math.ceil(length / spacing)
    return [(i + .5) * length / n for i in range(n)]

def make(diameter_mm, output):
    if diameter_mm not in RESOLUTIONS_MM:
        raise ValueError("diameter must be one of the predeclared resolutions")
    d = diameter_mm / 1000; r = d / 2
    x, z = .064, .036
    # A one-sphere-thick square lattice is explicit and reproducible.  It is
    # deliberately not inferred from the paper.
    lower = [(a, r, b) for a in centres(x, d) for b in centres(z, d)]
    upper = [(a, .050-r, b) for a in centres(x, d) for b in centres(z, d)]
    # Blades occur every 8 mm along x and extend into the bed by the reported
    # length.  Cross-z rows make them periodic with the cell and are rigidly
    # attached by construction in any future clump/wall implementation.
    blades = []
    for y0, sign, length in ((r, 1, .002), (.050-r, -1, .004)):
        for a in centres(x, .008):
            for b in centres(z, d):
                for k in range(1, math.ceil(length/d)+1):
                    blades.append((a, y0 + sign*k*d, b))
    output.mkdir(parents=True, exist_ok=True)
    with (output / 'wall_spheres.csv').open('w', newline='') as f:
        w=csv.writer(f); w.writerow(('x_m','y_m','z_m','part')); w.writerows((*p,'lower') for p in lower); w.writerows((*p,'upper') for p in upper); w.writerows((*p,'blade') for p in blades)
    manifest = {'label': f'predeclared square-lattice {diameter_mm:g}-mm candidate',
      'source_equivalent': False, 'selection_basis': 'predeclared_sensitivity',
      'diameter_mm': diameter_mm,
      'layout': 'one-sphere-thick square lattices; x-pitch 8-mm blade rows; lower/upper blade lengths 2/4 mm',
      'reference_observables_consulted': [], 'counts': {'lower':len(lower),'upper':len(upper),'blade':len(blades)},
      'source_constraints': {'periodic_cell_mm':[64,36],'blade_pitch_mm':8,'lower_blade_length_mm':2,'upper_blade_length_mm':4}}
    (output/'wall_realisation.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--diameter-mm',type=float,required=True); p.add_argument('--output',type=Path,required=True); a=p.parse_args()
    m=make(a.diameter_mm,a.output)
    print(f"CANDIDATE ONLY: {m['label']} ({sum(m['counts'].values())} spheres)")
if __name__=='__main__': main()
