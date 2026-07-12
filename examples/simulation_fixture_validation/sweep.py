#!/usr/bin/env python3
"""Validate public SimulationFixture measurements against a fixed contract."""
import csv
import math
import subprocess
from pathlib import Path

ROOT = Path(__file__).parent
REPO = ROOT.parent.parent
TOL = 1e-15
FIELDS = ("atom_rows", "nlocal", "natoms", "dem_rows", "materials", "pair_rows", "pair_columns", "dt")


def write_plot(labels, errors):
    """Emit the committed contract graph without a third-party Python package.

    The measured values are normally exact, so they are displayed at a visible
    1e-18 floor.  The dashed red line is the real 1e-15 failure threshold;
    values above it still make this driver exit non-zero below.
    """
    plot = ROOT / "plots" / "fixture_contract.svg"
    plot.parent.mkdir(exist_ok=True)
    width, height = 1240, 560
    left, right, top, bottom = 82, 28, 66, 164
    axis_bottom = height - bottom
    axis_height = axis_bottom - top
    x_step = (width - left - right) / len(errors)
    log_floor, log_ceiling = -18, 0

    def y_for(error):
        # Exact matches use a display floor only; the validation comparison
        # remains against the unmodified error values.
        exponent = math.log10(max(error, 1e-18)) if math.isfinite(error) else log_ceiling
        exponent = min(max(exponent, log_floor), log_ceiling)
        return top + (log_ceiling - exponent) / (log_ceiling - log_floor) * axis_height

    threshold_y = y_for(TOL)
    grid = []
    for exponent in (-18, -15, -12, -9, -6, -3, 0):
        y = y_for(10.0**exponent)
        grid.append(
            f'<line x1="{left}" y1="{y:.1f}" x2="{width-right}" y2="{y:.1f}" '
            f'stroke="#d9d9d9" stroke-width="1"/>'
            f'<text x="{left-10}" y="{y+5:.1f}" text-anchor="end" '
            f'font-family="sans-serif" font-size="12">1e{exponent}</text>'
        )

    points, ticks = [], []
    for index, (label, error) in enumerate(zip(labels, errors)):
        x = left + (index + 0.5) * x_step
        y = y_for(error)
        color = "#2f6f73" if error <= TOL else "#b22222"
        points.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{color}"/>')
        ticks.append(
            f'<text x="{x+3:.1f}" y="{axis_bottom+12}" transform="rotate(62 {x+3:.1f} {axis_bottom+12})" '
            f'font-family="sans-serif" font-size="10">{label.replace("&", "&amp;")}</text>'
        )

    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="white"/>
<text x="{width / 2}" y="30" text-anchor="middle" font-family="sans-serif" font-size="21">SimulationFixture structural contract through its public API</text>
<text x="{width / 2}" y="52" text-anchor="middle" font-family="sans-serif" font-size="13">20 measured scalar and CSR checks; exact matches are plotted at the 1e-18 floor</text>
{''.join(grid)}
<line x1="{left}" y1="{axis_bottom}" x2="{width-right}" y2="{axis_bottom}" stroke="black" stroke-width="1.5"/>
<line x1="{left}" y1="{top}" x2="{left}" y2="{axis_bottom}" stroke="black" stroke-width="1.5"/>
<line x1="{left}" y1="{threshold_y:.1f}" x2="{width-right}" y2="{threshold_y:.1f}" stroke="crimson" stroke-width="2" stroke-dasharray="7 5"/>
<text x="{width-right}" y="{threshold_y-8:.1f}" text-anchor="end" font-family="sans-serif" font-size="13" fill="crimson">PASS limit: 1e-15 relative error</text>
<text x="24" y="{(top + axis_bottom) / 2}" transform="rotate(-90 24 {(top + axis_bottom) / 2})" text-anchor="middle" font-family="sans-serif" font-size="14">relative error (log scale)</text>
{''.join(points)}
{''.join(ticks)}
<circle cx="{left}" cy="{height-28}" r="4" fill="#2f6f73"/><text x="{left+10}" y="{height-23}" font-family="sans-serif" font-size="13">measurement within the contract limit</text>
<circle cx="{left+350}" cy="{height-28}" r="4" fill="#b22222"/><text x="{left+360}" y="{height-23}" font-family="sans-serif" font-size="13">measurement above the contract limit (FAIL)</text>
</svg>'''
    plot.write_text(svg, encoding="utf-8")

def read_rows(text):
    return {row["case"]: row for row in csv.DictReader(text.splitlines())}

expected = read_rows((ROOT / "data" / "fixture_contract.csv").read_text())
result = subprocess.run(
    ["cargo", "run", "--quiet", "--example", "simulation_fixture_validation", "--no-default-features", "--features", "precision-double"],
    cwd=REPO, check=True, text=True, capture_output=True,
)
measured = read_rows(result.stdout)
if expected.keys() != measured.keys():
    raise SystemExit(f"case mismatch: expected {expected.keys()}, measured {measured.keys()}")

errors, labels = [], []
for case in expected:
    for field in FIELDS:
        want, got = float(expected[case][field]), float(measured[case][field])
        errors.append(abs(got - want) / max(abs(want), 1.0))
        labels.append(f"{case}\n{field}")
    for field in ("csr_offsets", "csr_indices"):
        errors.append(0.0 if measured[case][field] == expected[case][field] else math.inf)
        labels.append(f"{case}\n{field}")

failed = [label for label, error in zip(labels, errors) if error > TOL]
write_plot(labels, errors)
print(f"SIMULATION FIXTURE: {'PASS' if not failed else 'FAIL'} — {len(errors) - len(failed)}/{len(errors)} structural measurements match the committed contract")
if failed:
    raise SystemExit("failed: " + ", ".join(failed))
