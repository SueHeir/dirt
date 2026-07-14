#!/usr/bin/env python3
"""Execute declared DIRT scheduler contracts and render their actual outcome."""
import csv
import html
import subprocess
from pathlib import Path

ROOT = Path(__file__).parent
REPO = ROOT.parent.parent
PLOT = ROOT / "plots" / "typed_schedule_contract.svg"


def run_case(case):
    command = [
        "cargo", "test", "--quiet", "-p", case["crate"],
        "--no-default-features", "--features", "precision-double", case["test"],
    ]
    result = subprocess.run(command, cwd=REPO, text=True, capture_output=True)
    case["observed"] = "pass" if result.returncode == 0 else "fail"
    if result.returncode:
        print(result.stdout, result.stderr)


with (ROOT / "data" / "schedule_contract.csv").open(newline="") as source:
    cases = list(csv.DictReader(source))
for case in cases:
    run_case(case)

failed = [case["case"] for case in cases if case["observed"] != case["expected"]]
width, left, top, row = 1060, 280, 88, 72
height = top + row * len(cases) + 92
bars = []
for index, case in enumerate(cases):
    y = top + index * row
    observed = case["observed"]
    color = "#197a3d" if observed == "pass" else "#b42318"
    status = "PASS" if observed == "pass" else "FAIL"
    label = html.escape(case["case"].replace("_", " "))
    bars.append(
        f'<text x="24" y="{y + 28}" class="case">{label}</text>'
        f'<rect x="{left}" y="{y}" width="700" height="42" fill="#f1f5f9"/>'
        f'<rect x="{left}" y="{y}" width="700" height="42" fill="{color}"/>'
        f'<text x="{left + 720}" y="{y + 28}" class="status">{status}</text>'
    )
outcome = "PASS" if not failed else "FAIL"
svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<style>.title{{font:700 26px sans-serif;fill:#172554}}.subtitle{{font:16px sans-serif;fill:#334155}}.case{{font:16px sans-serif;fill:#0f172a}}.status{{font:700 17px sans-serif;fill:#0f172a}}.footer{{font:700 18px sans-serif;fill:#172554}}</style>
<rect width="100%" height="100%" fill="white"/>
<text x="24" y="35" class="title">DIRT typed scheduler contract matrix</text>
<text x="24" y="61" class="subtitle">Observed precision-double test outcome; every declared contract must PASS.</text>
{''.join(bars)}
<text x="24" y="{height - 28}" class="footer">Result: {outcome} — {len(cases) - len(failed)}/{len(cases)} declared contracts match</text>
</svg>'''
PLOT.parent.mkdir(exist_ok=True)
PLOT.write_text(svg)
print(f"TYPED SCHEDULE: {outcome} — {len(cases) - len(failed)}/{len(cases)} contracts match")
if failed:
    raise SystemExit("contract mismatch: " + ", ".join(failed))
