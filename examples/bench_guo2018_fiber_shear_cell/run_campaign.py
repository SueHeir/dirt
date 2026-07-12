#!/usr/bin/env python3
"""Run the source-derived periodic control-cell campaign.

Each case is materialized, topology-audited, and executed through ``run_case``;
the acceptance validator remains separate and fail-closed on missing histories.
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = ((651, 64), (1735, 64), (3470, 64), (651, 96), (1735, 96), (3470, 96))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    manifest = {"cases": [{"pressure_pa": p, "width_mm": w} for p, w in CASES],
                "ranks_per_case": args.ranks, "validator": "validate.py"}
    (args.output / "campaign_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    for pressure, width in CASES:
        subprocess.run([sys.executable, str(HERE / "run_case.py"),
                        "--pressure-pa", str(pressure), "--width-mm", str(width),
                        "--output", str(args.output / f"p{pressure}_w{width}"),
                        "--ranks", str(args.ranks),
                        *(["--prepare-only"] if args.prepare_only else [])], check=True)


if __name__ == "__main__":
    main()
