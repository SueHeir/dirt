#!/usr/bin/env python3
"""Execute the six predeclared Guo shear-cell solver cases under MPI.

This is deliberately an execution/provenance layer, not a result generator.
Each case is materialized by ``run_case.py`` and receives an immutable
case_manifest.json before DIRT starts.  ``validate.py`` remains the only
component allowed to declare the two-observable replication gate passed.
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
                        *( ["--prepare-only"] if args.prepare_only else [])], check=True)


if __name__ == "__main__":
    main()
