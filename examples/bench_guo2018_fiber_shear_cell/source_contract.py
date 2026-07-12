#!/usr/bin/env python3
"""Check that a DIRT case matches Guo et al.'s *numerical* control cell.

Guo et al. used a Schulze ring shear tester experimentally, but deliberately
did **not** simulate its annulus.  Their numerical comparison (Computational
Set-up, pp. 5--7, doi:10.1002/aic.16397) uses a small plane-driven cell with
periodic x/z boundaries, a vertically mobile loaded lid, and 4/2-mm lid/base
blades.  Conflating the apparatus with this published DEM reduction was the
previous branch's central source error.

This checker is intentionally a protocol gate, not a result gate.  It permits
solver execution only after the configured numerical boundary conditions are
source-compatible; it cannot turn a topology audit or a solver history into a
successful Fig. 6/7 replication.
"""
import argparse
import json
import tomllib
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
DOI = "10.1002/aic.16397"
CROSSREF = f"https://api.crossref.org/works/{DOI}"


def published_record() -> dict:
    """Fetch bibliographic identity only; geometry comes from the primary paper."""
    with urllib.request.urlopen(CROSSREF, timeout=30) as response:
        record = json.load(response)["message"]
    title = " ".join(record.get("title", []))
    if record.get("DOI", "").lower() != DOI or "Flexible Fiber Flows" not in title:
        raise ValueError("Crossref record does not identify Guo et al. doi:10.1002/aic.16397")
    return {"title": title, "doi": record["DOI"]}


def require_published_control_cell(config: Path) -> None:
    """Reject only deviations from the paper's stated DEM reduction.

    The current cylinder/fixed-boundary input is neither the physical SRST nor
    the periodic planar control cell used for Fig. 6/7.  Do not silently treat
    its output as an experimental comparison.
    """
    with config.open("rb") as stream:
        input_config = tomllib.load(stream)
    domain = input_config.get("domain", {})
    walls = input_config.get("wall", [])
    failures = []
    if domain.get("boundary_x") != "periodic" or domain.get("boundary_z") != "periodic":
        failures.append("the published numerical cell requires periodic x and z boundaries")
    if any(wall.get("type") == "cylinder" for wall in walls):
        failures.append("the published numerical cell has no cylindrical sidewall")
    names = {wall.get("name") for wall in walls}
    if not {"upper_blades", "lower_blades"} <= names:
        failures.append("the published numerical cell requires 4-mm upper and 2-mm lower blade arrays")
    if failures:
        raise RuntimeError("BLOCKED: current DIRT input is not the published periodic control cell: " + "; ".join(failures))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=HERE / "config.toml")
    parser.add_argument("--verify-doi", action="store_true", help="query Crossref for bibliographic identity")
    parser.add_argument("--require-runnable", action="store_true", help="raise unless the configured numerical cell matches the source protocol")
    args = parser.parse_args()
    if not args.verify_doi and not args.require_runnable:
        parser.error("choose --verify-doi and/or --require-runnable")
    if args.verify_doi:
        record = published_record()
        print(f"EXTERNAL BIBLIOGRAPHY: {record['doi']} — {record['title']}")
        print("PRIMARY-SOURCE METHOD: pp. 5--7 use periodic x/z planar control cell, loaded mobile lid, and 4/2-mm blade arrays.")
    if args.require_runnable:
        require_published_control_cell(args.config)


if __name__ == "__main__":
    main()
