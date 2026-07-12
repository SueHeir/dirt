#!/usr/bin/env python3
"""Audit the published apparatus contract before allowing a Guo-cell run.

The DOI metadata is an independent source: its abstract identifies the
apparatus as a *Schulze ring shear tester*.  A translating planar base in a
cylindrical cup is a different experiment, so it must never produce a history
that `validate.py` can present as a Guo et al. replication.
"""
import argparse
import json
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
DOI = "10.1002/aic.16397"
CROSSREF = f"https://api.crossref.org/works/{DOI}"


def published_apparatus() -> dict:
    """Fetch the publisher-indexed record, rather than trusting local prose."""
    with urllib.request.urlopen(CROSSREF, timeout=30) as response:
        record = json.load(response)["message"]
    title = " ".join(record.get("title", []))
    abstract = record.get("abstract", "")
    if DOI not in record.get("DOI", "").lower() or "Schulze ring shear tester" not in abstract:
        raise ValueError("Crossref record does not establish the reported ring-shear apparatus")
    return {"title": title, "doi": record["DOI"], "apparatus": "Schulze ring shear tester"}


def current_protocol_is_not_ring_shear(config: Path) -> None:
    """Document the incompatibility in the current DIRT wall capability.

    This is intentionally a *blocker check*, not a numerical acceptance test.
    It rejects the old plane-translation/cylindrical-cup surrogate even if it
    happens to yield values near a digitized figure.
    """
    text = config.read_text()
    has_translating_plane = 'name = "lower_drive"' in text and 'velocity = [0.0, 0.0, 0.0]' in text
    has_cylindrical_cup = 'name = "cylindrical_sidewall"' in text
    if not (has_translating_plane and has_cylindrical_cup):
        raise ValueError("expected the documented plane-driven cylindrical-cup surrogate")


def require_runnable_ring_protocol(config: Path) -> None:
    published_apparatus()
    current_protocol_is_not_ring_shear(config)
    raise RuntimeError(
        "BLOCKED: DIRT has no configured annular, rotational ring-shear drive. "
        "Do not run or compare the plane-driven cylindrical-cup surrogate to Guo et al."
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=HERE / "config.toml")
    parser.add_argument("--verify-doi", action="store_true", help="query Crossref and print the external apparatus evidence")
    parser.add_argument("--require-runnable", action="store_true", help="raise the current protocol blocker")
    args = parser.parse_args()
    if not args.verify_doi and not args.require_runnable:
        parser.error("choose --verify-doi and/or --require-runnable")
    if args.verify_doi:
        record = published_apparatus()
        print(f"EXTERNAL EVIDENCE: {record['doi']} — {record['title']}; apparatus: {record['apparatus']}")
    if args.require_runnable:
        require_runnable_ring_protocol(args.config)


if __name__ == "__main__":
    main()
