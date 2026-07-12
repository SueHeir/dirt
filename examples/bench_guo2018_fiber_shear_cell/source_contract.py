#!/usr/bin/env python3
"""Guard this unreceipted draft against being called a Guo replication.

Crossref can verify the DOI/title, but the primary article is not available in
this checkout.  Therefore no numerical-control-cell detail in ``config.toml``
is source-authenticated.  The historical draft must not execute or be used to
make a geometry claim until ``evidence_contract.py`` accepts a primary-source
receipt.  This is not a physics validation or a substitute for the two
published-observable comparisons required by the goal.
"""
import argparse
import json
import unittest
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
    """Always reject the historical draft before it reaches a solver."""
    del config
    raise RuntimeError(
        "BLOCKED: no primary-source receipt authenticates the control-cell geometry or wall construction; "
        "the historical input is not runnable as a Guo replication"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=HERE / "config.toml")
    parser.add_argument("--verify-doi", action="store_true", help="query Crossref for bibliographic identity")
    parser.add_argument("--require-runnable", action="store_true", help="raise unless the configured numerical cell matches the source protocol")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(SourceContractTests)
        raise SystemExit(not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful())
    if not args.verify_doi and not args.require_runnable:
        parser.error("choose --verify-doi and/or --require-runnable")
    if args.verify_doi:
        record = published_record()
        print(f"EXTERNAL BIBLIOGRAPHY: {record['doi']} — {record['title']}")
        print("No primary-source method detail is authenticated in this checkout.")
    if args.require_runnable:
        require_published_control_cell(args.config)


class SourceContractTests(unittest.TestCase):
    def test_rejects_unreceipted_historical_geometry(self):
        with self.assertRaisesRegex(RuntimeError, "no primary-source receipt"):
            require_published_control_cell(HERE / "config.toml")


if __name__ == "__main__":
    main()
