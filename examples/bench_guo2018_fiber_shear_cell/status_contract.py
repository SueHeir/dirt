#!/usr/bin/env python3
"""Fail closed on an unsupported Guo shear-cell replication claim."""
import argparse
import json
import tempfile
import unittest
from pathlib import Path

from reconstruction_readiness import load_ledger
from source_geometry_audit import load_contract

HERE = Path(__file__).resolve().parent
STATUS = HERE / "data" / "replication_status.json"


def load_status(path=STATUS):
    record = json.loads(Path(path).read_text())
    required = {"schema", "case", "doi", "classification", "source_equivalent_replication",
                "solver_campaign_present", "comparison_to_published_observables_present",
                "blocking_source_inputs", "reason"}
    if set(record) != required or record["schema"] != 1:
        raise ValueError("replication-status schema changed")
    if record["doi"] != "10.1002/aic.16397":
        raise ValueError("status is not bound to the Guo primary DOI")
    if record["classification"] != "source-audit-only":
        raise ValueError("this incomplete case cannot be classified as a replication")
    return record


def verify_status(path=STATUS):
    status, ledger, geometry = load_status(path), load_ledger(), load_contract()
    unresolved = set(ledger["unresolved_required_boundary_inputs"])
    # The Figure-2 scale entry is an evidentiary limitation, not an additional
    # wall input.  The diameter and layout remain the two inputs that block a
    # source-equivalent solver reconstruction.
    missing_geometry = {"wall_sphere_diameter_mm", "wall_sphere_layout"}
    if set(status["blocking_source_inputs"]) != missing_geometry or not missing_geometry <= unresolved:
        raise ValueError("status blockers must exactly match the independently recorded missing wall geometry")
    if status["source_equivalent_replication"]:
        raise ValueError("unreported wall geometry forbids a source-equivalent replication claim")
    if status["solver_campaign_present"] or status["comparison_to_published_observables_present"]:
        raise ValueError("a source-audit-only case cannot contain solver or published-observable result claims")
    return status


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        raise SystemExit(not unittest.main(argv=["status-contract"], exit=False).result.wasSuccessful())
    if not args.verify:
        parser.error("--verify is required")
    status = verify_status()
    print("STATUS OK: source-audit-only; no source-equivalent solver replication is claimed.")
    print("BLOCKED BY PRIMARY SOURCE:", ", ".join(status["blocking_source_inputs"]))


class StatusContractTests(unittest.TestCase):
    def test_committed_status_tracks_the_independent_ledgers(self):
        self.assertEqual(verify_status()["classification"], "source-audit-only")

    def test_solver_claim_cannot_be_added_without_resolving_source_geometry(self):
        status = load_status(); status["solver_campaign_present"] = True
        with tempfile.NamedTemporaryFile("w") as stream:
            json.dump(status, stream); stream.flush()
            with self.assertRaisesRegex(ValueError, "cannot contain solver"):
                verify_status(stream.name)

    def test_blocker_substitution_is_rejected(self):
        status = load_status(); status["blocking_source_inputs"] = ["upper_wall_mass_kg"]
        with tempfile.NamedTemporaryFile("w") as stream:
            json.dump(status, stream); stream.flush()
            with self.assertRaisesRegex(ValueError, "exactly match"):
                verify_status(stream.name)


if __name__ == "__main__":
    main()
