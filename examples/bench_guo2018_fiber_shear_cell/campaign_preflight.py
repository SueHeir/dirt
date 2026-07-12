#!/usr/bin/env python3
"""Pre-register a non-equivalent Guo wall-resolution sensitivity campaign.

The paper specifies sphere-built walls but not their sphere diameter or layout.
This is therefore *not* a source-equivalent replication input.  This check
prevents the common failure mode of choosing an unreported wall resolution by
matching the target curves: all candidate walls, normal loads, observables and
the steady-state window must be declared before solver results are accepted.
"""
import argparse
import json
import unittest
from pathlib import Path

from source_geometry_audit import validate_wall_realisation

HERE = Path(__file__).resolve().parent
CAMPAIGN = HERE / "data" / "non_equivalent_sensitivity_campaign.json"
NORMAL_STRESSES = [651, 1735, 3470]
OBSERVABLES = {"shear_stress_pa", "solid_fraction"}
ARTIFACTS = {"solver_input", "cell_history", "observable_summary", "provenance_receipt"}


def load_campaign(path=CAMPAIGN):
    record = json.loads(Path(path).read_text())
    required = {
        "schema", "campaign_label", "source_equivalent", "replication_claim",
        "selection_basis", "reference_observables_consulted_for_mesh_selection",
        "normal_stress_pa", "required_observables", "steady_window",
        "wall_realisations", "required_result_artifacts", "result_selection_rule",
    }
    if set(record) != required or record["schema"] != 1:
        raise ValueError("campaign schema changed")
    if record["source_equivalent"] or record["replication_claim"]:
        raise ValueError("an unreported wall mesh cannot claim source equivalence or replication")
    if record["selection_basis"] != "predeclared_sensitivity":
        raise ValueError("campaign must be a predeclared sensitivity study")
    if record["reference_observables_consulted_for_mesh_selection"]:
        raise ValueError("Fig. 6/7 values must not select a wall mesh")
    if record["normal_stress_pa"] != NORMAL_STRESSES:
        raise ValueError("campaign must retain all three digitized normal-stress cases")
    if set(record["required_observables"]) != OBSERVABLES:
        raise ValueError("campaign must retain both required external observables")
    window = record["steady_window"]
    if set(window) != {"start_shear_strain", "end_shear_strain", "basis"}:
        raise ValueError("steady window must be explicitly recorded")
    if not (0 <= window["start_shear_strain"] < window["end_shear_strain"]):
        raise ValueError("steady window must have a positive, ordered strain interval")
    if len(record["wall_realisations"]) < 3:
        raise ValueError("wall uncertainty requires at least three predeclared resolutions")
    walls = [validate_wall_realisation_file(item) for item in record["wall_realisations"]]
    if len({item["label"] for item in walls}) != len(walls):
        raise ValueError("wall labels must be unique")
    if len({item["diameter_mm"] for item in walls}) != len(walls):
        raise ValueError("wall resolutions must be distinct")
    if set(record["required_result_artifacts"]) != ARTIFACTS:
        raise ValueError("campaign must require solver inputs, histories, summaries, and provenance")
    return record


def validate_wall_realisation_file(item):
    """Use the same no-calibration schema as an external wall manifest."""
    import tempfile
    with tempfile.NamedTemporaryFile("w") as manifest:
        json.dump(item, manifest)
        manifest.flush()
        return validate_wall_realisation(manifest.name)


def require_results(campaign, result_root):
    """Refuse a result claim until every preregistered run is receipted.

    This checks existence only; it deliberately does not score a run against
    Fig. 6/7.  Physics acceptance remains the later two-observable comparison.
    """
    missing = []
    for wall in campaign["wall_realisations"]:
        for stress in campaign["normal_stress_pa"]:
            case = Path(result_root) / wall["label"] / f"p{stress}"
            for artifact in campaign["required_result_artifacts"]:
                if not (case / artifact).is_file():
                    missing.append(str(case / artifact))
    if missing:
        raise ValueError(
            "no solver-backed sensitivity result is available for every preregistered case; "
            f"missing {len(missing)} artifact(s), first: {missing[0]}"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, default=CAMPAIGN)
    parser.add_argument("--require-results", type=Path, metavar="RESULT_ROOT")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        raise SystemExit(not unittest.main(argv=["campaign-preflight"], exit=False).result.wasSuccessful())
    try:
        campaign = load_campaign(args.campaign)
        if args.require_results:
            require_results(campaign, args.require_results)
    except ValueError as error:
        parser.error(str(error))
    print("PASS: non-equivalent wall-resolution campaign is preregistered without Fig. 6/7 selection.")
    print("LIMIT: no solver artifacts were evaluated; this is not a replication result or PASS.")


class CampaignPreflightTests(unittest.TestCase):
    def test_committed_campaign_is_preregistered_and_complete(self):
        campaign = load_campaign()
        self.assertEqual(campaign["normal_stress_pa"], NORMAL_STRESSES)
        self.assertEqual(len(campaign["wall_realisations"]), 3)

    def test_reference_fitting_cannot_be_hidden_in_campaign(self):
        import tempfile
        campaign = load_campaign()
        campaign["reference_observables_consulted_for_mesh_selection"] = ["Fig. 6"]
        with tempfile.NamedTemporaryFile("w") as stream:
            json.dump(campaign, stream); stream.flush()
            with self.assertRaisesRegex(ValueError, "must not select"):
                load_campaign(stream.name)

    def test_partial_results_cannot_be_promoted(self):
        import tempfile
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "no solver-backed"):
                require_results(load_campaign(), directory)


if __name__ == "__main__":
    main()
