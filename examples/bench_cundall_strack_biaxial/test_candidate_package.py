import csv
import hashlib
import json
import pathlib
import tempfile
import unittest

from candidate_package import admit
from candidate_admission import REQUIRED_PROTOCOL_FIELDS


class CandidatePackageTests(unittest.TestCase):
    def write_artifact(self, root, name, contents):
        path = root / name
        path.write_text(contents)
        return {"path": name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}

    def package(self, root, ledger, response):
        input_artifact = self.write_artifact(root, "independent.in", "independent solver input\n")
        package = {
            "candidate": "independent walled quasi-2-D trajectory",
            "ledger": ledger,
            "input": input_artifact,
            "response": response,
        }
        path = root / "candidate.json"
        path.write_text(json.dumps(package))
        return path

    def test_complete_hashed_candidate_with_all_response_series_is_admissible(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            ledger_rows = [
                {"candidate": "independent", "field": field, "common": "yes",
                 "evidence": "archived independent input", "limitation": "none"}
                for field in REQUIRED_PROTOCOL_FIELDS
            ]
            ledger_path = root / "ledger.csv"
            with ledger_path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=ledger_rows[0])
                writer.writeheader(); writer.writerows(ledger_rows)
            ledger = {"path": "ledger.csv", "sha256": hashlib.sha256(ledger_path.read_bytes()).hexdigest()}
            response = self.write_artifact(root, "response.csv", "state,stress_ratio,volumetric_strain,fabric_anisotropy\n0,0.1,0,0\n1,0.2,0.01,0.02\n")
            self.assertTrue(admit(self.package(root, ledger, response)).eligible)

    def test_hash_tampering_is_rejected_before_candidate_is_admitted(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            ledger = self.write_artifact(root, "ledger.csv", "candidate,field,common,evidence,limitation\n")
            response = self.write_artifact(root, "response.csv", "state,stress_ratio,volumetric_strain,fabric_anisotropy\n0,0,0,0\n")
            package = self.package(root, ledger, response)
            (root / "response.csv").write_text("tampered\n")
            with self.assertRaisesRegex(RuntimeError, "response artifact hash mismatch"):
                admit(package)


if __name__ == "__main__":
    unittest.main()
