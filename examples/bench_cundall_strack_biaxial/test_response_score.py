import hashlib
import json
import pathlib
import tempfile
import unittest

from response_score import score


class ResponseScoreTests(unittest.TestCase):
    def write(self, root, name, body):
        path = root / name
        path.write_text(body)
        return {"path": name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}

    def package(self, root, dirt, reference, criteria):
        manifest = {"dirt_response": dirt, "reference_response": reference, "criteria": criteria}
        path = root / "score.json"
        path.write_text(json.dumps(manifest))
        return path

    def complete_case(self, root, reference_states="0,0.10,0.00,0.00\n1,0.20,0.01,0.02\n"):
        heading = "state,stress_ratio,volumetric_strain,fabric_anisotropy\n"
        dirt = self.write(root, "dirt.csv", heading + "0,0.11,0.01,0.01\n1,0.21,0.02,0.03\n")
        reference = self.write(root, "reference.csv", heading + reference_states)
        criteria = self.write(root, "criteria.json", json.dumps({
            name: {"max_abs_error": 0.02, "reference": "independent archived protocol"}
            for name in ("stress_ratio", "volumetric_strain", "fabric_anisotropy")
        }))
        return self.package(root, dirt, reference, criteria)

    def test_scores_all_goal_response_series_at_registered_states(self):
        with tempfile.TemporaryDirectory() as directory:
            report = score(self.complete_case(pathlib.Path(directory)))
        self.assertTrue(report["passed"])
        self.assertEqual(set(report["errors"]), {"stress_ratio", "volumetric_strain", "fabric_anisotropy"})

    def test_rejects_unregistered_reference_states_instead_of_interpolating(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(RuntimeError, "states must be identical"):
                score(self.complete_case(pathlib.Path(directory), "0.5,0.10,0.00,0.00\n1,0.20,0.01,0.02\n"))

    def test_rejects_uncited_or_nonpositive_tolerance(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            package = self.complete_case(root)
            criteria = root / "criteria.json"
            payload = json.loads(criteria.read_text())
            payload["stress_ratio"] = {"max_abs_error": 0, "reference": ""}
            criteria.write_text(json.dumps(payload))
            manifest = json.loads(package.read_text())
            manifest["criteria"]["sha256"] = hashlib.sha256(criteria.read_bytes()).hexdigest()
            package.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(RuntimeError, "positive bound and external reference"):
                score(package)


if __name__ == "__main__":
    unittest.main()
