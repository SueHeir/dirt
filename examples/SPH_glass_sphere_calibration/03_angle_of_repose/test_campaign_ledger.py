#!/usr/bin/env python3
"""Independent regression checks for hash-bound campaign-ledger admission."""

import importlib.util
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
SPEC = importlib.util.spec_from_file_location("repose_sweep", HERE / "sweep.py")
repose_sweep = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(repose_sweep)


class CampaignLedgerTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.old_sweep_dir = repose_sweep.SWEEP_DIR
        self.old_ledger = repose_sweep.CAMPAIGN_LEDGER
        repose_sweep.SWEEP_DIR = str(self.root / "sweep")
        repose_sweep.CAMPAIGN_LEDGER = str(self.root / "campaign_ledger.json")

    def tearDown(self):
        repose_sweep.SWEEP_DIR = self.old_sweep_dir
        repose_sweep.CAMPAIGN_LEDGER = self.old_ledger
        self.tempdir.cleanup()

    def _valid_ledger(self):
        entries = []
        for mu_r, rep in repose_sweep._expected_cases():
            case = Path(repose_sweep.case_dir(mu_r, rep))
            data = case / "data"
            data.mkdir(parents=True)
            seed = 1000 + rep + round(mu_r * 100)
            (case / "config.toml").write_text(f"seed = {seed}\n", encoding="utf-8")
            (data / "repose_results.csv").write_text("x,y,z,radius\n", encoding="utf-8")
            (data / "repose_qualification.json").write_text("{}\n", encoding="utf-8")
            (data / "repose_prelift.csv").write_text("x,y,z\n", encoding="utf-8")
            # Construct the receipt independently of sweep.py's writer so this
            # test exercises the verifier rather than round-tripping one helper.
            def digest(path):
                return hashlib.sha256(path.read_bytes()).hexdigest()
            entries.append({
                "mu_r": mu_r,
                "rep": rep,
                "seed": seed,
                "config_sha256": digest(case / "config.toml"),
                "results_sha256": digest(data / "repose_results.csv"),
                "qualification_sha256": digest(data / "repose_qualification.json"),
                "prelift_sha256": digest(data / "repose_prelift.csv"),
            })
        Path(repose_sweep.CAMPAIGN_LEDGER).write_text(
            json.dumps({"schema": 1, "cases": entries}), encoding="utf-8")

    def test_complete_hash_bound_ledger_is_accepted(self):
        self._valid_ledger()
        self.assertTrue(repose_sweep._campaign_ledger_ok())

    def test_seed_edit_is_rejected_even_if_digest_is_rewritten(self):
        self._valid_ledger()
        ledger_path = Path(repose_sweep.CAMPAIGN_LEDGER)
        ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
        case = Path(repose_sweep.case_dir(0.0, 0))
        config = case / "config.toml"
        config.write_text("seed = 999999\n", encoding="utf-8")
        ledger["cases"][0]["config_sha256"] = hashlib.sha256(config.read_bytes()).hexdigest()
        ledger_path.write_text(json.dumps(ledger), encoding="utf-8")
        self.assertFalse(repose_sweep._campaign_ledger_ok())


if __name__ == "__main__":
    unittest.main()
