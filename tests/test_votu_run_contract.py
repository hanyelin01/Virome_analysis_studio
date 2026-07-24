import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "helpers" / "votu_run_contract.py"


class VotuRunContractTest(unittest.TestCase):
    def command(self, contract: Path, manifest: Path, candidate: Path, min_length: int):
        return [
            sys.executable, str(SCRIPT), "--contract", str(contract),
            "--manifest", str(manifest), "--input", str(candidate),
            "--min-length", str(min_length), "--threads", "8", "--ani", "95",
            "--aligned-fraction", "85", "--read-identity", "95",
            "--read-aligned-percent", "75", "--covered-fraction", "10",
            "--importance-abundance", "5",
        ]

    def test_contract_records_effective_parameters_and_rejects_changes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = root / "manifest.tsv"
            candidate = root / "candidate.fna"
            contract = root / "contract.json"
            manifest.write_text("sample_id\nS01\n", encoding="utf-8")
            candidate.write_text(">v1\nAAAA\n", encoding="utf-8")
            created = subprocess.run(
                self.command(contract, manifest, candidate, 200),
                capture_output=True, text=True,
            )
            self.assertEqual(created.returncode, 0, created.stderr)
            self.assertEqual(json.loads(contract.read_text())["min_length"], 200)
            verified = subprocess.run(
                self.command(contract, manifest, candidate, 200),
                capture_output=True, text=True,
            )
            self.assertEqual(verified.returncode, 0, verified.stderr)
            changed = subprocess.run(
                self.command(contract, manifest, candidate, 1000),
                capture_output=True, text=True,
            )
            self.assertNotEqual(changed.returncode, 0)
            self.assertIn("min_length", changed.stderr)


if __name__ == "__main__":
    unittest.main()
