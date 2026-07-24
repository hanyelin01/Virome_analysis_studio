import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "helpers" / "adapter_evidence.py"
CATALOG = ROOT / "config" / "adapter_catalog.tsv"


def load_module():
    spec = importlib.util.spec_from_file_location("adapter_evidence", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(module)
    return module


class AdapterEvidenceTest(unittest.TestCase):
    def test_shipped_catalogue_is_valid(self):
        rows = load_module().load_catalog(CATALOG)
        self.assertEqual(rows[0]["profile_id"], "auto")
        self.assertTrue(any(
            row["profile_id"] == "illumina_truseq_nebnext" for row in rows
        ))

    def test_exact_sequence_matches_catalogue(self):
        module = load_module()
        matched, judgement = module.sequence_match(
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            module.load_catalog(CATALOG),
        )
        self.assertIn("illumina_truseq_nebnext", matched)
        self.assertEqual(judgement, "exact_catalogue_match")

    def test_fastp_json_becomes_auditable_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            report = root / "sample.fastp.json"
            output = root / "sample.adapter_evidence.tsv"
            report.write_text(json.dumps({
                "summary": {"before_filtering": {"total_reads": 1000}},
                "adapter_cutting": {
                    "adapter_trimmed_reads": 125,
                    "read1_adapter_sequence": "CTGTCTCTTATACACATCT",
                    "read2_adapter_sequence": "CTGTCTCTTATACACATCT",
                },
            }), encoding="utf-8")
            subprocess.run([
                sys.executable, str(SCRIPT), "summarize",
                "--catalog", str(CATALOG), "--profile", "auto",
                "--fastp-json", str(report), "--sample", "S01", "--output", str(output),
            ], check=True)
            with output.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["source_judgement"], "exact_catalogue_match")
            self.assertEqual(rows[0]["trimmed_read_fraction"], "0.125000")
            self.assertEqual(rows[0]["source_level"], "vendor")
            self.assertTrue(rows[0]["catalog_sha256"])

    def test_manual_profile_is_not_claimed_as_auto_detection(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            report = root / "sample.fastp.json"
            output = root / "sample.adapter_evidence.tsv"
            report.write_text(json.dumps({
                "summary": {"before_filtering": {"total_reads": 10}},
                "adapter_cutting": {
                    "adapter_trimmed_reads": 1,
                    "read1_adapter_sequence": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
                    "read2_adapter_sequence": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
                },
            }), encoding="utf-8")
            subprocess.run([
                sys.executable, str(SCRIPT), "summarize",
                "--catalog", str(CATALOG),
                "--profile", "illumina_truseq_nebnext",
                "--fastp-json", str(report), "--sample", "S02", "--output", str(output),
            ], check=True)
            text = output.read_text(encoding="utf-8")
            self.assertIn("reported_sequence_consistent_with_selected_profile", text)
            self.assertIn("configured_fallback_plus_fastp_auto", text)


if __name__ == "__main__":
    unittest.main()
