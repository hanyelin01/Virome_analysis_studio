import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "helpers" / "split_viral_candidates_by_sample.py"


def load_module():
    spec = importlib.util.spec_from_file_location("split_candidates", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(module)
    return module


class SplitViralCandidatesTest(unittest.TestCase):
    def test_exact_id_wins_over_suffix_fallback(self):
        module = load_module()
        table = {
            "sample__contig_1": {"sample_id": "A"},
            "sample__contig": {"sample_id": "B"},
        }
        sequence_id, row, method = module.resolve_provenance(table, "sample__contig_1")
        self.assertEqual(sequence_id, "sample__contig_1")
        self.assertEqual(row["sample_id"], "A")
        self.assertEqual(method, "exact")

    def test_checkv_provirus_child_restores_parent_and_audit_fields(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidate = root / "candidate.fna"
            provenance = root / "provenance.tsv"
            quality = root / "quality.tsv"
            taxonomy = root / "taxonomy.tsv"
            manifest = root / "manifest.tsv"
            output = root / "out"
            parent = "sample_A__k79_43359"
            child = parent + "_1"
            candidate.write_text(f">{child}\n{'A' * 1200}\n", encoding="utf-8")
            provenance.write_text(
                "sequence_id\tsample_id\toriginal_contig_id\tlength\n"
                f"{parent}\tsample_A\tk79_43359\t1500\n",
                encoding="utf-8",
            )
            quality.write_text(
                "contig_id\tprovirus\tcheckv_quality\tmiuvig_quality\tcompleteness\tcontamination\n"
                f"{parent}\tYes\tLow-quality\tGenome-fragment\t1.01\t1.0\n",
                encoding="utf-8",
            )
            taxonomy.write_text(
                "seq_name\ttaxonomy\tvirus_score\n"
                f"{parent}\tViruses;Caudoviricetes\t0.93\n",
                encoding="utf-8",
            )
            manifest.write_text(
                "sample_id\tread_type\traw_r1\traw_r2\tclean_r1\tclean_r2\tclean_single\tassembly_dir\n"
                "sample_A\tpe\t\t\t/r1\t/r2\t\t/a\n",
                encoding="utf-8",
            )
            completed = subprocess.run([
                sys.executable, str(SCRIPT), "--input", str(candidate),
                "--provenance", str(provenance), "--quality", str(quality),
                "--taxonomy", str(taxonomy), "--manifest", str(manifest),
                "--output-root", str(output), "--min-length", "1000",
            ], capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            metadata = output / "sample_A" / "01_candidates" / "candidate_metadata.tsv"
            with metadata.open(encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["sequence_id"], child)
            self.assertEqual(row["provenance_sequence_id"], parent)
            self.assertEqual(row["provenance_match"], "checkv_provirus_suffix")
            self.assertEqual(row["taxonomy"], "Viruses;Caudoviricetes")
            completion = json.loads(
                (output / "split_complete.json").read_text(encoding="utf-8")
            )
            self.assertEqual(completion["unassigned_candidate_count"], 0)
            self.assertEqual(completion["kept_candidate_count"], 1)
            self.assertEqual(
                (output / "unassigned_candidates.tsv").read_text(encoding="utf-8").count("\n"),
                1,
            )


if __name__ == "__main__":
    unittest.main()
