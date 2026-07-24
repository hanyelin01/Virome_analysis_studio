import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "helpers" / "scan_adapter_reference.py"
REFERENCE = ROOT / "config" / "adapter_sequence_reference.tsv"


class AdapterReferenceScanTest(unittest.TestCase):
    def test_finds_sispa_tag_without_enabling_trimming(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            r1 = root / "r1.fastq"
            r2 = root / "r2.fastq"
            tag = "GACCATCTAGCGACCTCCAC"
            r1.write_text(f"@r1\n{tag}ACGTACGT\n+\n{'I' * 28}\n", encoding="ascii")
            r2.write_text("@r1\nACGTACGT\n+\nIIIIIIII\n", encoding="ascii")
            output = root / "scan.tsv"
            subprocess.run([
                sys.executable, str(SCRIPT), "--reference", str(REFERENCE),
                "--sample", "S01", "--r1", str(r1), "--r2", str(r2),
                "--max-reads", "10", "--output", str(output),
            ], check=True)
            with output.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["sequence_id"], "sispa_k_tag")
            self.assertEqual(rows[0]["matches_at_5prime"], "1")
            self.assertEqual(rows[0]["trimming_action"], "protocol_specific")


if __name__ == "__main__":
    unittest.main()
