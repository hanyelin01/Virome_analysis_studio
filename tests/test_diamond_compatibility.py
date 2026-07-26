from __future__ import annotations

import csv
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "scripts" / "helpers" / "summarize_diamond_taxonomy.py"


class DiamondCompatibilityTest(unittest.TestCase):
    def test_current_diamond_fields_preserve_taxonkit_inputs(self) -> None:
        """The workflow requires a current DIAMOND with full taxonomy fields."""
        for script in ("05c_run_diamond_virus_discovery.sh", "10_run_diamond_taxonomy.sh"):
            source = (ROOT / "scripts" / script).read_text(encoding="utf-8")
            fields = re.search(r"^fields=\(([^)]*)\)$", source, flags=re.MULTILINE)
            self.assertIsNotNone(fields)
            self.assertIn("slineages", fields.group(1))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            hits = root / "hits.tsv"
            hits.write_text(
                "Q1\t200\t1\t180\t99\t180\t1e-30\t120\t2\t181\tref1\t10239\tViruses\tViruses\tViruses; Exampleviridae\n",
                encoding="utf-8",
            )
            lca = root / "lca.tsv"
            lca.write_text("Q1\t10239\t10239\n", encoding="utf-8")
            lineage = root / "lineage.tsv"
            lineage.write_text("Q1\t10239\t10239\tViruses\tViruses\tsuperkingdom\n", encoding="utf-8")
            output = root / "summary.tsv"
            subprocess.run(
                [sys.executable, str(SUMMARY), "--hits", str(hits), "--lca", str(lca), "--lineage", str(lineage), "--output", str(output)],
                check=True,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            with output.open(encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["best_staxids"], "10239")
            self.assertEqual(row["best_scientific_names"], "Viruses")
            self.assertEqual(row["best_lineages"], "Viruses; Exampleviridae")
            self.assertEqual(row["lca_lineage"], "Viruses")


if __name__ == "__main__":
    unittest.main()
