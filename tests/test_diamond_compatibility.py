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
CATALOGUE = ROOT / "scripts" / "helpers" / "build_virus_candidate_catalogue.py"


class DiamondCompatibilityTest(unittest.TestCase):
    def test_current_diamond_fields_preserve_taxonkit_inputs(self) -> None:
        """The workflow requires a current DIAMOND with full taxonomy fields."""
        for script in ("05c_run_diamond_virus_discovery.sh", "10_run_diamond_taxonomy.sh"):
            source = (ROOT / "scripts" / script).read_text(encoding="utf-8")
            fields = re.search(r"^fields=\(([^)]*)\)$", source, flags=re.MULTILINE)
            self.assertIsNotNone(fields)
            self.assertIn("slineages", fields.group(1))
        discovery_source = (ROOT / "scripts" / "05c_run_diamond_virus_discovery.sh").read_text(encoding="utf-8")
        self.assertIn('COMPLETE="$OUT/parameters.env"', discovery_source)
        runner_source = (ROOT / "scripts" / "run_virome_catalogue.sh").read_text(encoding="utf-8")
        self.assertIn("launch_megan_background", runner_source)
        self.assertIn("setsid bash -c", runner_source)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            hits = root / "hits.tsv"
            long_lineage = "Viruses; Exampleviridae; " + "X" * 5000
            hits.write_text(
                f"Q1\t200\t1\t180\t99\t180\t1e-30\t120\t2\t181\tref1\t10239\tViruses\tViruses\t{long_lineage}\n",
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
            self.assertTrue(row["best_lineages"].endswith("[truncated; see raw DIAMOND TSV]"))
            self.assertEqual(row["lca_lineage"], "Viruses")

    def test_catalogue_streams_a_long_raw_lineage_field(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prepared = root / "prepared.fna"
            prepared.write_text(">S1__contig\nACGTACGT\n", encoding="utf-8")
            provenance = root / "provenance.tsv"
            provenance.write_text("sequence_id\tsample_id\toriginal_contig_id\tlength\nS1__contig\tS1\tcontig\t8\n", encoding="utf-8")
            genomad = root / "genomad.fna"; genomad.write_text("", encoding="utf-8")
            virsorter = root / "virsorter.fna"; virsorter.write_text("", encoding="utf-8")
            hits = root / "hits.tsv"
            hits.write_text("S1__contig\t8\t1\t8\t99\t8\t1e-20\t80\tref\t10239\tViruses\t" + "X" * 150000 + "\n", encoding="utf-8")
            output = root / "catalogue"
            subprocess.run(
                [sys.executable, str(CATALOGUE), "--input-fasta", str(prepared), "--provenance", str(provenance), "--genomad-fasta", str(genomad), "--virsorter-fasta", str(virsorter), "--diamond-virus-hits", str(hits), "--output-dir", str(output)],
                check=True,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            with (output / "VC_discovery_evidence.tsv").open(encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["DIAMOND_NR_virus"], "yes")


if __name__ == "__main__":
    unittest.main()
