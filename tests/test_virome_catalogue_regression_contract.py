from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "virome_catalogue_v2"
HELPERS = ROOT / "scripts" / "helpers"


def table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run(script: str, *args: str) -> None:
    subprocess.run([sys.executable, str(HELPERS / script), *args], check=True)


class ViromeCatalogueRegressionContractTest(unittest.TestCase):
    def test_fixture_preserves_vc_vf_provenance_and_decision_contract(self) -> None:
        contract = json.loads((FIXTURE / "contract.json").read_text(encoding="utf-8"))["expected"]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            catalogue = root / "03_candidate_catalogue"
            run(
                "build_virus_candidate_catalogue.py",
                "--input-fasta", str(FIXTURE / "prepared.fna"),
                "--provenance", str(FIXTURE / "provenance.tsv"),
                "--genomad-fasta", str(FIXTURE / "genomad.fna"),
                "--virsorter-fasta", str(FIXTURE / "virsorter2.fna"),
                "--diamond-virus-hits", str(FIXTURE / "virus_discovery_hits.tsv"),
                "--output-dir", str(catalogue),
            )
            decisions = root / "04_nr_annotation"
            run(
                "resolve_viral_decision.py",
                "--catalogue-fasta", str(catalogue / "VC_catalogue.fna"),
                "--discovery-evidence", str(catalogue / "VC_discovery_evidence.tsv"),
                "--taxonomy", str(FIXTURE / "nr_taxonomy.tsv"),
                "--output-dir", str(decisions),
            )
            final = root / "07_final_catalogue"
            samples = root / "08_sample_results"
            run(
                "build_final_virome_catalogue.py",
                "--checkv-fasta", str(FIXTURE / "checkv.fna"),
                "--checkv-quality", str(FIXTURE / "checkv_quality.tsv"),
                "--decision", str(decisions / "viral_decision.tsv"),
                "--nr-taxonomy", str(FIXTURE / "nr_taxonomy.tsv"),
                "--source-mapping", str(catalogue / "VC_source_mapping.tsv"),
                "--ictv-hits", str(FIXTURE / "ictv_hits.tsv"),
                "--ictv-metadata", str(FIXTURE / "ictv_metadata.tsv"),
                "--manifest", str(FIXTURE / "manifest.tsv"),
                "--catalogue-dir", str(final),
                "--sample-dir", str(samples),
            )
            evidence = table(catalogue / "VC_discovery_evidence.tsv")
            decision_rows = table(decisions / "viral_decision.tsv")
            final_rows = table(final / "VF_catalogue.tsv")
            self.assertEqual(len(evidence), contract["candidate_count"])
            self.assertEqual(sum(row["decision"] == "confirmed_viral" for row in decision_rows), contract["confirmed_viral_count"])
            self.assertEqual(len(final_rows), contract["final_fragment_count"])
            self.assertEqual(final_rows[0]["parent_vc_id"], "VC_0000001")
            self.assertTrue((samples / "S1" / "references" / "VF_0000001.fna").is_file())
            report_root = root / "report"
            shutil.copytree(final.parent, report_root, dirs_exist_ok=True)
            (report_root / "03_candidate_catalogue").mkdir(parents=True, exist_ok=True)
            shutil.copy2(catalogue / "VC_discovery_evidence.tsv", report_root / "03_candidate_catalogue" / "VC_discovery_evidence.tsv")
            (report_root / "04_nr_annotation").mkdir(parents=True, exist_ok=True)
            shutil.copy2(decisions / "viral_decision.tsv", report_root / "04_nr_annotation" / "viral_decision.tsv")
            shutil.copytree(samples, report_root / "08_sample_results", dirs_exist_ok=True)
            run("build_virome_catalogue_report.py", "--output-dir", str(report_root))
            report = (report_root / "reports" / "virome_catalogue_dashboard.html").read_text(encoding="utf-8")
            self.assertIn("全局候选 VC", report)
            self.assertNotIn(contract["forbidden_result_identifier"], " ".join(row["vf_id"] for row in final_rows))
