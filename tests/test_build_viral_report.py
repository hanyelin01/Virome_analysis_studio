"""Smoke test for the offline viral dashboard data contract."""
from __future__ import annotations

import json
import importlib.util
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts" / "helpers" / "build_viral_report.py"


class ViralDashboardTest(unittest.TestCase):
    def test_versioned_ictv_reference_has_expected_attention_family(self) -> None:
        spec = importlib.util.spec_from_file_location("viral_report_builder", BUILDER)
        self.assertIsNotNone(spec)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(module)
        reference_path = ROOT / "config" / "ictv_family_genome_reference.tsv"
        reference = module.load_ictv_reference(reference_path)
        coronavirus = module.reference_for_family("Coronaviridae", reference)
        self.assertEqual(len(reference), 427)
        self.assertEqual(coronavirus["genome_group"], "RNA")
        self.assertEqual(coronavirus["genome_label"], "RNA 病毒")
        row = module.serialise_row({
            "votu_id": "votu_test",
            "representative_sequence_id": "seq_test",
            "representative_length": "4000",
            "covered_bases": "1200",
            "checkv_quality": "Complete",
            "completeness": "98",
            "taxonomy": "f__Coronaviridae;g__DemoVirus",
            "detected": "yes",
        })
        payload = module.sample_payload(
            "sample_test", "SUCCESS", [row],
            module.load_priority_reference(ROOT / "config" / "priority_review_taxa.tsv", set()), reference,
        )
        self.assertEqual(payload["rows"][0]["genome_group"], "RNA")
        self.assertEqual(payload["rows"][0]["covered_fraction"], 30.0)
        self.assertTrue(payload["rows"][0]["priority_family"])

    def test_builds_self_contained_dashboard_and_audit_data(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = root / "sample_manifest.tsv"
            manifest.write_text(
                "sample_id\tread_type\traw_r1\traw_r2\tclean_r1\tclean_r2\tclean_single\tassembly_dir\n"
                "sample_A\tpe\t\t\t\t\t\t/a\n"
                "sample_B\tpe\t\t\t\t\t/b\n",
                encoding="utf-8",
            )
            status = root / "04_sample_votu" / "sample_status.tsv"
            status.parent.mkdir(parents=True)
            status.write_text("sample_id\tstatus\nsample_A\tSUCCESS\nsample_B\tSUCCESS\n", encoding="utf-8")
            header = (
                "sample_id\tvotu_id\trepresentative_sequence_id\trepresentative_length\tmember_count\tcheckv_quality\tmiuvig_quality\tcompleteness\tcontamination\ttaxonomy\tvirus_score\trelative_abundance\tmean_coverage\tcovered_bases\tread_count\tdetected\timportance\n"
            )
            (root / "04_sample_votu" / "sample_A").mkdir()
            (root / "04_sample_votu" / "sample_A" / "votu_summary.tsv").write_text(
                header + "sample_A\tvotu_1\tseq_A\t4000\t2\tComplete\tHigh-quality\t98\t0\tf__DemoViridae;g__DemoVirus\t0.95\t12.5\t30\t1200\t80\tyes\t高优先级\n",
                encoding="utf-8",
            )
            (root / "04_sample_votu" / "sample_B").mkdir()
            (root / "04_sample_votu" / "sample_B" / "votu_summary.tsv").write_text(
                header + "sample_B\tvotu_2\tseq_B\t2100\t1\tMedium-quality\tGenome-fragment\t65\t0\tf__DemoViridae;g__DemoVirus\t0.70\t3.2\t8\t400\t17\tyes\t关注\n",
                encoding="utf-8",
            )
            completed = subprocess.run(
                [sys.executable, str(BUILDER), "--output-dir", str(root), "--manifest", str(manifest), "--overview-rank", "family"],
                capture_output=True,
                text=True,
                check=True,
            )
            self.assertIn("Offline dashboard written", completed.stdout)
            dashboard = root / "reports" / "virome_dashboard.html"
            self.assertTrue(dashboard.is_file())
            text = dashboard.read_text(encoding="utf-8")
            self.assertIn("病毒筛查报告中心", text)
            self.assertIn("DemoViridae", text)
            self.assertIn("mapped reads", text)
            self.assertIn("navigateToSampleV4", text)
            self.assertIn("backdrop-filter:blur", text)
            if shutil.which("node"):
                script = re.search(r"<script>(.*)</script>", text, flags=re.S)
                self.assertIsNotNone(script)
                js_path = root / "dashboard.js"
                js_path.write_text(script.group(1), encoding="utf-8")
                subprocess.run(["node", "--check", str(js_path)], check=True)
            self.assertTrue((root / "reports" / "data" / "report_metadata.json").is_file())
            payload = json.loads((root / "reports" / "data" / "samples" / "sample_A.json").read_text(encoding="utf-8"))
            self.assertEqual(payload["summary"]["detected_local_votu_count"], 1)
            self.assertEqual(payload["rows"][0]["read_count"], 80)
            self.assertIn("virome_dashboard.html", (root / "reports" / "samples" / "sample_A.html").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
