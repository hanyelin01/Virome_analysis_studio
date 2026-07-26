from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def write(path: Path, content: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run(script: str, *args: str) -> None:
    subprocess.run([sys.executable, str(ROOT / "scripts/helpers" / script), *args], check=True)


def test_global_catalogue_and_final_fragment_distribution(tmp_path: Path) -> None:
    prepared = write(tmp_path / "prepared.fna", ">S1__a\nACGTACGT\n>S2__b\nTTTTCCCC\n")
    provenance = write(tmp_path / "provenance.tsv", "sequence_id\tsample_id\toriginal_contig_id\tlength\nS1__a\tS1\ta\t8\nS2__b\tS2\tb\t8\n")
    genomad = write(tmp_path / "genomad.fna", ">S1__a\nACGTACGT\n")
    virsorter = write(tmp_path / "virsorter.fna", ">S1__a||full\nACGTACGT\n")
    virus_hits = write(tmp_path / "virus_hits.tsv", "qseqid\tqlen\tqstart\tqend\tpident\tlength\tevalue\tbitscore\tsseqid\tstaxids\tsscinames\tslineages\nS2__b\t8\t1\t8\t99\t8\t1e-20\t80\tx\t10239\tvirus\tViruses\n")
    catalogue = tmp_path / "catalogue"
    run("build_virus_candidate_catalogue.py", "--input-fasta", str(prepared), "--provenance", str(provenance), "--genomad-fasta", str(genomad), "--virsorter-fasta", str(virsorter), "--diamond-virus-hits", str(virus_hits), "--output-dir", str(catalogue))
    evidence = table(catalogue / "VC_discovery_evidence.tsv")
    assert len(evidence) == 2
    dual = next(row for row in evidence if row["supporting_method_count"] == "2")
    assert dual["discovery_pattern"] == "VirSorter2 + geNomad"

    taxonomy = write(tmp_path / "taxonomy.tsv", "query_id\tquery_taxids\tlca_taxid\tlca_lineage\tlca_name\tlca_rank\tbest_subject_id\tbest_evalue\tbest_bitscore\tbest_staxids\tbest_scientific_names\tbest_lineages\n" + f"{dual['vc_id']}\t10239\t10239\tViruses;f__Coronaviridae\tCoronaviridae\tfamily\tref\t1e-30\t100\t10239\tvirus\tViruses; Coronaviridae\n")
    decision_dir = tmp_path / "decision"
    run("resolve_viral_decision.py", "--catalogue-fasta", str(catalogue / "VC_catalogue.fna"), "--discovery-evidence", str(catalogue / "VC_discovery_evidence.tsv"), "--taxonomy", str(taxonomy), "--output-dir", str(decision_dir))
    decisions = table(decision_dir / "viral_decision.tsv")
    assert next(row for row in decisions if row["vc_id"] == dual["vc_id"])["decision"] == "confirmed_viral"

    checkv_fasta = write(tmp_path / "checkv.fna", f">{dual['vc_id']}_1\nACGTACGT\n")
    quality = write(tmp_path / "quality.tsv", f"contig_id\tcheckv_quality\tcompleteness\tcontamination\n{dual['vc_id']}\tHigh-quality\t90\t0\n")
    # Historical DIAMOND outfmt-6 files do not contain a header; the final
    # catalogue must remain able to read them after --header was added.
    ictv_hits = write(tmp_path / "ictv.tsv", f"{dual['vc_id']}_1\t8\t1\t8\t99\t8\t1e-20\t90\tREF1\n")
    ictv_meta = write(tmp_path / "ictv_meta.tsv", "reference_id\tfamily\tgenus\tspecies\tbaltimore_group\nREF1\tCoronaviridae\tBetacoronavirus\tSevere acute respiratory syndrome-related coronavirus\tIV\n")
    manifest = write(tmp_path / "manifest.tsv", "sample_id\nS1\nS2\n")
    final_dir, samples = tmp_path / "final", tmp_path / "samples"
    run("build_final_virome_catalogue.py", "--checkv-fasta", str(checkv_fasta), "--checkv-quality", str(quality), "--decision", str(decision_dir / "viral_decision.tsv"), "--nr-taxonomy", str(taxonomy), "--source-mapping", str(catalogue / "VC_source_mapping.tsv"), "--ictv-hits", str(ictv_hits), "--ictv-metadata", str(ictv_meta), "--manifest", str(manifest), "--catalogue-dir", str(final_dir), "--sample-dir", str(samples))
    final = table(final_dir / "VF_catalogue.tsv")
    assert final[0]["parent_vc_id"] == dual["vc_id"]
    assert final[0]["ictv_species"].startswith("Severe acute")
    assert (samples / "S1" / "references" / "VF_0000001.fna").is_file()


def test_report_labels_vc_and_vf_not_votu(tmp_path: Path) -> None:
    write(tmp_path / "03_candidate_catalogue/VC_discovery_evidence.tsv", "vc_id\tdiscovery_pattern\nVC_0000001\tgeNomad + VirSorter2\n")
    write(tmp_path / "04_nr_annotation/viral_decision.tsv", "vc_id\tdecision\nVC_0000001\tconfirmed_viral\n")
    write(tmp_path / "07_final_catalogue/VF_catalogue.tsv", "vf_id\tparent_vc_id\tdecision\tcheckv_quality\tnr_family\tictv_species\tbaltimore_group\tictv_pident\nVF_1\tVC_0000001\tconfirmed_viral\tHigh-quality\tCoronaviridae\tExample virus\tIV\t99\n")
    write(tmp_path / "08_sample_results/sample_fragment_presence.tsv", "sample_id\tvf_id\tparent_vc_id\nS1\tVF_1\tVC_0000001\n")
    run("build_virome_catalogue_report.py", "--output-dir", str(tmp_path))
    report = (tmp_path / "reports/virome_catalogue_dashboard.html").read_text(encoding="utf-8")
    assert "全局候选 VC" in report and "最终片段 VF" in report
    assert "vOTU" in report  # Explicit boundary statement, not a result identifier.


class FinalFragmentCoverageTest(unittest.TestCase):
    def test_joins_filtered_metrics_and_raw_counts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            tmp_path = Path(temporary)
            annotations = write(tmp_path / "annotations.tsv", "vf_id\tcheckv_quality\nVF_0000001\tHigh-quality\nVF_0000002\tNot-determined\n")
            coverage = write(tmp_path / "coverage.tsv", "Genome\tMean\tRelative Abundance (%)\tCovered Bases\n/path/VF_0000001.fna\t2.5\t17.3\t400\n")
            counts = write(tmp_path / "counts.tsv", "Genome\tCount\n/path/VF_0000001.fna\t12\n/path/VF_0000002.fna\t1\n")
            output = tmp_path / "quantified.tsv"
            run("join_fragment_coverage.py", "--annotations", str(annotations), "--coverage", str(coverage), "--counts", str(counts), "--output", str(output))
            result = {row["vf_id"]: row for row in table(output)}
            self.assertEqual(result["VF_0000001"], {"vf_id": "VF_0000001", "checkv_quality": "High-quality", "relative_abundance": "17.3", "mean_coverage": "2.5", "covered_bases": "400", "read_count": "12", "detected": "yes"})
            self.assertEqual(result["VF_0000002"]["read_count"], "1")
            self.assertEqual(result["VF_0000002"]["mean_coverage"], "0")


class ReportPresentationTest(unittest.TestCase):
    def test_global_report_explains_baltimore_groups_and_links_sample_reports(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write(root / "03_candidate_catalogue/VC_discovery_evidence.tsv", "vc_id\tgeNomad\tVirSorter2\tDIAMOND_NR_virus\tdiscovery_pattern\nVC_1\tyes\tyes\tno\tgeNomad + VirSorter2\n")
            write(root / "04_nr_annotation/viral_decision.tsv", "vc_id\tdecision\nVC_1\tconfirmed_viral\n")
            write(root / "07_final_catalogue/VF_catalogue.tsv", "vf_id\tparent_vc_id\tsource_sample_ids\tsource_contig_ids\tdecision\tcheckv_quality\tnr_family\tictv_species\tbaltimore_group\tictv_pident\nVF_1\tVC_1\tSample 1\tSample 1__contig_7\tconfirmed_viral\tHigh-quality\tReoviridae\tExample virus\tIII\t99\n")
            write(root / "08_sample_results/sample_fragment_presence.tsv", "sample_id\tvf_id\tparent_vc_id\nSample 1\tVF_1\tVC_1\n")
            write(root / "08_sample_results/Sample 1/viral_fragments_quantified.tsv", "vf_id\tcheckv_quality\tnr_family\tictv_species\tbaltimore_group\trelative_abundance\tmean_coverage\tread_count\tdetected\nVF_1\tHigh-quality\tReoviridae\tExample virus\tIII\t12.5\t3\t7\tyes\n")
            run("build_virome_catalogue_report.py", "--output-dir", str(root))
            global_report = (root / "reports/virome_catalogue_dashboard.html").read_text(encoding="utf-8")
            self.assertIn("III：双链 RNA（dsRNA）病毒", global_report)
            self.assertIn("打开单样本报告", global_report)
            self.assertIn("Sample 1__contig_7", global_report)
            self.assertIn("下载当前结果 TSV", global_report)
            self.assertIn("class='venn'", global_report)
            self.assertIn("log10", global_report)
            sample_pages = list((root / "reports/samples").glob("*.html"))
            self.assertEqual(len(sample_pages), 1)
            sample_page = sample_pages[0].read_text(encoding="utf-8")
            self.assertIn("单样本病毒分析报告", sample_page)
            self.assertIn("累计原始读数", sample_page)
