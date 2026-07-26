import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PREPARE = ROOT / "scripts" / "helpers" / "prepare_viral_contigs.py"
PREPARE_WRAPPER = ROOT / "scripts" / "04_prepare_viral_contigs.sh"


class PrepareViralContigsTest(unittest.TestCase):
    def make_project(self, root: Path) -> tuple[Path, Path]:
        assembly = root / "assembly"
        for sample, base in (("sample_A", "A"), ("sample_B", "C"), ("old_sample", "G")):
            sample_dir = assembly / sample
            sample_dir.mkdir(parents=True)
            (sample_dir / "final.contigs.fa").write_text(
                f">contig_1\n{base * 20}\n", encoding="utf-8"
            )
        manifest = root / "manifest.tsv"
        manifest.write_text(
            "sample_id\tread_type\traw_r1\traw_r2\tclean_r1\tclean_r2\tclean_single\tassembly_dir\n"
            f"sample_A\tpe\t\t\t/a1\t/a2\t\t{assembly / 'sample_A'}\n"
            f"sample_B\tpe\t\t\t/b1\t/b2\t\t{assembly / 'sample_B'}\n",
            encoding="utf-8",
        )
        return assembly, manifest

    def run_prepare(self, assembly: Path, manifest: Path, output: Path, *extra: str):
        return subprocess.run([
            sys.executable, str(PREPARE), "--assembly-dir", str(assembly),
            "--manifest", str(manifest), "--output-dir", str(output),
            "--min-length", "10", *extra,
        ], capture_output=True, text=True)

    def test_manifest_excludes_extra_assembly_sample_and_resume_is_verified(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            assembly, manifest = self.make_project(root)
            output = root / "prepared"
            first = self.run_prepare(assembly, manifest, output)
            self.assertEqual(first.returncode, 0, first.stderr)
            fasta = (output / "merged_assembled_contigs.fna").read_text(encoding="utf-8")
            self.assertIn(">sample_A__contig_1", fasta)
            self.assertIn(">sample_B__contig_1", fasta)
            self.assertNotIn("old_sample", fasta)
            with (output / "contig_preparation_summary.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                self.assertEqual(
                    [row["sample_id"] for row in csv.DictReader(handle, delimiter="\t")],
                    ["sample_A", "sample_B"],
                )
            fingerprint = json.loads(
                (output / "preparation_inputs.json").read_text(encoding="utf-8")
            )
            self.assertEqual([row["sample_id"] for row in fingerprint["inputs"]], ["sample_A", "sample_B"])
            resumed = self.run_prepare(assembly, manifest, output, "--resume")
            self.assertEqual(resumed.returncode, 0, resumed.stderr)
            self.assertIn("match the current manifest", resumed.stdout)

    def test_resume_rejects_a_different_manifest(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            assembly, manifest = self.make_project(root)
            output = root / "prepared"
            self.assertEqual(self.run_prepare(assembly, manifest, output).returncode, 0)
            manifest.write_text(
                "sample_id\tread_type\traw_r1\traw_r2\tclean_r1\tclean_r2\tclean_single\tassembly_dir\n"
                f"sample_A\tpe\t\t\t/a1\t/a2\t\t{assembly / 'sample_A'}\n",
                encoding="utf-8",
            )
            resumed = self.run_prepare(assembly, manifest, output, "--resume")
            self.assertNotEqual(resumed.returncode, 0)
            self.assertIn("different samples", resumed.stderr)

    def test_wrapper_archives_pre_fingerprint_output_then_rebuilds(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            assembly, manifest = self.make_project(root)
            report = root / "report"
            output = report / "01_prepared_contigs"
            self.assertEqual(self.run_prepare(assembly, manifest, output).returncode, 0)
            (output / "preparation_inputs.json").unlink()
            config = root / "pipeline.env"
            config.write_text(f"ALLOWED_DATA_ROOTS={root}\nVIRAL_MIN_CONTIG_LEN=10\n", encoding="utf-8")
            completed = subprocess.run(
                ["bash", str(PREPARE_WRAPPER), "--assembly-dir", str(assembly), "--manifest", str(manifest), "--output-dir", str(report), "--min-contig-length", "10", "--resume"],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env={**os.environ, "CONTIG_PIPELINE_CONFIG": str(config)},
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertIn("Archived legacy prepared contigs", completed.stdout)
            self.assertTrue((output / "preparation_inputs.json").is_file())
            archives = list((report / ".contig_pipeline" / "legacy_prepared_contigs").glob("*_01_prepared_contigs"))
            self.assertEqual(len(archives), 1)
            self.assertTrue((archives[0] / "merged_assembled_contigs.fna").is_file())

    def test_web_entrypoint_scripts_are_executable(self):
        for name in ("run_pipeline.sh", "run_viral_report.sh", "run_fine_annotation.sh"):
            path = ROOT / "scripts" / name
            self.assertTrue(os.access(path, os.X_OK), f"{path} must remain executable")


if __name__ == "__main__":
    unittest.main()
