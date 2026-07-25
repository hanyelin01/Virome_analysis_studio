from __future__ import annotations

import json
import os
import stat
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DIAGNOSTIC = ROOT / "scripts" / "diagnose_virome_environment.py"


class ViromeEnvironmentDiagnosticTest(unittest.TestCase):
    def create_ready_config(self, root: Path) -> Path:
        tools = root / "tools"
        tools.mkdir()
        for name in ("genomad", "virsorter", "checkv", "diamond", "taxonkit", "coverm", "daa2rma"):
            executable = tools / name
            executable.write_text("#!/bin/sh\necho 'test-tool 1.0'\n", encoding="utf-8")
            executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
        # The diagnostic checks the Python interpreter paired with VirSorter2,
        # rather than the caller's Python environment.
        virsorter_python = tools / "python"
        virsorter_python.write_text("#!/bin/sh\necho 'fixture-screed 1.0'\n", encoding="utf-8")
        virsorter_python.chmod(virsorter_python.stat().st_mode | stat.S_IXUSR)
        references = root / "references"
        for name in ("genomad", "checkv", "taxonkit"):
            (references / name).mkdir(parents=True)
        (references / "taxonkit" / "nodes.dmp").write_text("1\t|\t1\t|\tno rank\t|\n", encoding="utf-8")
        (references / "taxonkit" / "names.dmp").write_text("1\t|\troot\t|\t\t|\tscientific name\t|\n", encoding="utf-8")
        for name in ("nr.dmnd", "ictv.dmnd", "megan.map"):
            (references / name).write_text("fixture\n", encoding="utf-8")
        metadata = references / "ictv.tsv"
        metadata.write_text("reference_id\tfamily\tgenus\tspecies\tbaltimore_group\nREF1\tCoronaviridae\tBetacoronavirus\tExample virus\tIV\n", encoding="utf-8")
        data_root = root / "projects"
        data_root.mkdir()
        config = root / "pipeline.env"
        config.write_text(
            "\n".join(
                [
                    f"ALLOWED_DATA_ROOTS={data_root}",
                    "MAX_TOTAL_THREADS=16",
                    "MAX_THREADS_PER_VIRAL_TOOL=8",
                    "VIRAL_MIN_CONTIG_LEN=200",
                    f"GENOMAD_DB={references / 'genomad'}",
                    f"CHECKV_DB={references / 'checkv'}",
                    f"DIAMOND_NR_DB={references / 'nr.dmnd'}",
                    "DIAMOND_DEFAULT_TAXONLIST=10239",
                    "DIAMOND_EVALUE=1e-5",
                    f"TAXONKIT_DB={references / 'taxonkit'}",
                    f"ICTV_REFERENCE_DMND={references / 'ictv.dmnd'}",
                    f"ICTV_REFERENCE_METADATA={metadata}",
                    "ICTV_REFERENCE_VERSION=fixture.v1",
                    f"MEGAN_DAA2RMA={tools / 'daa2rma'}",
                    f"MEGAN_MAP_DB={references / 'megan.map'}",
                    f"VIRSORTER_COMMAND={tools / 'virsorter'}",
                    "VIRSORTER_USE_CONDA_OFF=0",
                    "COVERM_MIN_READ_PERCENT_IDENTITY=95",
                    "COVERM_MIN_READ_ALIGNED_PERCENT=75",
                    "COVERM_MIN_COVERED_FRACTION=10",
                ]
            ) + "\n",
            encoding="utf-8",
        )
        self.tools = tools
        return config

    def run_diagnostic(self, config: Path) -> tuple[subprocess.CompletedProcess[str], dict[str, object]]:
        environment = {**os.environ, "PATH": f"{self.tools}:{os.environ.get('PATH', '')}"}
        completed = subprocess.run(
            [sys.executable, str(DIAGNOSTIC), "--config", str(config), "--format", "json"],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=environment,
            check=False,
        )
        return completed, json.loads(completed.stdout)

    def test_ready_environment_has_machine_readable_passes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = self.create_ready_config(Path(temporary))
            completed, report = self.run_diagnostic(config)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(report["status"], "ready")
        checks = {item["check_id"]: item for item in report["checks"]}
        self.assertEqual(checks["tool.virsorter2"]["status"], "pass")
        self.assertEqual(checks["tool.virsorter2_screed"]["status"], "pass")
        self.assertEqual(checks["reference.ictv_metadata_schema"]["status"], "pass")
        self.assertEqual(checks["parameter.viral_min_contig_len"]["value"], "200")

    def test_invalid_ictv_metadata_blocks_start(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = self.create_ready_config(Path(temporary))
            metadata = next(line.split("=", 1)[1] for line in config.read_text(encoding="utf-8").splitlines() if line.startswith("ICTV_REFERENCE_METADATA="))
            Path(metadata).write_text("reference_id\tfamily\nREF1\tCoronaviridae\n", encoding="utf-8")
            completed, report = self.run_diagnostic(config)
        self.assertEqual(completed.returncode, 3)
        self.assertEqual(report["status"], "blocked")
        checks = {item["check_id"]: item for item in report["checks"]}
        self.assertEqual(checks["reference.ictv_metadata_schema"]["status"], "fail")
