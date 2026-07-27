import tempfile
import unittest
from pathlib import Path

from scripts.helpers import task_registry


class TaskRegistryTest(unittest.TestCase):
    def test_submission_is_bound_to_its_marked_run(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            tmp_path = Path(temporary)
            database = tmp_path / "registry.sqlite3"
            output = tmp_path / "output"
            task_id = task_registry.register_submission("qc_only", "质控", output, ["run"], 123, database)
            run = output / ".contig_pipeline" / "runs" / "20260727_120000_1"
            run.mkdir(parents=True)
            (run / "parameters.env").write_text(f"TASK_REGISTRY_ID={task_id}\nTASK=qc_only\n", encoding="utf-8")
            (run / "status").write_text("RUNNING\n", encoding="utf-8")
            (run / "pipeline.log").write_text("[STEP] fastp\n", encoding="utf-8")

            records = task_registry.refresh_tasks(database)

            self.assertEqual(records[0]["run_dir"], str(run))
            self.assertEqual(records[0]["status"], "RUNNING")
            self.assertEqual(records[0]["current_step"], "fastp")

    def test_import_history_is_idempotent(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            tmp_path = Path(temporary)
            database = tmp_path / "registry.sqlite3"
            run = tmp_path / "output" / ".contig_pipeline" / "runs" / "old"
            run.mkdir(parents=True)
            (run / "parameters.env").write_text("TASK=virome_catalogue_v2\n", encoding="utf-8")
            (run / "status").write_text("SUCCESS\n", encoding="utf-8")

            self.assertEqual(task_registry.import_history(tmp_path / "output", database), 1)
            self.assertEqual(task_registry.import_history(tmp_path / "output", database), 0)
            self.assertEqual(task_registry.refresh_tasks(database)[0]["workflow"], "virome_catalogue")

    def test_discovery_finds_history_without_an_output_path(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            tmp_path = Path(temporary)
            database = tmp_path / "registry.sqlite3"
            run = tmp_path / "projects" / "batch" / "report" / ".contig_pipeline" / "runs" / "old"
            run.mkdir(parents=True)
            (run / "parameters.env").write_text("TASK=full\n", encoding="utf-8")
            (run / "status").write_text("FAILED\n", encoding="utf-8")

            imported, locations, truncated = task_registry.discover_history([tmp_path / "projects"], database)

            self.assertEqual((imported, locations, truncated), (1, 1, False))
            self.assertEqual(task_registry.refresh_tasks(database)[0]["status"], "FAILED")
