import tempfile
import unittest
import subprocess
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

    def test_custom_display_name_can_be_saved(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            tmp_path = Path(temporary)
            database = tmp_path / "registry.sqlite3"
            task_id = task_registry.register_submission("qc_only", "质控", tmp_path / "output", ["run"], 123, database, "批次 A")
            self.assertEqual(task_registry.refresh_tasks(database)[0]["display_name"], "批次 A")
            task_registry.rename_task(task_id, "批次 A 重跑", database)
            self.assertEqual(task_registry.refresh_tasks(database)[0]["display_name"], "批次 A 重跑")

    def test_suggested_name_uses_batch_and_output_directories(self) -> None:
        output = Path("/home/hanyl/Projects/0WulabNGSData/2026BatCN_NHZY_Yunnan/04.Viral_report")
        self.assertEqual(
            task_registry.suggested_display_name("v2", output, "2026-07-27T00:00:00+00:00"),
            "2026BatCN_NHZY_Yunnan · 04.Viral_report · 首次运行",
        )

    def test_termination_only_targets_registered_process_group(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            tmp_path = Path(temporary)
            database = tmp_path / "registry.sqlite3"
            process = subprocess.Popen(["sleep", "30"], start_new_session=True)
            try:
                task_id = task_registry.register_submission("qc_only", "质控", tmp_path / "output", ["sleep", "30"], process.pid, database)
                success, _ = task_registry.terminate_task(task_id, database)
                self.assertTrue(success)
                process.wait(timeout=5)
                self.assertEqual(task_registry.refresh_tasks(database)[0]["status"], "CANCELLED")
            finally:
                if process.poll() is None:
                    process.kill()

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
