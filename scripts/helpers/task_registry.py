"""Persistent, local task registry for the Streamlit control panel.

The registry records only task metadata and paths already selected in the UI.
Pipeline scripts remain the source of truth for run status and logs.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from uuid import uuid4


RUN_STATUSES = {"RUNNING", "SUCCESS", "FAILED"}


def default_registry_path() -> Path:
    state_home = Path.home() / ".local" / "state"
    return Path(os.environ.get("XDG_STATE_HOME", state_home)) / "virome-contig-studio" / "tasks.sqlite3"


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def connect(database: Path | None = None) -> sqlite3.Connection:
    path = database or default_registry_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(path)
    connection.row_factory = sqlite3.Row
    connection.execute(
        """
        CREATE TABLE IF NOT EXISTS tasks (
            task_id TEXT PRIMARY KEY,
            workflow TEXT NOT NULL,
            workflow_label TEXT NOT NULL,
            state_base TEXT NOT NULL,
            command_json TEXT NOT NULL,
            submitted_at TEXT NOT NULL,
            pid INTEGER,
            run_dir TEXT,
            status TEXT NOT NULL,
            current_step TEXT,
            completed_at TEXT,
            last_seen_at TEXT NOT NULL
        )
        """
    )
    connection.execute("CREATE INDEX IF NOT EXISTS tasks_submitted_at ON tasks(submitted_at DESC)")
    return connection


def parse_parameters(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.is_file():
        return values
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" in line and not line.startswith("#"):
            key, value = line.split("=", 1)
            values[key] = value
    return values


def status_for_run(run_dir: Path) -> tuple[str, str | None]:
    status_file = run_dir / "status"
    raw = status_file.read_text(encoding="utf-8", errors="replace").strip().upper() if status_file.is_file() else "STARTING"
    status = raw if raw in RUN_STATUSES else "STARTING"
    log = run_dir / "pipeline.log"
    step: str | None = None
    if log.is_file():
        matches = re.findall(r"^\[STEP\]\s+(.+?)\s*$", log.read_text(encoding="utf-8", errors="replace"), re.MULTILINE)
        if matches:
            step = matches[-1]
    return status, step


def register_submission(
    workflow: str,
    workflow_label: str,
    state_base: Path,
    command: list[str],
    pid: int,
    database: Path | None = None,
) -> str:
    task_id = uuid4().hex
    timestamp = now()
    with connect(database) as connection:
        connection.execute(
            """INSERT INTO tasks(task_id, workflow, workflow_label, state_base, command_json,
                                 submitted_at, pid, status, last_seen_at)
               VALUES (?, ?, ?, ?, ?, ?, ?, 'STARTING', ?)""",
            (task_id, workflow, workflow_label, str(state_base), json.dumps(command), timestamp, pid, timestamp),
        )
    return task_id


def _find_run(state_base: Path, task_id: str) -> Path | None:
    runs = state_base / ".contig_pipeline" / "runs"
    if not runs.is_dir():
        return None
    for run_dir in sorted((item for item in runs.iterdir() if item.is_dir()), key=lambda item: item.stat().st_mtime, reverse=True):
        if parse_parameters(run_dir / "parameters.env").get("TASK_REGISTRY_ID") == task_id:
            return run_dir
    return None


def _pid_is_running(pid: int | None) -> bool:
    if not pid or pid < 1:
        return False
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    return True


def refresh_tasks(database: Path | None = None) -> list[dict[str, Any]]:
    with connect(database) as connection:
        records = [dict(row) for row in connection.execute("SELECT * FROM tasks ORDER BY submitted_at DESC")]
        for record in records:
            run_dir = Path(record["run_dir"]) if record["run_dir"] else _find_run(Path(record["state_base"]), record["task_id"])
            if run_dir and run_dir.is_dir():
                status, step = status_for_run(run_dir)
                completed = now() if status in {"SUCCESS", "FAILED"} and not record["completed_at"] else record["completed_at"]
                connection.execute(
                    "UPDATE tasks SET run_dir=?, status=?, current_step=?, completed_at=?, last_seen_at=? WHERE task_id=?",
                    (str(run_dir), status, step, completed, now(), record["task_id"]),
                )
                record.update(run_dir=str(run_dir), status=status, current_step=step, completed_at=completed)
            else:
                # A process that exits before it can create a run directory is
                # a failed submission (for example, a lock conflict or a missing
                # executable), not an indefinitely "starting" task.
                if record["status"] == "STARTING" and record["pid"] and not _pid_is_running(record["pid"]):
                    completed = now()
                    connection.execute(
                        "UPDATE tasks SET status='FAILED', completed_at=?, last_seen_at=? WHERE task_id=?",
                        (completed, completed, record["task_id"]),
                    )
                    record.update(status="FAILED", completed_at=completed)
                else:
                    record["status"] = record["status"] or "STARTING"
        connection.commit()
        return [dict(row) for row in connection.execute("SELECT * FROM tasks ORDER BY submitted_at DESC")]


def _workflow_from_parameters(values: dict[str, str]) -> tuple[str, str]:
    task = values.get("TASK", "legacy")
    mapping = {
        "qc_only": ("qc_only", "① 原始数据质控"),
        "assembly_only": ("assembly_only", "② MEGAHIT 拼接"),
        "full": ("full", "③ 完整拼接流程"),
        "viral_report": ("viral_report", "旧版病毒报告"),
        "virome_catalogue_v2": ("virome_catalogue", "④ 全局病毒发现、分类与精细注释"),
        "fine_annotation": ("fine_annotation", "⑤ 独立 DIAMOND 精细注释"),
    }
    return mapping.get(task, (task, f"历史任务：{task}"))


def import_history(state_base: Path, database: Path | None = None) -> int:
    runs = state_base / ".contig_pipeline" / "runs"
    if not runs.is_dir():
        return 0
    imported = 0
    with connect(database) as connection:
        for run_dir in (item for item in runs.iterdir() if item.is_dir()):
            values = parse_parameters(run_dir / "parameters.env")
            task_id = values.get("TASK_REGISTRY_ID") or "legacy-" + hashlib.sha256(str(run_dir.resolve()).encode()).hexdigest()[:24]
            workflow, label = _workflow_from_parameters(values)
            status, step = status_for_run(run_dir)
            submitted = datetime.fromtimestamp(run_dir.stat().st_mtime, timezone.utc).isoformat(timespec="seconds")
            result = connection.execute(
                """INSERT OR IGNORE INTO tasks(task_id, workflow, workflow_label, state_base, command_json,
                                                   submitted_at, run_dir, status, current_step, completed_at, last_seen_at)
                   VALUES (?, ?, ?, ?, '[]', ?, ?, ?, ?, ?, ?)""",
                (task_id, workflow, label, str(state_base), submitted, str(run_dir), status, step,
                 submitted if status in {"SUCCESS", "FAILED"} else None, now()),
            )
            imported += result.rowcount
    return imported


def discover_history(
    roots: list[Path],
    database: Path | None = None,
    *,
    max_depth: int = 8,
    max_directories: int = 25_000,
) -> tuple[int, int, bool]:
    """Find known run directories below approved data roots.

    Only directory names are traversed; raw FASTQ/FASTA files are never read.
    Limits keep an unusually large data root from delaying the web interface.
    Returns (new_tasks, output_locations, was_truncated).
    """
    imported = locations = visited = 0
    truncated = False
    for root in roots:
        if not root.is_dir():
            continue
        for current_text, directories, _ in os.walk(root, topdown=True, followlinks=False):
            current = Path(current_text)
            visited += 1
            if visited > max_directories:
                return imported, locations, True
            try:
                depth = len(current.relative_to(root).parts)
            except ValueError:
                directories[:] = []
                continue
            if current.name == ".contig_pipeline":
                if (current / "runs").is_dir():
                    imported += import_history(current.parent, database)
                    locations += 1
                directories[:] = []
            elif depth >= max_depth:
                directories[:] = []
    return imported, locations, truncated
