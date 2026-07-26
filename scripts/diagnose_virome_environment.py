#!/usr/bin/env python3
"""Read-only readiness checks for the virome catalogue v2 workflow.

The diagnostic intentionally verifies availability and provenance only.  It
does not run a classifier, mutate a reference database, or inspect research
data.  Its JSON output is also consumed by the Streamlit pre-run panel.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "config" / "pipeline.env"
SCHEMA_VERSION = "1.0"


@dataclass
class Check:
    check_id: str
    category: str
    status: str
    severity: str
    message: str
    remediation: str = ""
    value: str = ""


class Diagnostics:
    def __init__(self) -> None:
        self.checks: list[Check] = []

    def add(
        self,
        check_id: str,
        category: str,
        status: str,
        severity: str,
        message: str,
        remediation: str = "",
        value: str = "",
    ) -> None:
        self.checks.append(Check(check_id, category, status, severity, message, remediation, value))

    def path(
        self,
        check_id: str,
        category: str,
        raw_value: str,
        kind: str,
        *,
        required: bool = True,
        extension: str | None = None,
    ) -> Path | None:
        label = "目录" if kind == "dir" else "文件"
        if not raw_value:
            status, severity = ("fail", "required") if required else ("warn", "optional")
            self.add(
                check_id,
                category,
                status,
                severity,
                f"未配置{label}路径。",
                "在 config/pipeline.env 中填写绝对路径。",
            )
            return None
        path = Path(raw_value).expanduser()
        valid_kind = path.is_dir() if kind == "dir" else path.is_file()
        if not valid_kind or not os.access(path, os.R_OK):
            self.add(
                check_id,
                category,
                "fail",
                "required",
                f"{label}不可访问：{path}",
                "确认路径、挂载状态和运行账号的读取权限。",
                str(path),
            )
            return None
        if extension and path.suffix.lower() != extension:
            self.add(
                check_id,
                category,
                "warn",
                "recommended",
                f"文件可读取，但文件扩展名不是 {extension}：{path.name}",
                "确认此文件确为预期的参考数据库。",
                str(path),
            )
        else:
            self.add(check_id, category, "pass", "required", f"{label}可读取。", value=str(path))
        return path


def parse_env(path: Path) -> tuple[dict[str, str], str | None]:
    if not path.is_file():
        return {}, f"配置文件不存在：{path}"
    settings: dict[str, str] = {"PIPELINE_HOME": str(ROOT)}
    assignment = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*)=(.*)$")
    for line_number, raw in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        match = assignment.match(line)
        if not match:
            return {}, f"第 {line_number} 行不是 KEY=VALUE 格式。"
        key, value = match.groups()
        value = value.strip().strip('"').strip("'")
        value = os.path.expandvars(value.replace("$PIPELINE_HOME", str(ROOT)))
        settings[key] = value
    return settings, None


def executable_path(command: str) -> str | None:
    if not command:
        return None
    candidate = Path(command).expanduser()
    if "/" in command:
        return str(candidate) if candidate.is_file() and os.access(candidate, os.X_OK) else None
    return shutil.which(command)


def version_of(executable: str) -> str:
    try:
        result = subprocess.run(
            [executable, "--version"],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=4,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return "未报告"
    output = (result.stdout or result.stderr).strip().replace("\n", " ")
    return output[:180] if output else "未报告"


def check_command(diagnostics: Diagnostics, check_id: str, label: str, configured: str, remediation: str) -> None:
    executable = executable_path(configured)
    if executable is None:
        diagnostics.add(
            check_id,
            "tool",
            "fail",
            "required",
            f"未找到可执行程序：{label}（配置值：{configured or '未配置'}）。",
            remediation,
            configured,
        )
        return
    diagnostics.add(
        check_id,
        "tool",
        "pass",
        "required",
        f"{label} 可执行。",
        value=f"{executable}；版本：{version_of(executable)}",
    )


def check_diamond_version(diagnostics: Diagnostics, configured: str, minimum: str) -> None:
    executable = executable_path(configured)
    if executable is None:
        return
    observed = version_of(executable)
    found = re.search(r"\b(\d+(?:\.\d+)+)\b", observed)
    required = re.fullmatch(r"\d+(?:\.\d+)+", minimum or "")
    if not required:
        diagnostics.add("tool.diamond_version", "tool", "fail", "required", f"DIAMOND_MIN_VERSION 无效：{minimum or '未配置'}。", "使用形如 2.2.4 的版本号。", minimum)
        return
    if found is None:
        diagnostics.add("tool.diamond_version", "tool", "fail", "required", "无法读取 DIAMOND 版本。", "确认 diamond --version 可执行，且版本不低于要求。", observed)
        return
    actual_tuple = tuple(map(int, found.group(1).split(".")))
    required_tuple = tuple(map(int, minimum.split(".")))
    if actual_tuple >= required_tuple:
        diagnostics.add("tool.diamond_version", "tool", "pass", "required", f"DIAMOND 版本满足最低要求（>= {minimum}）。", value=found.group(1))
    else:
        diagnostics.add("tool.diamond_version", "tool", "fail", "required", f"DIAMOND {found.group(1)} 低于最低要求 {minimum}。", "在 contig-ui 环境执行 conda install -n contig-ui -c conda-forge -c bioconda diamond=2.2.4。", found.group(1))


def check_virsorter_python_dependency(diagnostics: Diagnostics, configured: str) -> None:
    """Verify a dependency from the interpreter that runs VirSorter2 itself."""
    executable = executable_path(configured)
    if executable is None:
        return
    interpreter = Path(executable).parent / "python"
    if not interpreter.is_file() or not os.access(interpreter, os.X_OK):
        diagnostics.add(
            "tool.virsorter2_screed",
            "tool",
            "fail",
            "required",
            f"无法定位与 VirSorter2 配套的 Python 解释器：{interpreter}。",
            "将 VIRSORTER_COMMAND 配置为独立环境中的 virsorter 绝对路径，或修复该环境。",
            str(interpreter),
        )
        return
    try:
        result = subprocess.run(
            [str(interpreter), "-c", "import screed; print(screed.__version__)"],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=4,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        diagnostics.add("tool.virsorter2_screed", "tool", "fail", "required", f"VirSorter2 依赖检查无法执行：{error}。", "修复 VirSorter2 独立环境后重新诊断。", str(interpreter))
        return
    if result.returncode == 0:
        diagnostics.add("tool.virsorter2_screed", "tool", "pass", "required", "VirSorter2 所需的 screed 包可导入。", value=result.stdout.strip() or "已安装")
    else:
        diagnostics.add("tool.virsorter2_screed", "tool", "fail", "required", "VirSorter2 环境缺少所需 Python 包 screed。", "执行 conda install -n virsorter2 -c conda-forge screed，然后重新诊断。", str(interpreter))


def check_ictv_metadata(diagnostics: Diagnostics, path: Path | None) -> None:
    if path is None:
        return
    required = {"reference_id", "family", "genus", "species", "baltimore_group"}
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            headers = set(reader.fieldnames or [])
            first = next(reader, None)
    except (OSError, csv.Error, UnicodeError) as error:
        diagnostics.add("reference.ictv_metadata_schema", "reference", "fail", "required", f"ICTV 元数据无法读取：{error}", "使用 UTF-8 TSV 并保留必需列。", str(path))
        return
    missing = sorted(required - headers)
    if missing:
        diagnostics.add("reference.ictv_metadata_schema", "reference", "fail", "required", f"ICTV 元数据缺少列：{', '.join(missing)}。", "使用 build_ictv_reference_db.sh 所要求的元数据格式。", str(path))
    elif first is None:
        diagnostics.add("reference.ictv_metadata_schema", "reference", "fail", "required", "ICTV 元数据只有表头，没有参考序列记录。", "重新审核并导出包含至少一条记录的 TSV。", str(path))
    elif any(not (first.get(column) or "").strip() for column in required):
        diagnostics.add("reference.ictv_metadata_schema", "reference", "fail", "required", "ICTV 元数据首条记录存在必填字段为空。", "补全 reference_id、分类学字段和 Baltimore 分类。", str(path))
    else:
        diagnostics.add("reference.ictv_metadata_schema", "reference", "pass", "required", "ICTV 元数据表头和首条记录符合 v2 契约。", value=str(path))


def positive_integer(value: str) -> bool:
    return bool(re.fullmatch(r"[1-9][0-9]*", value or ""))


def decimal_in_range(value: str, minimum: float, maximum: float) -> bool:
    try:
        parsed = float(value)
    except ValueError:
        return False
    return minimum <= parsed <= maximum


def check_parameter(
    diagnostics: Diagnostics,
    check_id: str,
    label: str,
    value: str,
    predicate: Callable[[str], bool],
    expected: str,
) -> None:
    if predicate(value):
        diagnostics.add(check_id, "parameter", "pass", "required", f"{label} 处于受支持范围。", value=value)
    else:
        diagnostics.add(check_id, "parameter", "fail", "required", f"{label} 无效：{value or '未配置'}。", f"应为：{expected}。", value)


def check_diamond_tmpdir(diagnostics: Diagnostics, value: str, block_size: str) -> None:
    path = Path(value).expanduser()
    if not path.is_dir() or not os.access(path, os.R_OK | os.W_OK):
        diagnostics.add("parameter.diamond_tmpdir", "parameter", "fail", "required", f"DIAMOND_TMPDIR 不存在或不可写：{path}。", "选择可写的高速临时目录；推荐容量充足时使用 /dev/shm。", str(path))
        return
    try:
        available = os.statvfs(path).f_bavail * os.statvfs(path).f_frsize
        available_gib = available / (1024 ** 3)
        required_gib = float(block_size) * 1.2
    except (OSError, ValueError):
        diagnostics.add("parameter.diamond_tmpdir", "parameter", "warn", "recommended", f"DIAMOND 临时目录可写，但无法评估可用空间：{path}。", "确认临时目录可用空间高于 block size，并预留额外空间。", str(path))
        return
    if available_gib < required_gib:
        diagnostics.add("parameter.diamond_tmpdir", "parameter", "warn", "recommended", f"DIAMOND 临时目录可写，但可用空间仅 {available_gib:.1f} GiB。", f"建议至少预留约 {required_gib:.1f} GiB；可降低 DIAMOND_BLOCK_SIZE 或改用更大临时盘。", str(path))
    else:
        diagnostics.add("parameter.diamond_tmpdir", "parameter", "pass", "required", f"DIAMOND 临时目录可写，可用空间 {available_gib:.1f} GiB。", value=str(path))


def run(config_path: Path) -> dict[str, object]:
    diagnostics = Diagnostics()
    settings, config_error = parse_env(config_path)
    if config_error:
        diagnostics.add("config.pipeline_env", "configuration", "fail", "required", config_error, "从 config/pipeline.env.example 复制后填写本机路径。", str(config_path))
        return render(config_path, diagnostics)
    diagnostics.add("config.pipeline_env", "configuration", "pass", "required", "服务器私有配置文件可读取。", value=str(config_path))

    roots = [Path(item).expanduser() for item in settings.get("ALLOWED_DATA_ROOTS", "").split(":") if item]
    if not roots:
        diagnostics.add("config.allowed_data_roots", "configuration", "fail", "required", "ALLOWED_DATA_ROOTS 未配置。", "限定可读写的数据根目录，例如 /home/hanyl/Projects。")
    elif all(root.is_dir() and os.access(root, os.R_OK | os.W_OK) for root in roots):
        diagnostics.add("config.allowed_data_roots", "configuration", "pass", "required", "全部数据根目录可读写。", value=":".join(map(str, roots)))
    else:
        diagnostics.add("config.allowed_data_roots", "configuration", "fail", "required", "至少一个数据根目录不存在或不可读写。", "确认目录、挂载状态和运行账号权限。", ":".join(map(str, roots)))

    for check_id, label, command, remediation in [
        ("tool.genomad", "geNomad", "genomad", "激活 contig-ui 环境，并确认 genomad 在 PATH 中。"),
        ("tool.virsorter2", "VirSorter2", settings.get("VIRSORTER_COMMAND", "virsorter"), "激活 contig-ui 环境；若使用独立环境，请配置 VIRSORTER_COMMAND 与 VIRSORTER_USE_CONDA_OFF。"),
        ("tool.checkv", "CheckV", "checkv", "激活 contig-ui 环境，并确认 checkv 在 PATH 中。"),
        ("tool.diamond", "DIAMOND", "diamond", "激活 contig-ui 环境，并确认 diamond 在 PATH 中。"),
        ("tool.taxonkit", "TaxonKit", "taxonkit", "激活 contig-ui 环境，并确认 taxonkit 在 PATH 中。"),
        ("tool.coverm", "CoverM", "coverm", "激活 contig-ui 环境，并确认 coverm 在 PATH 中。"),
        ("tool.megan_daa2rma", "MEGAN daa2rma", settings.get("MEGAN_DAA2RMA", ""), "在 pipeline.env 中填写可执行的 MEGAN_DAA2RMA 绝对路径，并确认许可证可用。"),
    ]:
        check_command(diagnostics, check_id, label, command, remediation)
    check_diamond_version(diagnostics, "diamond", settings.get("DIAMOND_MIN_VERSION", "2.2.4"))
    check_virsorter_python_dependency(diagnostics, settings.get("VIRSORTER_COMMAND", "virsorter"))

    genomad = diagnostics.path("reference.genomad", "reference", settings.get("GENOMAD_DB", ""), "dir")
    checkv = diagnostics.path("reference.checkv", "reference", settings.get("CHECKV_DB", ""), "dir")
    nr = diagnostics.path("reference.nr_diamond", "reference", settings.get("DIAMOND_NR_DB", ""), "file", extension=".dmnd")
    taxonkit = diagnostics.path("reference.taxonkit", "reference", settings.get("TAXONKIT_DB", ""), "dir")
    ictv_dmnd = diagnostics.path("reference.ictv_diamond", "reference", settings.get("ICTV_REFERENCE_DMND", ""), "file", extension=".dmnd")
    ictv_meta = diagnostics.path("reference.ictv_metadata", "reference", settings.get("ICTV_REFERENCE_METADATA", ""), "file", extension=".tsv")
    megan_map = diagnostics.path("reference.megan_map", "reference", settings.get("MEGAN_MAP_DB", ""), "file")
    _ = genomad, checkv, nr, ictv_dmnd, megan_map
    if taxonkit is not None:
        missing = [name for name in ("nodes.dmp", "names.dmp") if not (taxonkit / name).is_file()]
        if missing:
            diagnostics.add("reference.taxonkit_taxdump", "reference", "fail", "required", f"TaxonKit 数据库缺少：{', '.join(missing)}。", "重新下载完整 NCBI taxonomy dump。", str(taxonkit))
        else:
            diagnostics.add("reference.taxonkit_taxdump", "reference", "pass", "required", "TaxonKit taxonomy dump 包含 nodes.dmp 与 names.dmp。", value=str(taxonkit))
    check_ictv_metadata(diagnostics, ictv_meta)
    version = settings.get("ICTV_REFERENCE_VERSION", "")
    if version and version != "unconfigured":
        diagnostics.add("reference.ictv_version", "reference", "pass", "required", "ICTV 参考版本已记录。", value=version)
    else:
        diagnostics.add("reference.ictv_version", "reference", "fail", "required", "ICTV_REFERENCE_VERSION 未设置。", "填写 VMR/MSL 版本和本地构建版本，例如 VMR_MSL41.v1.20260721。")

    check_parameter(diagnostics, "parameter.viral_min_contig_len", "VIRAL_MIN_CONTIG_LEN", settings.get("VIRAL_MIN_CONTIG_LEN", ""), lambda value: positive_integer(value) and 200 <= int(value) <= 100000, "200–100000 的整数")
    check_parameter(diagnostics, "parameter.max_viral_threads", "MAX_THREADS_PER_VIRAL_TOOL", settings.get("MAX_THREADS_PER_VIRAL_TOOL", ""), positive_integer, "正整数")
    check_parameter(diagnostics, "parameter.total_threads", "MAX_TOTAL_THREADS", settings.get("MAX_TOTAL_THREADS", ""), positive_integer, "正整数")
    check_parameter(diagnostics, "parameter.diamond_taxonlist", "DIAMOND_DEFAULT_TAXONLIST", settings.get("DIAMOND_DEFAULT_TAXONLIST", ""), lambda value: bool(re.fullmatch(r"[1-9][0-9]*(,[1-9][0-9]*)*", value or "")), "以英文逗号分隔的正整数 TaxID")
    check_parameter(diagnostics, "parameter.diamond_evalue", "DIAMOND_EVALUE", settings.get("DIAMOND_EVALUE", ""), lambda value: decimal_in_range(value, 0.0, 1.0) and float(value) > 0.0, "0 到 1 之间的正数")
    diamond_threads = settings.get("DIAMOND_THREADS_PER_JOB", "64")
    diamond_block_size = settings.get("DIAMOND_BLOCK_SIZE", "4.0")
    diamond_index_chunks = settings.get("DIAMOND_INDEX_CHUNKS", "1")
    diamond_tmpdir = settings.get("DIAMOND_TMPDIR", "/dev/shm")
    check_parameter(diagnostics, "parameter.diamond_threads_per_job", "DIAMOND_THREADS_PER_JOB", diamond_threads, positive_integer, "正整数")
    check_parameter(diagnostics, "parameter.diamond_block_size", "DIAMOND_BLOCK_SIZE", diamond_block_size, lambda value: decimal_in_range(value, 0.01, 128.0), "0.01–128 的 GB 数值")
    check_parameter(diagnostics, "parameter.diamond_index_chunks", "DIAMOND_INDEX_CHUNKS", diamond_index_chunks, positive_integer, "正整数")
    check_diamond_tmpdir(diagnostics, diamond_tmpdir, diamond_block_size)
    for key in ("COVERM_MIN_READ_PERCENT_IDENTITY", "COVERM_MIN_READ_ALIGNED_PERCENT", "COVERM_MIN_COVERED_FRACTION"):
        check_parameter(diagnostics, f"parameter.{key.lower()}", key, settings.get(key, ""), lambda value: decimal_in_range(value, 0.0, 100.0), "0–100 的数值")
    check_parameter(diagnostics, "parameter.virsorter_conda_off", "VIRSORTER_USE_CONDA_OFF", settings.get("VIRSORTER_USE_CONDA_OFF", ""), lambda value: value in {"0", "1"}, "0 或 1")

    return render(config_path, diagnostics)


def render(config_path: Path, diagnostics: Diagnostics) -> dict[str, object]:
    counts = {status: sum(check.status == status for check in diagnostics.checks) for status in ("pass", "warn", "fail", "info")}
    status = "blocked" if counts["fail"] else ("ready_with_warnings" if counts["warn"] else "ready")
    return {
        "schema_version": SCHEMA_VERSION,
        "workflow": "virome_catalogue_v2",
        "config_file": str(config_path),
        "status": status,
        "summary": counts,
        "checks": [asdict(check) for check in diagnostics.checks],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Read-only readiness diagnostic for virome catalogue v2.")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG, help="pipeline.env path (default: project config/pipeline.env)")
    parser.add_argument("--format", choices=("json", "tsv", "text"), default="text")
    parser.add_argument("--write", type=Path, help="Optional path to also write the selected report format.")
    args = parser.parse_args()
    report = run(args.config)
    if args.format == "json":
        payload = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
    elif args.format == "tsv":
        columns = ("check_id", "category", "status", "severity", "message", "remediation", "value")
        rows = ["\t".join(columns)]
        rows.extend("\t".join(str(check.get(column, "")).replace("\t", " ").replace("\n", " ") for column in columns) for check in report["checks"])
        payload = "\n".join(rows) + "\n"
    else:
        summary = report["summary"]
        lines = [f"Virome catalogue v2 readiness: {report['status']}", f"pass={summary['pass']} warn={summary['warn']} fail={summary['fail']}"]
        lines.extend(f"[{check['status'].upper()}] {check['check_id']}: {check['message']}" for check in report["checks"] if check["status"] != "pass")
        payload = "\n".join(lines) + "\n"
    if args.write:
        args.write.parent.mkdir(parents=True, exist_ok=True)
        args.write.write_text(payload, encoding="utf-8")
    sys.stdout.write(payload)
    return 0 if report["status"] != "blocked" else 3


if __name__ == "__main__":
    raise SystemExit(main())
