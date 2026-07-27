from __future__ import annotations

import os
import csv
import io
import json
import shlex
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import streamlit as st
import streamlit.components.v1 as components


PIPELINE_HOME = Path(__file__).resolve().parent
HELPERS_DIR = PIPELINE_HOME / "scripts" / "helpers"
if str(HELPERS_DIR) not in sys.path:
    sys.path.insert(0, str(HELPERS_DIR))
import task_registry
PIPELINE_SCRIPT = PIPELINE_HOME / "scripts" / "run_pipeline.sh"
VIRAL_REPORT_SCRIPT = PIPELINE_HOME / "scripts" / "run_viral_report.sh"
VIROME_CATALOGUE_SCRIPT = PIPELINE_HOME / "scripts" / "run_virome_catalogue.sh"
FINE_ANNOTATION_SCRIPT = PIPELINE_HOME / "scripts" / "run_fine_annotation.sh"
VIROME_DIAGNOSTIC_SCRIPT = PIPELINE_HOME / "scripts" / "diagnose_virome_environment.py"
CONFIG_FILE = PIPELINE_HOME / "config" / "pipeline.env"


def load_settings() -> dict[str, str]:
    values: dict[str, str] = {}
    if CONFIG_FILE.is_file():
        for raw in CONFIG_FILE.read_text(encoding="utf-8", errors="replace").splitlines():
            line = raw.strip()
            if line and not line.startswith("#") and "=" in line:
                key, value = line.split("=", 1)
                values[key.strip()] = value.strip()
    return values


SETTINGS = load_settings()
ALLOWED_ROOTS = [Path(item).resolve() for item in SETTINGS.get("ALLOWED_DATA_ROOTS", "").split(":") if item]
MAX_TOTAL_THREADS = int(SETTINGS.get("MAX_TOTAL_THREADS", "96"))
MAX_VIRAL_THREADS = int(SETTINGS.get("MAX_THREADS_PER_VIRAL_TOOL", "32"))
ADAPTER_CATALOG = Path(
    SETTINGS.get("ADAPTER_CATALOG", str(PIPELINE_HOME / "config" / "adapter_catalog.tsv"))
    .replace("$PIPELINE_HOME", str(PIPELINE_HOME))
)
ADAPTER_SEQUENCE_REFERENCE = Path(
    SETTINGS.get(
        "ADAPTER_SEQUENCE_REFERENCE",
        str(PIPELINE_HOME / "config" / "adapter_sequence_reference.tsv"),
    ).replace("$PIPELINE_HOME", str(PIPELINE_HOME))
)


def load_adapter_profiles() -> list[dict[str, str]]:
    with ADAPTER_CATALOG.open(encoding="utf-8", newline="") as handle:
        return [
            row for row in csv.DictReader(handle, delimiter="\t")
            if row.get("status") in {"active", "review"}
        ]


ADAPTER_PROFILES = load_adapter_profiles()


def load_adapter_reference() -> list[dict[str, str]]:
    with ADAPTER_SEQUENCE_REFERENCE.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


ADAPTER_REFERENCE_ROWS = load_adapter_reference()


def validate_path(text: str, label: str, *, exists: bool, directory: bool = True) -> Path:
    if not text.strip() or not text.strip().startswith("/"):
        raise ValueError(f"{label}必须填写 Linux 服务器上的绝对路径。")
    path = Path(text.strip()).resolve(strict=False)
    if ALLOWED_ROOTS and not any(path == root or root in path.parents for root in ALLOWED_ROOTS):
        raise ValueError(f"{label}不在 ALLOWED_DATA_ROOTS 允许范围内：{path}")
    if exists:
        if directory and not path.is_dir():
            raise ValueError(f"{label}不存在或不是目录：{path}")
        if not directory and not path.is_file():
            raise ValueError(f"{label}不存在或不是文件：{path}")
    elif not path.parent.is_dir():
        raise ValueError(f"{label}的父目录不存在：{path.parent}")
    return path


def launch(command: list[str], task_id: str) -> int:
    process = subprocess.Popen(
        command,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        shell=False,
        start_new_session=True,
        env={**os.environ, "LANG": "C", "LC_ALL": "C", "CONTIG_PIPELINE_TASK_ID": task_id},
    )
    return process.pid


def run_virome_diagnostic() -> tuple[dict[str, object] | None, str | None]:
    """Run the read-only v2 readiness check and return its JSON contract."""
    try:
        result = subprocess.run(
            [sys.executable, str(VIROME_DIAGNOSTIC_SCRIPT), "--config", str(CONFIG_FILE), "--format", "json"],
            cwd=PIPELINE_HOME,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=20,
            check=False,
        )
        report = json.loads(result.stdout)
    except (OSError, subprocess.TimeoutExpired, json.JSONDecodeError) as error:
        return None, f"诊断程序未能生成有效结果：{error}"
    if not isinstance(report, dict) or "checks" not in report:
        return None, "诊断程序返回的结果格式不符合预期。"
    return report, None


def show_virome_diagnostic_panel() -> None:
    key = "virome_environment_diagnostic"
    st.caption("只读检查：不会启动分析、下载数据库或读取你的测序数据。")
    if st.button("检查 v2 环境与参考数据库", key="run_virome_environment_diagnostic", use_container_width=True):
        report, error = run_virome_diagnostic()
        st.session_state[key] = {"report": report, "error": error}
    state = st.session_state.get(key)
    if not state:
        st.info("提交前建议运行一次诊断；它会检查工具、数据库、ICTV 元数据和关键参数。")
        return
    if state["error"]:
        st.error(state["error"])
        return
    report = state["report"]
    status = report["status"]
    summary = report["summary"]
    if status == "ready":
        st.success("v2 运行条件已就绪。仍建议先使用小型验收数据完成一次完整运行。")
    elif status == "ready_with_warnings":
        st.warning("v2 可以启动，但存在建议处理的警告。")
    else:
        st.error("v2 当前不可启动；请先修复下方标为“失败”的检查项。")
    c1, c2, c3 = st.columns(3)
    c1.metric("通过", summary["pass"])
    c2.metric("警告", summary["warn"])
    c3.metric("失败", summary["fail"])
    labels = {"pass": "通过", "warn": "警告", "fail": "失败", "info": "信息"}
    rows = [{
        "状态": labels.get(item["status"], item["status"]),
        "类别": item["category"],
        "检查项": item["check_id"],
        "结果": item["message"],
        "修复建议": item["remediation"],
        "路径/版本": item["value"],
    } for item in report["checks"]]
    st.dataframe(rows, use_container_width=True, hide_index=True)
    st.download_button(
        "下载诊断 JSON",
        data=json.dumps(report, ensure_ascii=False, indent=2),
        file_name="virome_catalogue_v2_readiness.json",
        mime="application/json",
        use_container_width=True,
    )


def latest_run(state_base: Path) -> Path | None:
    runs = state_base / ".contig_pipeline" / "runs"
    if not runs.is_dir():
        return None
    candidates = [item for item in runs.iterdir() if item.is_dir()]
    return max(candidates, key=lambda item: item.stat().st_mtime) if candidates else None


def tail(path: Path, limit: int = 24_000) -> str:
    if not path.is_file():
        return "日志尚未创建。"
    size = path.stat().st_size
    with path.open("rb") as handle:
        handle.seek(max(0, size - limit))
        text = handle.read().decode("utf-8", errors="replace")
    return ("…\n" if size > limit else "") + text


STATUS_LABELS = {"STARTING": "等待建立运行目录", "RUNNING": "运行中", "SUCCESS": "已完成", "FAILED": "失败", "CANCELLED": "已终止"}

MODULES = {
    "preflight": ("输入与样本清单检查", "检查路径、FASTQ/contig 配对和本次样本 manifest。"),
    "fastp": ("fastp 质控与去接头", "生成 cleandata、质控 HTML/JSON 与接头证据。"),
    "megahit": ("MEGAHIT 拼接", "从 clean reads 进行样本级 de novo 拼接。"),
    "check_contigs": ("拼接 contig 检查", "汇总每个样本的拼接数量、长度和基础统计。"),
    "prepare_contigs": ("contig 准备", "仅按本次 manifest 合并 contig，并应用长度阈值。"),
    "genomad": ("geNomad 病毒发现", "基于 geNomad 识别潜在病毒 contig。"),
    "virsorter2": ("VirSorter2 病毒发现", "以 VirSorter2 补充病毒候选识别。"),
    "diamond_virus_discovery": ("DIAMOND 病毒发现", "以 NR 病毒范围的 DIAMOND 命中补充候选。"),
    "build_candidate_catalogue": ("候选序列目录", "合并三种发现证据并构建去冗余 VC catalogue。"),
    "diamond_nr_taxonomy": ("DIAMOND + TaxonKit 分类", "以完整 NR 进行分类比对，并生成 TaxonKit LCA。"),
    "diamond_taxonomy": ("DIAMOND + TaxonKit 分类", "对已有候选执行 NR 分类和 TaxonKit LCA。"),
    "diamond_megan": ("DIAMOND + MEGAN 辅助文件", "生成供本地人工查看的 DAA/RMA6，不阻塞主报告。"),
    "resolve_viral_evidence": ("病毒证据判定", "综合发现证据与 NR 分类，筛选进入质量控制的病毒片段。"),
    "checkv": ("CheckV 质量评估", "评估病毒片段完整度与宿主污染，并保留质量摘要。"),
    "select_ictv_candidates": ("ICTV 候选筛选", "提取已分类到科、适合进入 ICTV 精细注释的片段。"),
    "ictv_refinement": ("ICTV 精细注释", "以本地 ICTV 参考库对已分类到科的片段进行精细比对。"),
    "build_final_catalogue": ("最终病毒片段目录", "整合 CheckV、NR 与 ICTV 信息，形成 VF catalogue。"),
    "quantify_fragments": ("样本分发与丰度", "将最终片段分发回样本并使用 reads 定量。"),
    "votu_abundance": ("vOTU 与丰度", "旧版流程的 vOTU 聚类、样本分发和丰度计算。"),
    "prepare_custom_input": ("自定义候选输入准备", "合并并规范化用户提供的 FASTA 候选序列。"),
    "refresh_main_report": ("主报告刷新", "将独立注释结果回写并刷新旧版主报告。"),
    "custom_annotation_report": ("独立注释报告", "为自定义候选生成独立的注释报告。"),
    "report": ("报告生成", "生成批次总览、单样本页面及可下载结果表。"),
}

WORKFLOW_MODULES = {
    "qc_only": ("preflight", "fastp"),
    "assembly_only": ("preflight", "megahit", "check_contigs"),
    "full": ("preflight", "fastp", "megahit", "check_contigs"),
    "viral_report": ("preflight", "prepare_contigs", "genomad", "virsorter2", "checkv", "votu_abundance", "report"),
    "virome_catalogue": ("preflight", "prepare_contigs", "genomad", "virsorter2", "diamond_virus_discovery", "build_candidate_catalogue", "diamond_nr_taxonomy", "diamond_megan", "resolve_viral_evidence", "checkv", "select_ictv_candidates", "ictv_refinement", "build_final_catalogue", "quantify_fragments", "report"),
    "fine_annotation": ("prepare_custom_input", "diamond_megan", "diamond_taxonomy", "refresh_main_report", "custom_annotation_report"),
}


def task_is_allowed(record: dict[str, object]) -> bool:
    try:
        path = Path(str(record["state_base"])).resolve(strict=False)
    except (KeyError, OSError):
        return False
    return not ALLOWED_ROOTS or any(path == root or root in path.parents for root in ALLOWED_ROOTS)


def task_log_paths(record: dict[str, object], module: str | None = None) -> list[Path]:
    run_text = record.get("run_dir")
    if not run_text or not task_is_allowed(record):
        return []
    run = Path(str(run_text))
    output = Path(str(record["state_base"]))
    paths: list[Path] = []
    if module:
        candidates = [run / f"{module}.log"]
        module_terms = {
            "diamond_virus_discovery": ("02c_diamond_virus",), "diamond_nr_taxonomy": ("04_nr_annotation",),
            "diamond_taxonomy": ("04_nr_annotation", "04_nr_taxonomy"), "diamond_megan": ("04_nr_megan",), "ictv_refinement": ("06_ictv_refinement",),
            "quantify_fragments": ("09_abundance",), "report": ("reports",), "checkv": ("05_checkv", "03_checkv"),
            "genomad": ("02_genomad",), "virsorter2": ("02b_virsorter2",),
            "prepare_contigs": ("01_prepared_contigs",), "fastp": ("fastp_report",), "megahit": ("megahit",),
            "votu_abundance": ("04_sample_votu",), "refresh_main_report": ("reports",),
            "custom_annotation_report": ("reports",), "check_contigs": (".contig_pipeline/reports",),
        }.get(module, ())
        for name in module_terms:
            stage = output / name
            if stage.is_dir():
                candidates.extend(path for path in stage.rglob("*") if path.is_file() and path.suffix in {".log", ".err", ".out"})
                candidates.extend(path for path in stage.glob("background.*") if path.is_file())
                candidates.extend(path for path in stage.glob("*command*.sh") if path.is_file())
        paths = candidates
    else:
        paths = [run / "pipeline.log", run / "parameters.env"]
    seen: set[Path] = set()
    return [path for path in paths if path.is_file() and not (path in seen or seen.add(path))][:36]


def module_parameters(run_dir: Path | None, module: str) -> str:
    if run_dir is None:
        return ""
    values = task_registry.parse_parameters(run_dir / "parameters.env")
    keys = {
        "preflight": ("TASK", "ASSEMBLY_DIR", "CLEAN_DIR", "RAW_DIR"),
        "fastp": ("QC_PARALLEL", "QC_THREADS", "ADAPTER_PROFILE"),
        "megahit": ("ASSEMBLY_PARALLEL", "ASSEMBLY_THREADS", "MIN_CONTIG_LEN"),
        "prepare_contigs": ("MIN_CONTIG_LENGTH",),
        "genomad": ("THREADS", "GENOMAD_DB"), "virsorter2": ("THREADS",),
        "diamond_virus_discovery": ("DIAMOND_THREADS", "DIAMOND_BLOCK_SIZE", "DIAMOND_INDEX_CHUNKS", "DIAMOND_TMPDIR"),
        "diamond_nr_taxonomy": ("DIAMOND_THREADS", "DIAMOND_BLOCK_SIZE", "DIAMOND_INDEX_CHUNKS", "DIAMOND_TMPDIR", "DIAMOND_NR_MAX_TARGET_SEQS"),
        "diamond_taxonomy": ("THREADS", "BLOCK_SIZE", "INDEX_CHUNKS", "TMPDIR", "MAX_TARGET_SEQS", "TAXON_SCOPE", "TAXONLIST"),
        "diamond_megan": ("DIAMOND_THREADS", "DIAMOND_BLOCK_SIZE", "DIAMOND_INDEX_CHUNKS", "DIAMOND_TMPDIR"),
        "checkv": ("THREADS", "CHECKV_DB"), "ictv_refinement": ("DIAMOND_THREADS", "ICTV_REFERENCE_VERSION", "ICTV_REFERENCE_DMND"),
        "quantify_fragments": ("THREADS",), "report": ("GROUPS_FILE",),
    }.get(module, ())
    return "\n".join(f"{key}={values[key]}" for key in keys if values.get(key))


def task_display_name(record: dict[str, object]) -> str:
    stored = str(record.get("display_name") or "").strip()
    if stored:
        return stored
    run_dir = Path(str(record["run_dir"])) if record.get("run_dir") else None
    return task_registry.suggested_display_name(
        str(record["workflow_label"]), Path(str(record["state_base"])), str(record["submitted_at"]), run_dir,
    )


def modules_for_task(record: dict[str, object]) -> list[str]:
    modules = list(WORKFLOW_MODULES.get(str(record["workflow"]), ()))
    run_dir = Path(str(record["run_dir"])) if record.get("run_dir") else None
    parameters = task_registry.parse_parameters(run_dir / "parameters.env") if run_dir else {}
    if record["workflow"] == "fine_annotation":
        if parameters.get("SOURCE") == "checkv":
            modules = [item for item in modules if item not in {"prepare_custom_input", "custom_annotation_report"}]
        elif parameters.get("SOURCE") == "custom":
            modules = [item for item in modules if item != "refresh_main_report"]
    if run_dir and run_dir.is_dir():
        for log in run_dir.glob("*.log"):
            stage = log.stem
            if stage in MODULES and stage not in modules:
                modules.append(stage)
    return modules


def show_task_detail(record: dict[str, object]) -> None:
    st.markdown("<div class='section-title'>任务详情与软件日志</div>", unsafe_allow_html=True)
    status = str(record.get("status", "STARTING"))
    headline = f"{task_display_name(record)} · {STATUS_LABELS.get(status, status)}"
    if status == "RUNNING": st.warning(headline)
    elif status == "SUCCESS": st.success(headline)
    elif status == "FAILED": st.error(headline)
    else: st.info(headline)
    st.caption(f"工作流：{record['workflow_label']} ｜ 提交时间：{record['submitted_at']} ｜ 输出位置：{record['state_base']}")
    rename_col, save_col = st.columns([4, 1])
    with rename_col:
        renamed = st.text_input("任务名称", value=task_display_name(record), key=f"rename_task_{record['task_id']}")
    with save_col:
        st.write("")
        if st.button("保存名称", key=f"save_task_name_{record['task_id']}", use_container_width=True):
            task_registry.rename_task(str(record["task_id"]), renamed)
            st.success("任务名称已保存。")
    if record.get("run_dir"):
        st.caption(f"运行目录：{record['run_dir']} ｜ 当前/最后阶段：{record.get('current_step') or '尚未写入阶段'}")
    else:
        st.caption("后台进程已提交，正在等待后端创建运行目录。")
    if status in {"STARTING", "RUNNING"}:
        with st.expander("终止运行中的任务", expanded=False):
            st.warning("此操作会向该任务及其子进程发送安全终止请求。已有结果文件会保留，但本次任务将标记为“已终止”，不能作为完整结果使用。")
            confirmed = st.checkbox("我已确认要终止这项运行中的任务", key=f"terminate_confirm_{record['task_id']}")
            phrase = st.text_input("请输入“终止”以确认", key=f"terminate_phrase_{record['task_id']}")
            if st.button("确认终止任务", type="primary", disabled=not (confirmed and phrase.strip() == "终止"), key=f"terminate_task_{record['task_id']}", use_container_width=True):
                success, message = task_registry.terminate_task(str(record["task_id"]))
                if success:
                    st.success(message)
                else:
                    st.error(message)
    general_logs = task_log_paths(record)
    if general_logs:
        selected = st.selectbox("总运行记录", general_logs, format_func=lambda item: item.name, key=f"general_log_{record['task_id']}")
        st.code(tail(selected), language="text")
        st.download_button("下载当前日志/参数文件", selected.read_bytes(), file_name=selected.name, use_container_width=True, key=f"download_general_{record['task_id']}")
    st.caption("以下仅显示该工作流相关的软件阶段；状态由该次运行目录的阶段日志和输出文件判定。")
    for module in modules_for_task(record):
        title, explanation = MODULES[module]
        paths = task_log_paths(record, module)
        run_dir = Path(str(record["run_dir"])) if record.get("run_dir") else None
        current = record.get("current_step") == module and status == "RUNNING"
        available = bool(paths)
        module_status = "运行中" if current else ("已产生记录" if available else ("本次未执行或不适用" if status in {"SUCCESS", "FAILED"} else "等待执行"))
        with st.expander(f"{title} · {module_status}"):
            st.caption(explanation)
            parameters = module_parameters(run_dir, module)
            if parameters:
                st.caption("本次任务的有效参数")
                st.code(parameters, language="text")
            if not paths:
                st.info("此任务目前没有该软件的日志或辅助文件。")
                continue
            selected = st.selectbox("查看该模块文件", paths, format_func=lambda item: str(item.relative_to(Path(str(record['state_base'])))) if Path(str(record['state_base'])) in item.parents else item.name, key=f"module_log_{record['task_id']}_{module}")
            st.code(tail(selected), language="text")
            st.download_button("下载此文件", selected.read_bytes(), file_name=selected.name, use_container_width=True, key=f"download_module_{record['task_id']}_{module}")


def show_task_center() -> None:
    st.markdown("<div class='section-title'>任务历史中心</div>", unsafe_allow_html=True)
    scan_key = "task_history_auto_discovery_complete"
    force_scan = st.button("刷新任务状态与历史记录", use_container_width=True, key="refresh_tasks")
    if ALLOWED_ROOTS and (force_scan or not st.session_state.get(scan_key)):
        with st.spinner("正在从已配置的数据目录发现历史任务（仅读取目录名、状态与日志元数据）…"):
            imported, locations, truncated = task_registry.discover_history(ALLOWED_ROOTS)
        st.session_state[scan_key] = True
        message = f"已扫描 {locations} 个含运行记录的输出位置；新增登记 {imported} 个任务。"
        st.caption(message + (" 扫描达到安全上限，未继续遍历更深目录。" if truncated else ""))
    elif not ALLOWED_ROOTS:
        st.info("配置 ALLOWED_DATA_ROOTS 后，任务历史中心会自动发现其中的历史运行记录。")
    records = [record for record in task_registry.refresh_tasks() if task_is_allowed(record)]
    categories = [("全部", records), ("运行中", [r for r in records if r["status"] in {"STARTING", "RUNNING"}]), ("失败", [r for r in records if r["status"] == "FAILED"]), ("已完成", [r for r in records if r["status"] == "SUCCESS"]), ("已终止", [r for r in records if r["status"] == "CANCELLED"])]
    tabs = st.tabs([item[0] for item in categories])

    def date_group(record: dict[str, object]) -> str:
        try:
            submitted = datetime.fromisoformat(str(record["submitted_at"]).replace("Z", "+00:00"))
            local_date = submitted.astimezone(ZoneInfo("Asia/Shanghai")).date()
            today = datetime.now(ZoneInfo("Asia/Shanghai")).date()
        except (TypeError, ValueError):
            return "日期未知"
        if local_date == today:
            return f"今天 · {local_date.isoformat()}"
        if (today - local_date).days == 1:
            return f"昨天 · {local_date.isoformat()}"
        return local_date.isoformat()

    def show_record(record: dict[str, object], category: str) -> None:
        left, middle, right = st.columns([4, 3, 1])
        left.write(f"**{task_display_name(record)}**")
        run_name = Path(str(record["run_dir"])).name if record.get("run_dir") else "等待建立运行目录"
        middle.caption(f"{record['workflow_label']} · {STATUS_LABELS.get(record['status'], record['status'])} · {run_name}")
        right.button("进入工作流", key=f"open_task_{category}_{record['task_id']}", use_container_width=True,
                     on_click=open_task_from_center, args=(record,))

    for tab, (category, subset) in zip(tabs, categories):
        with tab:
            if not subset:
                st.caption("暂无任务记录。")
                continue
            if category == "运行中":
                for record in subset[:100]:
                    show_record(record, category)
                continue
            grouped: dict[str, list[dict[str, object]]] = {}
            for record in subset[:200]:
                grouped.setdefault(date_group(record), []).append(record)
            for group_name, day_records in grouped.items():
                with st.expander(f"{group_name}（{len(day_records)} 条任务）", expanded=group_name.startswith("今天")):
                    for record in day_records:
                        show_record(record, category)


def show_selected_task_for_workflow(workflow: str) -> None:
    selected_id = st.session_state.get("selected_task_id")
    if not selected_id:
        return
    records = [record for record in task_registry.refresh_tasks() if task_is_allowed(record)]
    record = next((item for item in records if item["task_id"] == selected_id), None)
    if record is None:
        return
    related = {workflow}
    if workflow == "virome_catalogue":
        related.add("viral_report")
    if record["workflow"] not in related:
        return
    show_task_detail(record)
    if record["workflow"] in {"virome_catalogue", "viral_report"}:
        show_report_center(Path(str(record["state_base"])))


def show_adapter_evidence(clean_root: Path | None) -> None:
    if clean_root is None or not clean_root.is_dir():
        return
    evidence_files = sorted(clean_root.glob("*/fastp_report/*.adapter_evidence.tsv"))
    scan_files = sorted(clean_root.glob("*/fastp_report/*.adapter_reference_scan.tsv"))
    if evidence_files or scan_files:
        st.markdown("<div class='section-title'>接头识别与来源判断</div>", unsafe_allow_html=True)
    if evidence_files:
        rows: list[dict[str, str]] = []
        for evidence_file in evidence_files:
            with evidence_file.open(encoding="utf-8", newline="") as handle:
                rows.extend(csv.DictReader(handle, delimiter="\t"))
        st.caption("“自动识别”来自 fastp 的双端重叠/PE识别；手动方案表示按所选建库序列剪切。来源判断不等同于仅凭序列确定建库试剂盒。")
        visible = [{
            "样本": row["sample_id"], "读段": row["read"],
            "fastp报告序列": row["fastp_reported_sequence"],
            "目录匹配": row["matched_catalogue_profiles"],
            "参考序列匹配": row.get("matched_reference_sequences", "旧版证据表未提供"),
            "判断": row["source_judgement"], "剪切reads": row["adapter_trimmed_reads"],
            "剪切比例": row["trimmed_read_fraction"], "来源等级": row["source_level"],
            "来源": row["source_title"],
        } for row in rows]
        st.dataframe(visible, use_container_width=True, hide_index=True)
        buffer = io.StringIO()
        writer = csv.DictWriter(buffer, fieldnames=rows[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
        st.download_button("下载接头识别证据表", data=buffer.getvalue(), file_name="adapter_evidence.tsv", mime="text/tab-separated-values")
    if scan_files:
        scan_rows: list[dict[str, str]] = []
        for scan_file in scan_files:
            with scan_file.open(encoding="utf-8", newline="") as handle:
                scan_rows.extend(csv.DictReader(handle, delimiter="\t"))
        if scan_rows:
            st.markdown("**原始 reads 参考序列抽样命中**")
            st.caption("精确匹配用于发现SISPA标签、完整PCR引物和flow-cell伪影；参考或拒绝条目不会自动参与剪切。")
            st.dataframe(scan_rows, use_container_width=True, hide_index=True)


def show_report_center(report_root: Path | None) -> None:
    """Embed the self-contained offline dashboard without exposing a file URL."""
    if report_root is None:
        return
    dashboard = report_root / "reports" / "virome_catalogue_dashboard.html"
    if not dashboard.is_file():
        dashboard = report_root / "reports" / "virome_dashboard.html"
    if not dashboard.is_file():
        st.info("报告完成后，这里将显示“批次总览”和可切换的单样本判读页面。")
        return
    st.markdown("<div class='section-title'>报告中心</div>", unsafe_allow_html=True)
    st.caption(f"离线交互式报告：{dashboard}")
    report_html = dashboard.read_text(encoding="utf-8", errors="replace")
    key = f"dashboard_open_{abs(hash(str(dashboard)))}"
    col_open, col_download = st.columns([1, 1])
    with col_open:
        if st.button("在网页内打开交互式报告", use_container_width=True):
            st.session_state[key] = not st.session_state.get(key, False)
    with col_download:
        st.download_button(
            "下载离线 HTML 报告",
            data=report_html,
            file_name="virome_dashboard.html",
            mime="text/html",
            use_container_width=True,
        )
    if st.session_state.get(key, False):
        components.html(report_html, height=1450, scrolling=True)


def taxon_scope_ui() -> tuple[str, str]:
    scope_label = st.radio(
        "NR 检索范围",
        ["仅病毒（默认，NCBI TaxID 10239）", "不限制分类范围（检索完整 NR）", "自定义 NCBI TaxID"],
        help="这是 DIAMOND 的数据库检索过滤条件，不是更换数据库。三个选项都使用管理员配置的完整 NR 数据库。",
    )
    if scope_label.startswith("仅病毒"):
        return "virus", SETTINGS.get("DIAMOND_DEFAULT_TAXONLIST", "10239")
    if scope_label.startswith("不限制"):
        return "none", ""
    ids = st.text_input("自定义 TaxID 列表", placeholder="例如：10239,2157", help="填写一个或多个 NCBI TaxID，使用英文逗号分隔。")
    return "custom", ids.strip()


st.set_page_config(page_title="Virome Contig Studio", page_icon="🧬", layout="wide")
st.markdown(
    """<style>
    .block-container{max-width:1420px;padding-top:1.6rem;padding-bottom:3rem}
    [data-testid='stAppViewContainer']{background:#f5f8fb;color:#172b3a}
    [data-testid='stSidebar']{background:#102a43}
    [data-testid='stSidebar'] *{color:#edf6f7}
    h1{color:#0b3954!important;letter-spacing:-.03em}.hero{background:linear-gradient(120deg,#0b3954,#146b7a);border-radius:18px;padding:24px 28px;color:#fff;margin-bottom:20px}.hero h2{color:#fff;margin:0 0 6px 0}.hero p{margin:0;color:#d8f0ee}.section-title{font-size:1.15rem;font-weight:700;color:#0b3954;margin:1.35rem 0 .55rem}.task-note{padding:.4rem 0 .6rem;color:#466375}.status-card{background:#fff;border:1px solid #d9e6ec;border-radius:12px;padding:12px 14px}.stButton>button[kind='primary']{background:#0f897f;border-color:#0f897f}.stButton>button[kind='primary']:hover{background:#0a7069;border-color:#0a7069}.stRadio>div{gap:.35rem}.stRadio label{padding:.22rem .35rem}.stNumberInput label,.stTextInput label{font-weight:600}
    </style>""",
    unsafe_allow_html=True,
)
st.markdown("""<div class='hero'><h2>🧬 Virome Contig Studio</h2><p>从原始测序数据到病毒多样性解读与 DIAMOND 精细注释的中文科研工作台。</p></div>""", unsafe_allow_html=True)

if not all(script.is_file() for script in [PIPELINE_SCRIPT, VIRAL_REPORT_SCRIPT, VIROME_CATALOGUE_SCRIPT, FINE_ANNOTATION_SCRIPT, VIROME_DIAGNOSTIC_SCRIPT]):
    st.error("未找到必需脚本。请重新同步完整的 contig_pipeline 软件目录。")
    st.stop()
if not ALLOWED_ROOTS:
    st.warning("尚未配置路径白名单。请先复制 config/pipeline.env.example 为 config/pipeline.env，并设置 ALLOWED_DATA_ROOTS。")

TASKS = {
    "① 原始数据质控": "qc_only",
    "② MEGAHIT 拼接": "assembly_only",
    "③ 完整拼接流程": "full",
    "④ 病毒发现注释": "virome_catalogue",
    "⑤ 任务历史中心": "task_center",
}
VIROME_MAIN_OPERATION = "病毒发现注释（v2 全流程）"
VIROME_DIAMOND_OPERATION = "对已有候选进行独立 DIAMOND 精细注释"


def open_task_from_center(record: dict[str, object]) -> None:
    """Open a historical task in the workflow that created it."""
    workflow = str(record["workflow"])
    target = {
        "qc_only": "① 原始数据质控",
        "assembly_only": "② MEGAHIT 拼接",
        "full": "③ 完整拼接流程",
        "virome_catalogue": "④ 病毒发现注释",
        "viral_report": "④ 病毒发现注释",
        "fine_annotation": "④ 病毒发现注释",
    }.get(workflow, "⑤ 任务历史中心")
    st.session_state["selected_task_id"] = record["task_id"]
    st.session_state["workflow_navigation"] = target
    if workflow in {"fine_annotation"}:
        st.session_state["virome_operation"] = VIROME_DIAMOND_OPERATION
    elif workflow in {"virome_catalogue", "viral_report"}:
        st.session_state["virome_operation"] = VIROME_MAIN_OPERATION


with st.sidebar:
    st.markdown("## 工作流导航")
    chosen = st.radio("选择任务", list(TASKS), label_visibility="collapsed", key="workflow_navigation")
    task = TASKS[chosen]
    st.markdown("---")
    st.caption("网页仅调度固定的后端脚本；不会执行任意 Shell 命令。")
    st.caption(f"线程上限：{MAX_TOTAL_THREADS}；病毒工具上限：{MAX_VIRAL_THREADS}")

if task == "task_center":
    show_task_center()
    st.stop()

if task == "virome_catalogue":
    operation = st.radio(
        "第④工作流中的操作",
        [VIROME_MAIN_OPERATION, VIROME_DIAMOND_OPERATION],
        horizontal=True,
        key="virome_operation",
    )
    if operation == VIROME_DIAMOND_OPERATION:
        task = "fine_annotation"

st.markdown(f"<div class='section-title'>{chosen}</div>", unsafe_allow_html=True)
task_copy = {
    "qc_only": "将双端 rawdata 进行 fastp 质控，生成标准化 cleandata。",
    "assembly_only": "从已有 cleandata 启动 PE 或 SE MEGAHIT 拼接。",
    "full": "顺序执行 fastp 质控、MEGAHIT 拼接与 contig 检查。",
    "virome_catalogue": "以 geNomad、VirSorter2 和 DIAMOND-NR-virus 建立全局潜在病毒序列池；随后进行完整 NR 分类、全局 CheckV、ICTV 精细注释、样本分发与 reads 定量。",
    "fine_annotation": "对 CheckV 默认候选或自定义候选 contigs，追加 NR DIAMOND、MEGAN RMA6 和 TaxonKit LCA 注释。",
}
st.markdown(f"<div class='task-note'>{task_copy[task]}</div>", unsafe_allow_html=True)
show_selected_task_for_workflow(task)

layout_labels = {"每个样本一个子文件夹": "sample_subdirs", "所有 FASTQ 位于同一文件夹": "flat"}
state_base: Path | None = None
adapter_evidence_root: Path | None = None
command: list[str] | None = None

try:
    if task in {"qc_only", "assembly_only", "full"}:
        needs_raw = task in {"qc_only", "full"}
        needs_assembly = task in {"assembly_only", "full"}
        raw_text = clean_text = assembly_text = ""
        raw_layout = clean_layout = "sample_subdirs"
        read_type = "pe"
        if needs_raw:
            with st.container(border=True):
                st.subheader("输入：rawdata")
                raw_text = st.text_input("rawdata 路径", placeholder="/data/project/rawdata")
                raw_layout = layout_labels[st.radio("测序文件存放方式", list(layout_labels), horizontal=True)]
        with st.container(border=True):
            st.subheader("cleandata")
            clean_text = st.text_input("cleandata 路径", placeholder="/data/project/cleandata", help="质控任务中为输出路径；拼接任务中为输入路径。")
            if task == "assembly_only":
                clean_layout = layout_labels[st.radio("clean reads 存放方式", list(layout_labels), horizontal=True)]
                read_type = {"双端 PE": "pe", "单端 SE": "se"}[st.radio("测序类型", ["双端 PE", "单端 SE"], horizontal=True)]
            elif task == "full":
                st.caption("完整流程固定为双端 PE；fastp 输出会自动采用每样本一个子文件夹的标准布局。")
        if needs_assembly:
            with st.container(border=True):
                st.subheader("输出：assembly")
                assembly_text = st.text_input("assembly 输出路径", placeholder="/data/project/assembly")
        st.markdown("<div class='section-title'>资源设置</div>", unsafe_allow_html=True)
        c1, c2, c3, c4 = st.columns(4)
        with c1: qc_parallel = st.number_input("fastp 并发样本数", 1, 64, 4, disabled=not needs_raw)
        with c2: qc_threads = st.number_input("fastp 每样本线程", 1, 128, 8, disabled=not needs_raw)
        with c3: asm_parallel = st.number_input("MEGAHIT 并发样本数", 1, 32, 2, disabled=not needs_assembly)
        with c4: asm_threads = st.number_input("MEGAHIT 每样本线程", 1, 256, 30, disabled=not needs_assembly)
        min_len = st.number_input("最短 contig 长度", 1, 100000, 40, disabled=not needs_assembly)
        adapter_profile = SETTINGS.get("FASTP_DEFAULT_ADAPTER_PROFILE", "auto")
        if needs_raw:
            st.markdown("<div class='section-title'>去接头策略</div>", unsafe_allow_html=True)
            profile_by_label = {
                f"{row['display_name']}｜{row['source_level']}": row for row in ADAPTER_PROFILES
            }
            selected_label = st.selectbox(
                "建库接头方案",
                list(profile_by_label),
                help="建库类型未知时请选择自动识别。手动方案只应在实验记录能够确认建库试剂盒时使用。",
            )
            selected_adapter = profile_by_label[selected_label]
            adapter_profile = selected_adapter["profile_id"]
            if selected_adapter["r1_sequence"]:
                st.code(
                    f"R1: {selected_adapter['r1_sequence']}\n"
                    f"R2: {selected_adapter['r2_sequence']}",
                    language="text",
                )
            st.caption(
                f"来源：{selected_adapter['source_title']}；核验日期："
                f"{selected_adapter['verified_on']}；状态：{selected_adapter['status']}"
            )
            with st.expander("查看接头目录与维护信息"):
                st.dataframe(ADAPTER_PROFILES, use_container_width=True, hide_index=True)
                st.markdown("**识别参考序列（不会自动作为剪切方案）**")
                st.dataframe(ADAPTER_REFERENCE_ROWS, use_container_width=True, hide_index=True)
        st.info(f"fastp 峰值：{qc_parallel * qc_threads if needs_raw else 0} 线程；MEGAHIT 峰值：{asm_parallel * asm_threads if needs_assembly else 0} 线程。")
        raw = validate_path(raw_text, "rawdata 路径", exists=True) if needs_raw else None
        clean = validate_path(clean_text, "cleandata 路径", exists=task == "assembly_only")
        assembly = validate_path(assembly_text, "assembly 输出路径", exists=False) if needs_assembly else None
        if needs_raw and qc_parallel * qc_threads > MAX_TOTAL_THREADS: raise ValueError("fastp 总线程数超过配置上限。")
        if needs_assembly and asm_parallel * asm_threads > MAX_TOTAL_THREADS: raise ValueError("MEGAHIT 总线程数超过配置上限。")
        command = [str(PIPELINE_SCRIPT), "--task", task, "--cleandata-dir", str(clean), "--qc-parallel", str(qc_parallel), "--qc-threads", str(qc_threads), "--assembly-parallel", str(asm_parallel), "--assembly-threads", str(asm_threads), "--adapter-profile", adapter_profile, "--min-contig-len", str(min_len)]
        if raw: command += ["--rawdata-dir", str(raw), "--raw-layout", raw_layout]
        if assembly: command += ["--clean-layout", clean_layout, "--read-type", read_type, "--assembly-dir", str(assembly)]
        state_base = assembly if assembly else clean
        adapter_evidence_root = clean if needs_raw else None

    elif task == "virome_catalogue":
        with st.expander("运行前：环境与参考数据库诊断", expanded=True):
            show_virome_diagnostic_panel()
        with st.container(border=True):
            st.subheader("输入：已有 cleandata 与 assembly")
            assembly_text = st.text_input("assembly 路径", placeholder="/data/project/assembly")
            clean_text = st.text_input("cleandata 路径", placeholder="/data/project/cleandata")
            clean_layout = layout_labels[st.radio("clean reads 存放方式", list(layout_labels), horizontal=True)]
            read_type = {"双端 PE": "pe", "单端 SE": "se"}[st.radio("测序类型", ["双端 PE", "单端 SE"], horizontal=True)]
        with st.container(border=True):
            st.subheader("输出：全局病毒目录报告")
            output_text = st.text_input("virome_catalogue 输出路径", placeholder="/data/project/virome_catalogue")
            st.caption("先建立全局完全去冗余 VC 序列目录，再生成 CheckV 修正后的 VF 病毒片段；VC/VF 不是 vOTU，也不直接等同于病毒物种。")
        c1, c2 = st.columns(2)
        with c1: threads = st.number_input("病毒分析线程数", 1, MAX_VIRAL_THREADS, min(8, MAX_VIRAL_THREADS))
        with c2: viral_min = st.number_input("进入三工具发现的最短 contig 长度", 200, 100000, max(200, int(SETTINGS.get("VIRAL_MIN_CONTIG_LEN", "1000"))))
        with st.expander("DIAMOND 性能参数", expanded=True):
            st.caption("这些参数同时用于病毒发现、全 NR 分类、后台 DAA 生成和 ICTV 精细比对；每次运行都会写入参数与 DIAMOND 命令审计文件。")
            d1, d2, d3, d4 = st.columns(4)
            with d1:
                diamond_threads = st.number_input("DIAMOND 每任务线程数", 1, MAX_VIRAL_THREADS, min(MAX_VIRAL_THREADS, int(SETTINGS.get("DIAMOND_THREADS_PER_JOB", "64"))), key="v2_diamond_threads")
            with d2:
                block_options = ["0.5", "1", "2", "4", "6", "8", "12", "16"]
                configured_block = SETTINGS.get("DIAMOND_BLOCK_SIZE", "4.0").rstrip("0").rstrip(".") or "4"
                diamond_block_size = st.selectbox("--block-size（GB）", block_options, index=block_options.index(configured_block) if configured_block in block_options else block_options.index("4"), key="v2_diamond_block_size")
            with d3:
                chunk_options = [1, 2, 4, 8]
                configured_chunks = int(SETTINGS.get("DIAMOND_INDEX_CHUNKS", "1"))
                diamond_index_chunks = st.selectbox("--index-chunks", chunk_options, index=chunk_options.index(configured_chunks) if configured_chunks in chunk_options else 0, key="v2_diamond_index_chunks")
            with d4:
                diamond_tmpdir = st.text_input("DIAMOND 临时目录（-t）", SETTINGS.get("DIAMOND_TMPDIR", "/dev/shm"), key="v2_diamond_tmpdir", help="推荐 /dev/shm；必须存在、可写且具有足够可用空间。")
        with st.expander("查看 v2 分析结构与有效参数"):
            st.dataframe([
                {"阶段": "发现", "参数": "geNomad + VirSorter2 + DIAMOND-NR-virus", "值": "均启用", "来源": "v2 固定流程"},
                {"阶段": "分类", "参数": "DIAMOND", "值": "完整 NR（无 TaxID 限制）", "来源": "v2 固定流程"},
                {"阶段": "质量", "参数": "CheckV", "值": "全局候选目录运行一次", "来源": "v2 固定流程"},
                {"阶段": "精细注释", "参数": "ICTV 本地参考库", "值": SETTINGS.get("ICTV_REFERENCE_VERSION", "未配置"), "来源": "pipeline.env"},
                {"参数": "CoverM read identity (%)", "值": SETTINGS.get("COVERM_MIN_READ_PERCENT_IDENTITY", "95"), "来源": "pipeline.env"},
                {"参数": "CoverM read aligned (%)", "值": SETTINGS.get("COVERM_MIN_READ_ALIGNED_PERCENT", "75"), "来源": "pipeline.env"},
                {"参数": "CoverM covered fraction (%)", "值": SETTINGS.get("COVERM_MIN_COVERED_FRACTION", "10"), "来源": "pipeline.env"},
            ], use_container_width=True, hide_index=True)
        required_v2 = ["GENOMAD_DB", "CHECKV_DB", "DIAMOND_NR_DB", "TAXONKIT_DB", "MEGAN_DAA2RMA", "MEGAN_MAP_DB", "ICTV_REFERENCE_DMND", "ICTV_REFERENCE_METADATA"]
        missing_v2 = [key for key in required_v2 if not SETTINGS.get(key)]
        if missing_v2:
            st.warning("v2 尚未配置：" + "、".join(missing_v2) + "。请完成 config/pipeline.env 后再提交。")
        assembly = validate_path(assembly_text, "assembly 路径", exists=True)
        clean = validate_path(clean_text, "cleandata 路径", exists=True)
        output = validate_path(output_text, "virome_catalogue 输出路径", exists=False)
        command = [str(VIROME_CATALOGUE_SCRIPT), "--assembly-dir", str(assembly), "--cleandata-dir", str(clean), "--clean-layout", clean_layout, "--read-type", read_type, "--output-dir", str(output), "--threads", str(threads), "--diamond-threads", str(diamond_threads), "--diamond-block-size", diamond_block_size, "--diamond-index-chunks", str(diamond_index_chunks), "--diamond-tmpdir", diamond_tmpdir.strip(), "--min-contig-length", str(viral_min)]
        state_base = output

    else:
        with st.container(border=True):
            st.subheader("① 候选 contigs 来源")
            source_label = st.radio("输入来源", ["使用已有 viral_report 中 CheckV 筛选后的病毒候选 contigs（默认）", "使用自定义候选 contigs"], help="默认来源可回写主病毒多样性报告；自定义来源生成独立补充报告，不会与已有 vOTU 丰度强行合并。")
            source = "checkv" if source_label.startswith("使用已有") else "custom"
            if source == "checkv":
                report_root_text = st.text_input("已有 viral_report 路径", placeholder="/data/project/viral_report")
                custom_input = custom_type = custom_output = None
                st.caption("输入固定为：viral_report/03_checkv/viral_candidates_checkv.fna")
            else:
                custom_type = {"单个 FASTA 文件（推荐）": "file", "文件夹中的多个 FASTA 文件，自动合并": "directory"}[st.radio("自定义输入形式", ["单个 FASTA 文件（推荐）", "文件夹中的多个 FASTA 文件，自动合并"], horizontal=True)]
                custom_input_text = st.text_input("自定义候选 contigs 路径", placeholder="/data/project/candidates.fna")
                custom_output_text = st.text_input("独立注释输出路径", placeholder="/data/project/custom_annotation_001")
                report_root = None
        with st.container(border=True):
            st.subheader("② DIAMOND 检索策略")
            scope, taxonlist = taxon_scope_ui()
            mode = {"同时生成 DAA/RMA6 和 outfmt 6/LCA（推荐）": "both", "仅生成 DAA 与 MEGAN RMA6": "megan", "仅生成 outfmt 6 与 TaxonKit LCA": "taxonomy"}[st.radio("输出内容", ["同时生成 DAA/RMA6 和 outfmt 6/LCA（推荐）", "仅生成 DAA 与 MEGAN RMA6", "仅生成 outfmt 6 与 TaxonKit LCA"], horizontal=True)]
            st.caption("默认使用完整 NR 数据库并以 TaxID 10239 限制至病毒。选择“完整 NR”时不传递 --taxonlist，但仍可通过 TaxonKit 生成 LCA。")
        with st.container(border=True):
            st.subheader("③ 资源与可追溯性")
            threads = st.number_input("DIAMOND 每任务线程数", 1, MAX_VIRAL_THREADS, min(MAX_VIRAL_THREADS, int(SETTINGS.get("DIAMOND_THREADS_PER_JOB", "64"))), key="annotation_diamond_threads")
            max_hits = st.number_input("每条 contig 最大保留命中数", 1, 200, int(SETTINGS.get("DIAMOND_NR_MAX_TARGET_SEQS", "25")), help="TaxonKit LCA 需要保留多个命中；不建议设为 1。")
            d1, d2, d3 = st.columns(3)
            with d1:
                block_options = ["0.5", "1", "2", "4", "6", "8", "12", "16"]
                configured_block = SETTINGS.get("DIAMOND_BLOCK_SIZE", "4.0").rstrip("0").rstrip(".") or "4"
                diamond_block_size = st.selectbox("--block-size（GB）", block_options, index=block_options.index(configured_block) if configured_block in block_options else block_options.index("4"), key="annotation_diamond_block_size")
            with d2:
                chunk_options = [1, 2, 4, 8]
                configured_chunks = int(SETTINGS.get("DIAMOND_INDEX_CHUNKS", "1"))
                diamond_index_chunks = st.selectbox("--index-chunks", chunk_options, index=chunk_options.index(configured_chunks) if configured_chunks in chunk_options else 0, key="annotation_diamond_index_chunks")
            with d3:
                diamond_tmpdir = st.text_input("DIAMOND 临时目录（-t）", SETTINGS.get("DIAMOND_TMPDIR", "/dev/shm"), key="annotation_diamond_tmpdir", help="推荐 /dev/shm；必须存在且可写。")
            if mode in {"megan", "both"} and (not SETTINGS.get("MEGAN_DAA2RMA") or not SETTINGS.get("MEGAN_MAP_DB")):
                st.warning("选择 RMA6 时需要在配置中填写 MEGAN_DAA2RMA 与 MEGAN_MAP_DB。")
            if mode in {"taxonomy", "both"} and not SETTINGS.get("TAXONKIT_DB"):
                st.warning("选择 LCA 时需要在配置中填写 TAXONKIT_DB。")
        if scope == "custom" and not taxonlist:
            raise ValueError("自定义 NCBI TaxID 不能为空。")
        command = [str(FINE_ANNOTATION_SCRIPT), "--source", source, "--mode", mode, "--threads", str(threads), "--taxon-scope", scope, "--max-target-seqs", str(max_hits), "--block-size", diamond_block_size, "--index-chunks", str(diamond_index_chunks), "--tmpdir", diamond_tmpdir.strip()]
        if taxonlist: command += ["--taxonlist", taxonlist]
        if source == "checkv":
            report_root = validate_path(report_root_text, "已有 viral_report 路径", exists=True)
            command += ["--viral-report-dir", str(report_root)]
            state_base = report_root
        else:
            custom_input = validate_path(custom_input_text, "自定义候选 contigs 路径", exists=True, directory=custom_type == "directory")
            custom_output = validate_path(custom_output_text, "独立注释输出路径", exists=False)
            command += ["--custom-input", str(custom_input), "--custom-input-type", str(custom_type), "--output-dir", str(custom_output)]
            state_base = custom_output

    resume = st.checkbox("跳过已完整完成的步骤（resume）", value=True)
    st.caption("若改变影响结果的参数，已有输出只有在参数契约一致时才会续跑；不一致时请使用新的输出目录。")
    if resume and command is not None:
        command.append("--resume")
    if command is not None:
        st.markdown("<div class='section-title'>执行确认</div>", unsafe_allow_html=True)
        workflow_label = chosen if task != "fine_annotation" else f"{chosen} · 独立 DIAMOND 精细注释"
        inferred_batch = state_base.parent.name or state_base.name
        inferred_task = state_base.name or str(state_base)
        st.caption("任务名称按“样本批次 · 执行任务 · 运行性质”组合；前两项默认按输出路径的倒数两层识别。")
        n1, n2, n3 = st.columns(3)
        with n1:
            batch_choice = st.selectbox("样本批次", [inferred_batch, "手动填写…"], key=f"task_batch_choice_{task}")
            batch_name = st.text_input("自定义样本批次", key=f"task_batch_custom_{task}") if batch_choice == "手动填写…" else batch_choice
        with n2:
            task_choice = st.selectbox("执行任务", [inferred_task, "手动填写…"], key=f"task_operation_choice_{task}")
            operation_name = st.text_input("自定义执行任务", key=f"task_operation_custom_{task}") if task_choice == "手动填写…" else task_choice
        with n3:
            run_kind = st.selectbox("运行性质", ["首次运行", "优化重跑"], key=f"task_run_kind_{task}")
        if not batch_name.strip() or not operation_name.strip():
            raise ValueError("手动命名时，“样本批次”和“执行任务”不能为空。")
        name_parts = [batch_name.strip(), operation_name.strip(), run_kind]
        extra_note = st.text_input("补充说明（可选）", placeholder="例如：提高 DIAMOND block size", key=f"task_note_{task}")
        if extra_note.strip():
            name_parts.append(extra_note.strip())
        task_name = " · ".join(name_parts)
        st.info(f"任务名称预览：{task_name}")
        st.code(" \\\n  ".join(shlex.quote(part) for part in command), language="bash")
        c_run, c_refresh = st.columns([1, 1])
        with c_run:
            if st.button("提交后台任务", type="primary", use_container_width=True):
                task_id = task_registry.register_submission(
                    task, workflow_label, state_base, command, 0,
                    display_name=task_name,
                )
                try:
                    pid = launch(command, task_id)
                except OSError as error:
                    with task_registry.connect() as connection:
                        connection.execute("UPDATE tasks SET status='FAILED', completed_at=? WHERE task_id=?", (task_registry.now(), task_id))
                    st.error(f"后台任务未能启动：{error}")
                else:
                    # The process PID becomes available only after Popen succeeds.
                    # Re-registering is intentionally avoided: update the existing row instead.
                    with task_registry.connect() as connection:
                        connection.execute("UPDATE tasks SET pid=? WHERE task_id=?", (pid, task_id))
                    st.session_state["selected_task_id"] = task_id
                    st.success(f"任务已提交至 Linux 后台（启动 PID：{pid}）。已登记到任务历史中心，可随时切换页面后回访。")
        with c_refresh:
            st.button("刷新运行中心", use_container_width=True)
except ValueError as error:
    st.info(str(error))

show_adapter_evidence(adapter_evidence_root)
