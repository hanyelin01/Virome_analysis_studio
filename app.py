from __future__ import annotations

import os
import csv
import io
import shlex
import subprocess
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components


PIPELINE_HOME = Path(__file__).resolve().parent
PIPELINE_SCRIPT = PIPELINE_HOME / "scripts" / "run_pipeline.sh"
VIRAL_REPORT_SCRIPT = PIPELINE_HOME / "scripts" / "run_viral_report.sh"
FINE_ANNOTATION_SCRIPT = PIPELINE_HOME / "scripts" / "run_fine_annotation.sh"
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


def launch(command: list[str]) -> int:
    process = subprocess.Popen(
        command,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        shell=False,
        start_new_session=True,
        env={**os.environ, "LANG": "C", "LC_ALL": "C"},
    )
    return process.pid


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


def show_status(state_base: Path | None) -> None:
    st.markdown("<div class='section-title'>运行中心</div>", unsafe_allow_html=True)
    if state_base is None:
        st.info("填写有效输出路径后，可在这里查看任务状态与实时日志尾部。")
        return
    run = latest_run(state_base)
    if run is None:
        st.info("当前输出位置还没有提交过任务。")
        return
    status_file = run / "status"
    status = status_file.read_text(encoding="utf-8", errors="replace").strip() if status_file.is_file() else "STARTING"
    if status == "RUNNING":
        st.warning(f"正在运行：{run.name}")
    elif status == "SUCCESS":
        st.success(f"已完成：{run.name}")
    elif status == "FAILED":
        st.error(f"运行失败：{run.name}")
    else:
        st.info(f"状态：{status or 'STARTING'}；运行编号：{run.name}")
    st.caption(f"运行目录：{run}")
    st.code(tail(run / "pipeline.log"), language="text")


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

if not all(script.is_file() for script in [PIPELINE_SCRIPT, VIRAL_REPORT_SCRIPT, FINE_ANNOTATION_SCRIPT]):
    st.error("未找到必需脚本。请重新同步完整的 contig_pipeline 软件目录。")
    st.stop()
if not ALLOWED_ROOTS:
    st.warning("尚未配置路径白名单。请先复制 config/pipeline.env.example 为 config/pipeline.env，并设置 ALLOWED_DATA_ROOTS。")

TASKS = {
    "① 原始数据质控": "qc_only",
    "② MEGAHIT 拼接": "assembly_only",
    "③ 完整拼接流程": "full",
    "④ 逐样本病毒报告与批次总览": "viral_report",
    "⑤ DIAMOND 精细注释": "fine_annotation",
}
with st.sidebar:
    st.markdown("## 工作流导航")
    chosen = st.radio("选择任务", list(TASKS), label_visibility="collapsed")
    task = TASKS[chosen]
    st.markdown("---")
    st.caption("网页仅调度固定的后端脚本；不会执行任意 Shell 命令。")
    st.caption(f"线程上限：{MAX_TOTAL_THREADS}；病毒工具上限：{MAX_VIRAL_THREADS}")

st.markdown(f"<div class='section-title'>{chosen}</div>", unsafe_allow_html=True)
task_copy = {
    "qc_only": "将双端 rawdata 进行 fastp 质控，生成标准化 cleandata。",
    "assembly_only": "从已有 cleandata 启动 PE 或 SE MEGAHIT 拼接。",
    "full": "顺序执行 fastp 质控、MEGAHIT 拼接与 contig 检查。",
    "viral_report": "为每个样本独立生成病毒解读报告，并自动建立可点击的批次病毒检出总览。",
    "fine_annotation": "对 CheckV 默认候选或自定义候选 contigs，追加 NR DIAMOND、MEGAN RMA6 和 TaxonKit LCA 注释。",
}
st.markdown(f"<div class='task-note'>{task_copy[task]}</div>", unsafe_allow_html=True)

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

    elif task == "viral_report":
        with st.container(border=True):
            st.subheader("已有流程结果")
            assembly_text = st.text_input("assembly 路径", placeholder="/data/project/assembly")
            clean_text = st.text_input("cleandata 路径", placeholder="/data/project/cleandata")
            clean_layout = layout_labels[st.radio("clean reads 存放方式", list(layout_labels), horizontal=True)]
            read_type = {"双端 PE": "pe", "单端 SE": "se"}[st.radio("测序类型", ["双端 PE", "单端 SE"], horizontal=True)]
        with st.container(border=True):
            st.subheader("报告输出与批次总览")
            output_text = st.text_input("viral_report 输出路径", placeholder="/data/project/viral_report")
            overview_rank = {"病毒科（推荐）": "family", "病毒属": "genus", "病毒种（仅适合高可信注释）": "species"}[st.radio("批次总览默认分类层级", ["病毒科（推荐）", "病毒属", "病毒种（仅适合高可信注释）"], horizontal=True)]
            st.caption("每个样本独立进行候选去冗余、reads 回贴和出报告；批次页面仅汇总“哪些分类单元出现在哪些样本”，不进行跨样本生态统计。")
        c1, c2 = st.columns(2)
        with c1: threads = st.number_input("病毒分析线程数", 1, MAX_VIRAL_THREADS, min(8, MAX_VIRAL_THREADS))
        with c2: viral_min = st.number_input("geNomad 输入最短 contig 长度", 200, 100000, max(200, int(SETTINGS.get("VIRAL_MIN_CONTIG_LEN", "1000"))))
        unified_minimum = st.checkbox(
            "CheckV后报告与vOTU沿用同一长度阈值",
            value=True,
            help="推荐保持勾选。取消后可为CheckV输出单独设置更严格的长度下限。",
        )
        post_checkv_min = st.number_input(
            "CheckV后报告/vOTU最短病毒片段长度",
            200,
            100000,
            int(viral_min) if unified_minimum else max(
                200, int(SETTINGS.get("VOTU_POST_CHECKV_MIN_LEN", "1000"))
            ),
            disabled=unified_minimum,
        )
        if unified_minimum:
            post_checkv_min = viral_min
        st.caption(
            f"本次有效阈值：geNomad输入 ≥ {viral_min} bp；"
            f"CheckV后进入报告/vOTU ≥ {post_checkv_min} bp。"
        )
        with st.expander("查看本次其余有效病毒分析参数"):
            st.dataframe([
                {"参数": "vOTU ANI (%)", "值": SETTINGS.get("VOTU_ANI", "95"), "来源": "pipeline.env"},
                {"参数": "vOTU aligned fraction (%)", "值": SETTINGS.get("VOTU_ALIGNED_FRACTION", "85"), "来源": "pipeline.env"},
                {"参数": "CoverM read identity (%)", "值": SETTINGS.get("COVERM_MIN_READ_PERCENT_IDENTITY", "95"), "来源": "pipeline.env"},
                {"参数": "CoverM read aligned (%)", "值": SETTINGS.get("COVERM_MIN_READ_ALIGNED_PERCENT", "75"), "来源": "pipeline.env"},
                {"参数": "CoverM covered fraction (%)", "值": SETTINGS.get("COVERM_MIN_COVERED_FRACTION", "10"), "来源": "pipeline.env"},
                {"参数": "重要性相对丰度 (%)", "值": SETTINGS.get("VOTU_IMPORTANCE_RELATIVE_ABUNDANCE", "5"), "来源": "pipeline.env"},
                {"参数": "批次默认分类层级", "值": overview_rank, "来源": "当前页面"},
            ], use_container_width=True, hide_index=True)
        enable_vs2 = st.checkbox("启用 VirSorter2 交叉验证（可选，耗时更高）")
        if not SETTINGS.get("GENOMAD_DB") or not SETTINGS.get("CHECKV_DB"):
            st.warning("GENOMAD_DB 或 CHECKV_DB 尚未配置；提交前请先完成 config/pipeline.env 设置。")
        assembly = validate_path(assembly_text, "assembly 路径", exists=True)
        clean = validate_path(clean_text, "cleandata 路径", exists=True)
        output = validate_path(output_text, "viral_report 输出路径", exists=False)
        command = [str(VIRAL_REPORT_SCRIPT), "--assembly-dir", str(assembly), "--cleandata-dir", str(clean), "--clean-layout", clean_layout, "--read-type", read_type, "--output-dir", str(output), "--threads", str(threads), "--min-contig-length", str(viral_min), "--post-checkv-min-length", str(post_checkv_min), "--overview-rank", overview_rank]
        if enable_vs2: command.append("--enable-virsorter2")
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
            threads = st.number_input("DIAMOND 线程数", 1, MAX_VIRAL_THREADS, min(8, MAX_VIRAL_THREADS))
            max_hits = st.number_input("每条 contig 最大保留命中数", 1, 200, int(SETTINGS.get("DIAMOND_NR_MAX_TARGET_SEQS", "25")), help="TaxonKit LCA 需要保留多个命中；不建议设为 1。")
            if mode in {"megan", "both"} and (not SETTINGS.get("MEGAN_DAA2RMA") or not SETTINGS.get("MEGAN_MAP_DB")):
                st.warning("选择 RMA6 时需要在配置中填写 MEGAN_DAA2RMA 与 MEGAN_MAP_DB。")
            if mode in {"taxonomy", "both"} and not SETTINGS.get("TAXONKIT_DB"):
                st.warning("选择 LCA 时需要在配置中填写 TAXONKIT_DB。")
        if scope == "custom" and not taxonlist:
            raise ValueError("自定义 NCBI TaxID 不能为空。")
        command = [str(FINE_ANNOTATION_SCRIPT), "--source", source, "--mode", mode, "--threads", str(threads), "--taxon-scope", scope, "--max-target-seqs", str(max_hits)]
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
        st.code(" \\\n  ".join(shlex.quote(part) for part in command), language="bash")
        c_run, c_refresh = st.columns([1, 1])
        with c_run:
            if st.button("提交后台任务", type="primary", use_container_width=True):
                pid = launch(command)
                st.success(f"任务已提交至 Linux 后台（启动 PID：{pid}）。可刷新或重新打开页面查看状态。")
        with c_refresh:
            st.button("刷新运行中心", use_container_width=True)
except ValueError as error:
    st.info(str(error))

show_status(state_base)
show_adapter_evidence(adapter_evidence_root)
if task == "viral_report":
    show_report_center(state_base)
