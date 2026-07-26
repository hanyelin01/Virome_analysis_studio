#!/usr/bin/env python3
"""Create offline global and per-sample reports for the virome catalogue."""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import re
from collections import Counter
from pathlib import Path


BALTIMORE_LABELS = {
    "I": "I：双链 DNA（dsDNA）病毒",
    "II": "II：单链 DNA（ssDNA）病毒",
    "III": "III：双链 RNA（dsRNA）病毒",
    "IV": "IV：正链单链 RNA［(+)ssRNA］病毒",
    "V": "V：负链单链 RNA［(−)ssRNA］病毒",
    "VI": "VI：单链 RNA 逆转录（ssRNA-RT）病毒",
    "VII": "VII：双链 DNA 逆转录（dsDNA-RT）病毒",
    "UNCLASSIFIED": "未分类：ICTV 参考库未明确映射",
    "": "未分类：ICTV 参考库未明确映射",
}


def rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def esc(value: object) -> str:
    return html.escape(str(value or ""))


def baltimore_label(value: str) -> str:
    value = (value or "").strip().upper()
    return BALTIMORE_LABELS.get(value, f"{value}：未识别的巴尔的摩组别")


def report_filename(sample_id: str) -> str:
    stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", sample_id).strip("._") or "sample"
    digest = hashlib.sha1(sample_id.encode("utf-8")).hexdigest()[:10]
    return f"{stem}_{digest}.html"


def numeric(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def metric_cards(metrics: list[tuple[str, object]]) -> str:
    return "".join(f"<div class='metric'><small>{esc(label)}</small><strong>{esc(value)}</strong></div>" for label, value in metrics)


def table(headers: list[str], body: str, empty: str = "无结果") -> str:
    if not body:
        body = f"<tr><td colspan='{len(headers)}'>{esc(empty)}</td></tr>"
    head = "".join(f"<th>{esc(header)}</th>" for header in headers)
    return f"<div class='tablewrap'><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>"


STYLE = """body{margin:0;background:#f5f8fb;color:#183142;font:14px/1.55 Arial,'Microsoft YaHei',sans-serif}main{max-width:1420px;margin:auto;padding:30px}h1,h2{color:#0b3954}.sub{color:#587182}.metrics{display:grid;grid-template-columns:repeat(5,1fr);gap:12px;margin:22px 0}.metric,.panel{background:#fff;border:1px solid #d8e4ea;border-radius:13px;padding:16px;box-shadow:0 4px 14px #163a5f10}.metric small{color:#587182}.metric strong{display:block;color:#0b6d78;font-size:27px}.grid{display:grid;grid-template-columns:1fr 1fr;gap:16px;margin:16px 0}table{border-collapse:collapse;width:100%}td,th{padding:8px;border-bottom:1px solid #e5edf1;text-align:left;vertical-align:top}th{color:#587182}.bar{display:inline-block;height:11px;background:#2479aa;border-radius:99px;margin-right:8px;vertical-align:middle;min-width:2px}.teal{background:#11877f}.tablewrap{overflow:auto;max-height:530px}.note{border-left:4px solid #d79532;background:#fffaf0;padding:12px 15px;border-radius:7px}a{color:#0b6d78;text-decoration:none;font-weight:600}a:hover{text-decoration:underline}.venn-wrap{display:flex;justify-content:center;overflow:auto}.venn{width:min(100%,1100px);height:auto}.venn-bg{fill:url(#venn-bg);stroke:#d7e4eb}.venn-a{fill:#1d80b7;fill-opacity:.52;stroke:#14608a;stroke-width:2}.venn-b{fill:#8366bf;fill-opacity:.52;stroke:#573f8d;stroke-width:2}.venn-c{fill:#12958d;fill-opacity:.42;stroke:#087069;stroke-width:2}.venn-main-label{font:600 21px Arial,sans-serif;fill:#075c57}.venn-main-value{font:700 27px Arial,sans-serif;fill:#075c57}.venn-leader{stroke:#587182;stroke-width:1.5;fill:none}.venn-key{stroke-width:1.5}.venn-a-key{fill:#e7f3fa;stroke:#14608a}.venn-b-key{fill:#f0ebf8;stroke:#573f8d}.venn-key-name{font:600 15px Arial,sans-serif;fill:#183142}.venn-key-value{font:500 14px Arial,sans-serif;fill:#183142}.venn-exact{font:500 14px Arial,sans-serif;fill:#183142}.evidence-controls{display:flex;gap:10px;align-items:end;flex-wrap:wrap;margin-bottom:10px}.evidence-controls label{display:grid;gap:3px;color:#587182;font-size:12px}.evidence-controls input,.evidence-controls select,.evidence-controls button{padding:7px 9px;border:1px solid #c8d7df;border-radius:7px;background:white;color:#183142}.evidence-controls button{cursor:pointer;background:#e9f4f3;color:#075c57;font-weight:600}#final-evidence-table th[data-sort]{cursor:pointer;white-space:nowrap;user-select:none}#final-evidence-table tr[hidden]{display:none}@media(max-width:800px){.metrics,.grid{grid-template-columns:1fr}.venn-exact{font-size:12px}}"""


def document(title: str, content: str) -> str:
    return f"<!doctype html><html lang='zh-CN'><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>{esc(title)}</title><style>{STYLE}</style><main>{content}</main></html>"


DISCOVERY_TOOLS = ("geNomad", "VirSorter2", "DIAMOND-NR-virus")


def discovery_tools(row: dict[str, str]) -> frozenset[str]:
    columns = {"geNomad": "geNomad", "VirSorter2": "VirSorter2", "DIAMOND-NR-virus": "DIAMOND_NR_virus"}
    if any(column in row for column in columns.values()):
        return frozenset(tool for tool, column in columns.items() if row.get(column, "").lower() == "yes")
    return frozenset(tool for tool in DISCOVERY_TOOLS if tool in row.get("discovery_pattern", ""))


def venn_diagram(discovery: list[dict[str, str]]) -> str:
    regions = Counter(discovery_tools(row) for row in discovery)
    totals = Counter(tool for methods in regions for tool in methods for _ in range(regions[methods]))
    total_vc = len(discovery)
    def count(*tools: str) -> int:
        return regions[frozenset(tools)]
    maximum = max(totals.values(), default=1)
    def ellipse(tool: str) -> tuple[float, float]:
        # All ellipses have the same aspect ratio, so their areas are strictly
        # proportional to the number of VC records detected by the method.
        scale = (totals[tool] / maximum) ** 0.5
        return 465 * scale, 132 * scale
    diamond_rx, diamond_ry = ellipse("DIAMOND-NR-virus")
    virsorter_rx, virsorter_ry = ellipse("VirSorter2")
    genomad_rx, genomad_ry = ellipse("geNomad")
    def legend(tool: str, x: int, y: int, style: str) -> str:
        percentage = totals[tool] / total_vc * 100 if total_vc else 0
        return f"<g transform='translate({x} {y})'><rect width='238' height='54' rx='12' class='venn-key {style}'/><text x='15' y='22' class='venn-key-name'>{esc(tool)}</text><text x='15' y='42' class='venn-key-value'>{totals[tool]:,} VC（{percentage:.2f}%）</text></g>"
    exact = f"DIAMOND 独有 {count('DIAMOND-NR-virus'):,}　·　VirSorter2 独有 {count('VirSorter2'):,}　·　geNomad 独有 {count('geNomad'):,}　·　D∩V {count('DIAMOND-NR-virus', 'VirSorter2'):,}　·　D∩G {count('DIAMOND-NR-virus', 'geNomad'):,}　·　三者共有 {count(*DISCOVERY_TOOLS):,}"
    return f"""<div class='venn-wrap'><svg class='venn' viewBox='0 0 1160 440' role='img' aria-label='按各方法检出 VC 数量比例缩放的 geNomad、VirSorter2 与 DIAMOND-NR-virus 横向韦恩图'><title>按检出量比例缩放的三工具病毒候选韦恩图</title><desc>椭圆面积按各工具检出的完全去冗余 VC 总数严格成比例。由于 DIAMOND-NR-virus 的检出量远大于另外两种方法，较小集合以放大标注辅助阅读；所有交集数量以底部文字精确列出。</desc><defs><filter id='venn-shadow' x='-20%' y='-20%' width='140%' height='140%'><feDropShadow dx='0' dy='7' stdDeviation='7' flood-opacity='.18'/></filter><linearGradient id='venn-bg' x1='0' x2='1'><stop stop-color='#edf7fb'/><stop offset='1' stop-color='#f4f0fb'/></linearGradient></defs><rect x='16' y='16' width='1128' height='408' rx='26' class='venn-bg'/><g filter='url(#venn-shadow)'><ellipse cx='544' cy='214' rx='{diamond_rx:.2f}' ry='{diamond_ry:.2f}' class='venn-c'/><ellipse cx='1003' cy='252' rx='{virsorter_rx:.2f}' ry='{virsorter_ry:.2f}' class='venn-b'/><ellipse cx='976' cy='199' rx='{genomad_rx:.2f}' ry='{genomad_ry:.2f}' class='venn-a'/></g><text x='544' y='203' class='venn-main-label' text-anchor='middle'>DIAMOND-NR-virus</text><text x='544' y='232' class='venn-main-value' text-anchor='middle'>{totals['DIAMOND-NR-virus']:,} VC</text><path d='M 862 93 L 949 191' class='venn-leader'/><path d='M 902 326 L 986 258' class='venn-leader'/>{legend('geNomad', 60, 60, 'venn-a-key')}{legend('VirSorter2', 60, 306, 'venn-b-key')}<text x='580' y='382' class='venn-exact' text-anchor='middle'>{esc(exact)}</text></svg></div>"""


def final_evidence_table(final: list[dict[str, str]]) -> str:
    headers = ["VF ID", "来源 contig（样本标签）", "VC 父记录", "判定", "CheckV", "NR 科", "ICTV 物种候选", "遗传物质类型", "ICTV identity"]
    decisions = sorted({row.get("decision", "") for row in final if row.get("decision", "")})
    samples = sorted({sample.strip() for row in final for sample in row.get("source_sample_ids", "").split(";") if sample.strip()})
    groups = sorted({baltimore_label(row.get("baltimore_group", "")) for row in final})
    options = lambda values: "".join(f"<option value='{esc(value)}'>{esc(value)}</option>" for value in values)
    body = "".join(
        "<tr data-decision='{decision}' data-samples='{samples}' data-baltimore='{baltimore}'>".format(decision=esc(row.get("decision", "")), samples=esc(row.get("source_sample_ids", "")), baltimore=esc(baltimore_label(row.get("baltimore_group", ""))))
        + "".join(f"<td>{esc(value)}</td>" for value in (row.get("vf_id", ""), row.get("source_contig_ids", ""), row.get("parent_vc_id", ""), row.get("decision", ""), row.get("checkv_quality", ""), row.get("nr_family", ""), row.get("ictv_species", ""), baltimore_label(row.get("baltimore_group", "")), row.get("ictv_pident", "")))
        + "</tr>"
        for row in final
    )
    head = "".join(f"<th scope='col' data-sort='{index}'>{esc(header)} ↕</th>" for index, header in enumerate(headers))
    return f"""<div class='evidence-controls'><label>关键词检索<input id='evidence-search' type='search' placeholder='VF、样本、contig、分类…'></label><label>判定<select id='evidence-decision'><option value=''>全部</option>{options(decisions)}</select></label><label>来源样本<select id='evidence-sample'><option value=''>全部</option>{options(samples)}</select></label><label>遗传物质类型<select id='evidence-baltimore'><option value=''>全部</option>{options(groups)}</select></label><button id='evidence-reset' type='button'>重置筛选</button><button id='evidence-download' type='button'>下载当前结果 TSV</button></div><p id='evidence-status' class='sub'></p><div class='tablewrap'><table id='final-evidence-table'><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div><script>(() => {{ const table=document.getElementById('final-evidence-table'), body=table.tBodies[0], search=document.getElementById('evidence-search'), decision=document.getElementById('evidence-decision'), sample=document.getElementById('evidence-sample'), baltimore=document.getElementById('evidence-baltimore'), status=document.getElementById('evidence-status'); let column=-1, ascending=true; const apply=() => {{ let shown=0; for (const row of body.rows) {{ const text=row.innerText.toLowerCase(), sampleValues=row.dataset.samples.split(';').map(value=>value.trim()); const visible=(!search.value || text.includes(search.value.toLowerCase())) && (!decision.value || row.dataset.decision===decision.value) && (!sample.value || sampleValues.includes(sample.value)) && (!baltimore.value || row.dataset.baltimore===baltimore.value); row.hidden=!visible; if (visible) shown++; }} status.textContent=`当前显示 ${{shown}} / ${{body.rows.length}} 条最终片段记录`; }}; [search,decision,sample,baltimore].forEach(control=>control.addEventListener('input',apply)); document.getElementById('evidence-reset').addEventListener('click',()=>{{search.value='';decision.value='';sample.value='';baltimore.value='';apply();}}); for (const header of table.tHead.rows[0].cells) header.addEventListener('click',()=>{{const next=Number(header.dataset.sort); ascending=next===column ? !ascending : true; column=next; [...body.rows].sort((left,right)=>left.cells[column].innerText.localeCompare(right.cells[column].innerText,undefined,{{numeric:true,sensitivity:'base'}})*(ascending?1:-1)).forEach(row=>body.append(row)); apply();}}); document.getElementById('evidence-download').addEventListener('click',()=>{{const lines=[[...table.tHead.rows[0].cells].map(cell=>cell.innerText.replace(' ↕','')).join('\\t')]; for (const row of body.rows) if (!row.hidden) lines.push([...row.cells].map(cell=>cell.innerText.replace(/[\\t\\r\\n]+/g,' ')).join('\\t')); const blob=new Blob([lines.join('\\n')+'\\n'],{{type:'text/tab-separated-values;charset=utf-8'}}), link=document.createElement('a'); link.href=URL.createObjectURL(blob); link.download='final_virome_evidence_filtered.tsv'; link.click(); URL.revokeObjectURL(link.href);}}); apply(); }})();</script>"""


def sample_report(sample_id: str, data: list[dict[str, str]], global_href: str) -> str:
    detected = [row for row in data if row.get("detected", "").lower() == "yes"]
    families = Counter((row.get("ictv_family") or row.get("nr_family") or "未分类") for row in detected)
    baltimore = Counter(baltimore_label(row.get("baltimore_group", "")) for row in detected)
    total_reads = sum(numeric(row.get("read_count", "")) for row in detected)
    family_rows = "".join(f"<tr><td>{esc(name)}</td><td>{count}</td></tr>" for name, count in families.most_common())
    baltimore_rows = "".join(f"<tr><td>{esc(name)}</td><td>{count}</td></tr>" for name, count in baltimore.most_common())
    detail_rows = "".join(
        "<tr>" + "".join(f"<td>{esc(value)}</td>" for value in (
            row.get("vf_id", ""), row.get("source_contig_ids", ""), row.get("checkv_quality", ""), row.get("nr_family", ""),
            row.get("ictv_species", ""), baltimore_label(row.get("baltimore_group", "")),
            row.get("relative_abundance", "0"), row.get("mean_coverage", "0"), row.get("read_count", "0"), row.get("detected", "no"),
        )) + "</tr>"
        for row in data
    )
    return document(
        f"病毒样本报告：{sample_id}",
        f"<p><a href='{esc(global_href)}'>← 返回全局总览</a></p><h1>单样本病毒分析报告</h1><p class='sub'>样本：<strong>{esc(sample_id)}</strong>。定量仅针对最终 CheckV 修正病毒片段（VF）；“检出”表示该 VF 的 CoverM 原始比对读数大于 0。</p>"
        f"<div class='metrics'>{metric_cards([('分发的最终 VF', len(data)), ('检出 VF', len(detected)), ('累计原始读数', f'{total_reads:.0f}'), ('检出病毒科', len(families)), ('检出巴尔的摩类型', len(baltimore))])}</div>"
        f"<div class='grid'><section class='panel'><h2>检出病毒科</h2>{table(['病毒科（ICTV 优先）', '检出 VF 数'], family_rows)}</section><section class='panel'><h2>检出病毒的遗传物质类型</h2>{table(['巴尔的摩分类（可读名称）', '检出 VF 数'], baltimore_rows)}</section></div>"
        f"<section class='panel'><h2>该样本的最终病毒片段及定量</h2>{table(['VF ID', '来源 contig（样本标签）', 'CheckV', 'NR 科', 'ICTV 物种候选', '遗传物质类型', '相对丰度 (%)', '平均覆盖度', '原始读数', '检出'], detail_rows)}<p class='note'>ICTV 物种候选仅是序列相似性证据；应结合覆盖区域、阴阳性对照和实验背景复核，不能单独作为临床诊断。</p></section>",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    root = args.output_dir
    discovery = rows(root / "03_candidate_catalogue" / "VC_discovery_evidence.tsv")
    decisions = rows(root / "04_nr_annotation" / "viral_decision.tsv")
    final = rows(root / "07_final_catalogue" / "VF_catalogue.tsv")
    final_by_vf = {row.get("vf_id", ""): row for row in final if row.get("vf_id", "")}
    presence = rows(root / "08_sample_results" / "sample_fragment_presence.tsv")
    patterns = Counter(row.get("discovery_pattern", "未知") for row in discovery)
    decision_counts = Counter(row.get("decision", "未知") for row in decisions)
    family = Counter((row.get("ictv_family") or row.get("nr_family") or "未分类") for row in final)
    baltimore = Counter(baltimore_label(row.get("baltimore_group", "")) for row in final)
    samples = sorted({row.get("sample_id", "") for row in presence if row.get("sample_id", "")})
    reports = root / "reports"
    sample_reports = reports / "samples"
    sample_reports.mkdir(parents=True, exist_ok=True)

    sample_rows = []
    for sample_id in samples:
        quantified = rows(root / "08_sample_results" / sample_id / "viral_fragments_quantified.tsv")
        if not quantified:
            quantified = rows(root / "08_sample_results" / sample_id / "viral_fragments.tsv")
        quantified = [{**row, **{key: value for key, value in final_by_vf.get(row.get("vf_id", ""), {}).items() if key in ("source_sample_ids", "source_contig_ids")}} for row in quantified]
        filename = report_filename(sample_id)
        (sample_reports / filename).write_text(sample_report(sample_id, quantified, "../virome_catalogue_dashboard.html"), encoding="utf-8")
        detected = sum(row.get("detected", "").lower() == "yes" for row in quantified)
        sample_rows.append(f"<tr><td><a href='samples/{esc(filename)}'>{esc(sample_id)}</a></td><td>{len(quantified)}</td><td>{detected}</td><td><a href='samples/{esc(filename)}'>打开单样本报告</a></td></tr>")

    max_family = max(family.values(), default=1)
    family_rows = "".join(f"<tr><td>{esc(name)}</td><td><span class='bar teal' style='width:{count / max_family * 100:.1f}%'></span>{count}</td></tr>" for name, count in family.most_common(30))
    decision_rows = "".join(f"<tr><td>{esc(name)}</td><td>{count}</td></tr>" for name, count in decision_counts.most_common())
    baltimore_rows = "".join(f"<tr><td>{esc(name)}</td><td>{count}</td></tr>" for name, count in baltimore.most_common())
    metrics = [("全局候选 VC", len(discovery)), ("确认病毒", decision_counts["confirmed_viral"]), ("新颖候选", decision_counts["putative_novel_virus"]), ("CheckV 最终片段 VF", len(final)), ("检出样本", len(samples))]
    report = document(
        "Virome Catalogue v2",
        f"<h1>病毒全局目录分析报告 <span class='sub'>v2</span></h1><p class='sub'>从 cleandata 组装 contigs 经多工具发现、全 NR 分类、CheckV 和 ICTV 本地参考库精细注释生成。VC 是完全去冗余候选序列，VF 是 CheckV 修正后的最终病毒片段；两者均不是 vOTU 或病毒物种。</p><div class='metrics'>{metric_cards(metrics)}</div>"
        f"<section class='panel'><h2>① 三工具潜在病毒序列池：按检出量比例缩放的横向韦恩图</h2><p class='sub'>椭圆面积严格按每种方法检出的完全去冗余 VC 总数缩放。由于 DIAMOND-NR-virus 的检出量远大于其他方法，小集合通过引线标签呈现；底部文字给出全部精确交集计数。</p>{venn_diagram(discovery)}</section>"
        f"<div class='grid'><section class='panel'><h2>② 全 NR 证据判定</h2>{table(['结论', 'VC 数'], decision_rows)}<p class='note'>缺少 NR 命中但有至少两种发现方法支持的序列保留为“新颖候选”，不会因数据库中没有近缘参考而被删除。</p></section><section class='panel'><h2>③ CheckV 后病毒科</h2>{table(['病毒科（ICTV 优先）', 'VF 数'], family_rows)}</section></div>"
        f"<div class='grid'><section class='panel'><h2>④ 巴尔的摩分类：按遗传物质类型解释</h2>{table(['巴尔的摩分类（可读名称）', 'VF 数'], baltimore_rows)}<p class='sub'>I–VII 对应病毒基因组类型与复制策略；标签只来自版本化 ICTV 参考元数据。没有明确映射时保留“未分类”。</p></section><section class='panel'><h2>⑤ 单样本报告</h2><p class='sub'>每个样本都有独立的最终 VF、定量和分类结果页面。</p>{table(['样本', '分发 VF', '检出 VF', '报告'], ''.join(sample_rows), '没有可分发的最终病毒片段')}</section></div>"
        f"<section class='panel'><h2>⑥ 最终病毒片段证据表</h2><p class='sub'>来源 contig 为拼接阶段保留的“样本 ID__原始 contig ID”，可通过表头排序、筛选器或关键词快速定位；下载仅导出当前筛选后的记录。</p>{final_evidence_table(final)}<p class='sub'>“ICTV 物种候选”是序列相似性证据，需同时结合覆盖度、支持区域、阴阳性对照及实验背景复核；不构成临床诊断。</p></section>",
    )
    (reports / "virome_catalogue_dashboard.html").write_text(report, encoding="utf-8")
    print(f"[INFO] Virome catalogue report: {reports / 'virome_catalogue_dashboard.html'}; sample reports: {len(samples)}")


if __name__ == "__main__":
    main()
