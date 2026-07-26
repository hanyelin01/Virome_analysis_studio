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


STYLE = """body{margin:0;background:#f5f8fb;color:#183142;font:14px/1.55 Arial,'Microsoft YaHei',sans-serif}main{max-width:1420px;margin:auto;padding:30px}h1,h2{color:#0b3954}.sub{color:#587182}.metrics{display:grid;grid-template-columns:repeat(5,1fr);gap:12px;margin:22px 0}.metric,.panel{background:#fff;border:1px solid #d8e4ea;border-radius:13px;padding:16px;box-shadow:0 4px 14px #163a5f10}.metric small{color:#587182}.metric strong{display:block;color:#0b6d78;font-size:27px}.grid{display:grid;grid-template-columns:1fr 1fr;gap:16px;margin:16px 0}table{border-collapse:collapse;width:100%}td,th{padding:8px;border-bottom:1px solid #e5edf1;text-align:left;vertical-align:top}th{color:#587182}.bar{display:inline-block;height:11px;background:#2479aa;border-radius:99px;margin-right:8px;vertical-align:middle;min-width:2px}.teal{background:#11877f}.dots b{color:#0b6d78;margin-right:6px}.dots i{color:#ccd8df;margin-right:6px;font-style:normal}.tablewrap{overflow:auto;max-height:530px}.note{border-left:4px solid #d79532;background:#fffaf0;padding:12px 15px;border-radius:7px}.tag{display:inline-block;background:#e9f4f3;color:#075c57;border-radius:99px;padding:2px 8px;font-size:12px}a{color:#0b6d78;text-decoration:none;font-weight:600}a:hover{text-decoration:underline}@media(max-width:800px){.metrics,.grid{grid-template-columns:1fr}}"""


def document(title: str, content: str) -> str:
    return f"<!doctype html><html lang='zh-CN'><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>{esc(title)}</title><style>{STYLE}</style><main>{content}</main></html>"


def sample_report(sample_id: str, data: list[dict[str, str]], global_href: str) -> str:
    detected = [row for row in data if row.get("detected", "").lower() == "yes"]
    families = Counter((row.get("ictv_family") or row.get("nr_family") or "未分类") for row in detected)
    baltimore = Counter(baltimore_label(row.get("baltimore_group", "")) for row in detected)
    total_reads = sum(numeric(row.get("read_count", "")) for row in detected)
    family_rows = "".join(f"<tr><td>{esc(name)}</td><td>{count}</td></tr>" for name, count in families.most_common())
    baltimore_rows = "".join(f"<tr><td>{esc(name)}</td><td>{count}</td></tr>" for name, count in baltimore.most_common())
    detail_rows = "".join(
        "<tr>" + "".join(f"<td>{esc(value)}</td>" for value in (
            row.get("vf_id", ""), row.get("checkv_quality", ""), row.get("nr_family", ""),
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
        f"<section class='panel'><h2>该样本的最终病毒片段及定量</h2>{table(['VF ID', 'CheckV', 'NR 科', 'ICTV 物种候选', '遗传物质类型', '相对丰度 (%)', '平均覆盖度', '原始读数', '检出'], detail_rows)}<p class='note'>ICTV 物种候选仅是序列相似性证据；应结合覆盖区域、阴阳性对照和实验背景复核，不能单独作为临床诊断。</p></section>",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    root = args.output_dir
    discovery = rows(root / "03_candidate_catalogue" / "VC_discovery_evidence.tsv")
    decisions = rows(root / "04_nr_annotation" / "viral_decision.tsv")
    final = rows(root / "07_final_catalogue" / "VF_catalogue.tsv")
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
        filename = report_filename(sample_id)
        (sample_reports / filename).write_text(sample_report(sample_id, quantified, "../virome_catalogue_dashboard.html"), encoding="utf-8")
        detected = sum(row.get("detected", "").lower() == "yes" for row in quantified)
        sample_rows.append(f"<tr><td><a href='samples/{esc(filename)}'>{esc(sample_id)}</a></td><td>{len(quantified)}</td><td>{detected}</td><td><a href='samples/{esc(filename)}'>打开单样本报告</a></td></tr>")

    max_pattern = max(patterns.values(), default=1)
    max_family = max(family.values(), default=1)
    upset = "".join(f"<tr><td>{esc(pattern)}</td><td><span class='dots'>{''.join('<b>●</b>' if tool in pattern else '<i>○</i>' for tool in ('geNomad', 'VirSorter2', 'DIAMOND-NR-virus'))}</span></td><td><span class='bar' style='width:{count / max_pattern * 100:.1f}%'></span>{count}</td></tr>" for pattern, count in patterns.most_common())
    family_rows = "".join(f"<tr><td>{esc(name)}</td><td><span class='bar teal' style='width:{count / max_family * 100:.1f}%'></span>{count}</td></tr>" for name, count in family.most_common(30))
    decision_rows = "".join(f"<tr><td>{esc(name)}</td><td>{count}</td></tr>" for name, count in decision_counts.most_common())
    baltimore_rows = "".join(f"<tr><td>{esc(name)}</td><td>{count}</td></tr>" for name, count in baltimore.most_common())
    detail_rows = "".join("<tr>" + "".join(f"<td>{esc(value)}</td>" for value in (row.get("vf_id", ""), row.get("parent_vc_id", ""), row.get("decision", ""), row.get("checkv_quality", ""), row.get("nr_family", ""), row.get("ictv_species", ""), baltimore_label(row.get("baltimore_group", "")), row.get("ictv_pident", ""))) + "</tr>" for row in final[:500])
    metrics = [("全局候选 VC", len(discovery)), ("确认病毒", decision_counts["confirmed_viral"]), ("新颖候选", decision_counts["putative_novel_virus"]), ("CheckV 最终片段 VF", len(final)), ("检出样本", len(samples))]
    report = document(
        "Virome Catalogue v2",
        f"<h1>病毒全局目录分析报告 <span class='sub'>v2</span></h1><p class='sub'>从 cleandata 组装 contigs 经多工具发现、全 NR 分类、CheckV 和 ICTV 本地参考库精细注释生成。VC 是完全去冗余候选序列，VF 是 CheckV 修正后的最终病毒片段；两者均不是 vOTU 或病毒物种。</p><div class='metrics'>{metric_cards(metrics)}</div>"
        f"<section class='panel'><h2>① 三工具潜在病毒序列池：发现交集</h2><p class='sub'>主图采用 UpSet 风格组合表；实心点依次表示 geNomad、VirSorter2、DIAMOND-NR-virus 的支持。</p>{table(['发现组合', '工具支持', '完全去冗余 VC 数'], upset)}</section>"
        f"<div class='grid'><section class='panel'><h2>② 全 NR 证据判定</h2>{table(['结论', 'VC 数'], decision_rows)}<p class='note'>缺少 NR 命中但有至少两种发现方法支持的序列保留为“新颖候选”，不会因数据库中没有近缘参考而被删除。</p></section><section class='panel'><h2>③ CheckV 后病毒科</h2>{table(['病毒科（ICTV 优先）', 'VF 数'], family_rows)}</section></div>"
        f"<div class='grid'><section class='panel'><h2>④ 巴尔的摩分类：按遗传物质类型解释</h2>{table(['巴尔的摩分类（可读名称）', 'VF 数'], baltimore_rows)}<p class='sub'>I–VII 对应病毒基因组类型与复制策略；标签只来自版本化 ICTV 参考元数据。没有明确映射时保留“未分类”。</p></section><section class='panel'><h2>⑤ 单样本报告</h2><p class='sub'>每个样本都有独立的最终 VF、定量和分类结果页面。</p>{table(['样本', '分发 VF', '检出 VF', '报告'], ''.join(sample_rows), '没有可分发的最终病毒片段')}</section></div>"
        f"<section class='panel'><h2>⑥ 最终病毒片段证据表</h2>{table(['VF ID', 'VC 父记录', '判定', 'CheckV', 'NR 科', 'ICTV 物种候选', '遗传物质类型', 'ICTV identity'], detail_rows)}<p class='sub'>“ICTV 物种候选”是序列相似性证据，需同时结合覆盖度、支持区域、阴阳性对照及实验背景复核；不构成临床诊断。</p></section>",
    )
    (reports / "virome_catalogue_dashboard.html").write_text(report, encoding="utf-8")
    print(f"[INFO] Virome catalogue report: {reports / 'virome_catalogue_dashboard.html'}; sample reports: {len(samples)}")


if __name__ == "__main__":
    main()
