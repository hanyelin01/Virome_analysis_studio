#!/usr/bin/env python3
"""Create a compact offline report for the global viral-catalogue workflow."""
from __future__ import annotations
import argparse
import csv
import html
from collections import Counter, defaultdict
from pathlib import Path


def rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0: return []
    with path.open(encoding="utf-8", errors="replace", newline="") as handle: return list(csv.DictReader(handle, delimiter="\t"))


def esc(value: object) -> str: return html.escape(str(value or ""))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args(); root = args.output_dir
    discovery = rows(root / "03_candidate_catalogue" / "VC_discovery_evidence.tsv")
    decisions = rows(root / "04_nr_annotation" / "viral_decision.tsv")
    final = rows(root / "07_final_catalogue" / "VF_catalogue.tsv")
    presence = rows(root / "08_sample_results" / "sample_fragment_presence.tsv")
    patterns = Counter(row.get("discovery_pattern", "未知") for row in discovery)
    decision_counts = Counter(row.get("decision", "未知") for row in decisions)
    family = Counter((row.get("ictv_family") or row.get("nr_family") or "未分类") for row in final)
    baltimore = Counter(row.get("baltimore_group") or "UNCLASSIFIED" for row in final)
    sample_counts = Counter(row.get("sample_id", "") for row in presence)
    max_pattern = max(patterns.values(), default=1); max_family = max(family.values(), default=1)
    upset = "".join(f"<tr><td>{esc(pattern)}</td><td><span class='dots'>{''.join('<b>●</b>' if tool in pattern else '<i>○</i>' for tool in ('geNomad','VirSorter2','DIAMOND-NR-virus'))}</span></td><td><span class='bar' style='width:{count/max_pattern*100:.1f}%'></span>{count}</td></tr>" for pattern, count in patterns.most_common())
    family_rows = "".join(f"<tr><td>{esc(name)}</td><td><span class='bar teal' style='width:{count/max_family*100:.1f}%'></span>{count}</td></tr>" for name, count in family.most_common(30))
    decision_rows = "".join(f"<tr><td>{esc(name)}</td><td>{count}</td></tr>" for name, count in decision_counts.most_common())
    sample_rows = "".join(f"<tr><td>{esc(name)}</td><td>{count}</td></tr>" for name, count in sorted(sample_counts.items())) or "<tr><td colspan='2'>没有可分发的最终病毒片段</td></tr>"
    detail_rows = "".join("<tr>" + "".join(f"<td>{esc(row.get(field,''))}</td>" for field in ('vf_id','parent_vc_id','decision','checkv_quality','nr_family','ictv_species','baltimore_group','ictv_pident')) + "</tr>" for row in final[:500])
    metrics = [("全局候选 VC", len(discovery)), ("确认病毒", decision_counts['confirmed_viral']), ("新颖候选", decision_counts['putative_novel_virus']), ("CheckV 最终片段 VF", len(final)), ("检出样本", len(sample_counts))]
    metric_html = "".join(f"<div class='metric'><small>{label}</small><strong>{count}</strong></div>" for label, count in metrics)
    report = f"""<!doctype html><html lang='zh-CN'><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>Virome Catalogue v2</title>
<style>body{{margin:0;background:#f5f8fb;color:#183142;font:14px/1.55 Arial,'Microsoft YaHei',sans-serif}}main{{max-width:1420px;margin:auto;padding:30px}}h1,h2{{color:#0b3954}}.sub{{color:#587182}}.metrics{{display:grid;grid-template-columns:repeat(5,1fr);gap:12px;margin:22px 0}}.metric,.panel{{background:#fff;border:1px solid #d8e4ea;border-radius:13px;padding:16px;box-shadow:0 4px 14px #163a5f10}}.metric small{{color:#587182}}.metric strong{{display:block;color:#0b6d78;font-size:30px}}.grid{{display:grid;grid-template-columns:1fr 1fr;gap:16px;margin:16px 0}}table{{border-collapse:collapse;width:100%}}td,th{{padding:8px;border-bottom:1px solid #e5edf1;text-align:left}}th{{color:#587182}}.bar{{display:inline-block;height:11px;background:#2479aa;border-radius:99px;margin-right:8px;vertical-align:middle;min-width:2px}}.teal{{background:#11877f}}.dots b{{color:#0b6d78;margin-right:6px}}.dots i{{color:#ccd8df;margin-right:6px;font-style:normal}}.tablewrap{{overflow:auto;max-height:530px}}.note{{border-left:4px solid #d79532;background:#fffaf0;padding:12px 15px;border-radius:7px}}@media(max-width:800px){{.metrics,.grid{{grid-template-columns:1fr}}}}</style>
<main><h1>病毒全局目录分析报告 <span class='sub'>v2</span></h1><p class='sub'>从 cleandata 组装 contigs 经多工具发现、全 NR 分类、CheckV 和 ICTV 本地参考库精细注释生成。VC 是完全去冗余候选序列，VF 是 CheckV 修正后的最终病毒片段；两者均不是 vOTU 或病毒物种。</p><div class='metrics'>{metric_html}</div>
<section class='panel'><h2>① 三工具潜在病毒序列池：发现交集</h2><p class='sub'>主图采用 UpSet 风格组合表；实心点依次表示 geNomad、VirSorter2、DIAMOND-NR-virus 的支持。点击原始 TSV 可复核每一条调用。</p><div class='tablewrap'><table><thead><tr><th>发现组合</th><th>工具支持</th><th>完全去冗余 VC 数</th></tr></thead><tbody>{upset}</tbody></table></div></section>
<div class='grid'><section class='panel'><h2>② 全 NR 证据判定</h2><table><thead><tr><th>结论</th><th>VC 数</th></tr></thead><tbody>{decision_rows}</tbody></table><p class='note'>缺少 NR 命中但有至少两种发现方法支持的序列保留为“新颖候选”，不会因数据库中没有近缘参考而被删除。</p></section><section class='panel'><h2>③ CheckV 后病毒科</h2><div class='tablewrap'><table><thead><tr><th>病毒科（ICTV 优先）</th><th>VF 数</th></tr></thead><tbody>{family_rows}</tbody></table></div></section></div>
<div class='grid'><section class='panel'><h2>④ 巴尔的摩分类</h2><table><thead><tr><th>组别</th><th>VF 数</th></tr></thead><tbody>{''.join(f'<tr><td>{esc(name)}</td><td>{count}</td></tr>' for name,count in baltimore.most_common())}</tbody></table><p class='sub'>标签只来自版本化 ICTV 参考元数据；没有明确映射时保留 UNCLASSIFIED。</p></section><section class='panel'><h2>⑤ 按原始来源分发到样本</h2><table><thead><tr><th>样本</th><th>最终 VF 数</th></tr></thead><tbody>{sample_rows}</tbody></table></section></div>
<section class='panel'><h2>⑥ 最终病毒片段证据表</h2><div class='tablewrap'><table><thead><tr><th>VF ID</th><th>VC 父记录</th><th>判定</th><th>CheckV</th><th>NR 科</th><th>ICTV 物种候选</th><th>巴尔的摩组</th><th>ICTV identity</th></tr></thead><tbody>{detail_rows}</tbody></table></div><p class='sub'>“ICTV 物种候选”是序列相似性证据，需同时结合覆盖度、支持区域、阴阳性对照及实验背景复核；不构成临床诊断。</p></section></main></html>"""
    reports = root / "reports"; reports.mkdir(exist_ok=True)
    (reports / "virome_catalogue_dashboard.html").write_text(report, encoding="utf-8")
    print(f"[INFO] Virome catalogue report: {reports / 'virome_catalogue_dashboard.html'}")


if __name__ == "__main__": main()
