#!/usr/bin/env python3
"""Build an offline, evidence-aware viral screening dashboard.

The dashboard is deliberately a *navigator*, not a cross-sample ecological
analysis.  Each sample keeps its own local vOTUs and read-mapping evidence;
the batch view only answers which annotated taxa have read-supported local
vOTUs in which samples.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import quote


SCHEMA_VERSION = "2.0"
QUALITY_ORDER = ["Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_sample_groups(path: Path | None, samples: list[str]) -> dict[str, str]:
    if path is None:
        return {}
    rows = read_tsv(path)
    groups: dict[str, str] = {}
    known = set(samples)
    for row in rows:
        sample = (row.get("sample_id") or "").strip()
        group = (row.get("group") or "").strip()
        if not sample or not group:
            raise SystemExit("Groups file requires non-empty sample_id and group columns")
        if sample not in known:
            raise SystemExit(f"Groups file contains a sample outside the manifest: {sample}")
        if sample in groups:
            raise SystemExit(f"Groups file contains a duplicate sample: {sample}")
        groups[sample] = group
    return groups


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def number(value: str | float | int | None) -> float:
    try:
        return float(value or 0)
    except (TypeError, ValueError):
        return 0.0


def integer(value: str | float | int | None) -> int | None:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def safe_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "_", value)


GROUP_LABELS = {
    "DNA": "DNA 病毒",
    "RNA": "RNA 病毒",
    "RT": "逆转录病毒",
    "MIXED": "家族内异质",
    "UNCLASSIFIED": "未归类",
}


def normalise_taxon_name(value: str) -> str:
    """Match taxonomy names conservatively across common annotation formats."""
    return re.sub(r"\s+", " ", (value or "").strip()).casefold()


def source_fingerprint(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {"path": str(path), "available": False}
    return {
        "path": str(path),
        "available": True,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def load_ictv_reference(path: Path) -> dict[str, dict[str, str]]:
    """Load a versioned family-level ICTV Genome mapping, without inference."""
    reference: dict[str, dict[str, str]] = {}
    for row in read_tsv(path):
        family = (row.get("family") or "").strip()
        group = (row.get("genome_group") or "UNCLASSIFIED").strip().upper()
        if family and group in GROUP_LABELS:
            reference[normalise_taxon_name(family)] = {
                "genome_group": group,
                "genome_label": GROUP_LABELS[group],
                "genome_detail": (row.get("genome_detail") or GROUP_LABELS[group]).strip(),
                "dictionary_status": (row.get("dictionary_status") or "reference_mapped").strip(),
                "source_release": (row.get("source_release") or "ICTV MSL").strip(),
                "source_url": (row.get("source_url") or "https://ictv.global/msl").strip(),
            }
    return reference


def load_priority_reference(path: Path, fallback: set[str]) -> dict[str, dict[str, str]]:
    """Load concise, user-maintained attention taxa for navigation only."""
    records: dict[str, dict[str, str]] = {}
    for row in read_tsv(path):
        taxon = (row.get("taxon_name") or "").strip()
        enabled = (row.get("enabled") or "yes").strip().lower()
        if taxon and enabled in {"yes", "true", "1", "on"}:
            records[normalise_taxon_name(taxon)] = {
                "taxon_name": taxon,
                "taxon_rank": (row.get("taxon_rank") or "family").strip(),
                "display_badge": (row.get("display_badge") or "关注科").strip(),
                "review_order": (row.get("review_order") or "9999").strip(),
            }
    if not records:
        for family in fallback:
            records[normalise_taxon_name(family)] = {
                "taxon_name": family,
                "taxon_rank": "family",
                "display_badge": "关注科",
                "review_order": "9999",
            }
    return records


def reference_for_family(family: str, reference: dict[str, dict[str, str]]) -> dict[str, str]:
    return reference.get(normalise_taxon_name(family), {
        "genome_group": "UNCLASSIFIED",
        "genome_label": GROUP_LABELS["UNCLASSIFIED"],
        "genome_detail": "未在本次 ICTV 科级参考字典中匹配；不由 contig 自动推断。",
        "dictionary_status": "not_matched",
        "source_release": "not_matched",
        "source_url": "https://ictv.global/msl",
    })


def normalise_family_set(text: str) -> set[str]:
    return {item.strip().lower() for item in (text or "").split(",") if item.strip()}


def taxonomy_rank(taxonomy: str, rank: str) -> str:
    """Return a taxon name from common geNomad/TaxonKit lineage formats."""
    text = taxonomy or ""
    prefix = {"family": "f", "genus": "g", "species": "s"}[rank]
    match = re.search(rf"(?:^|;)\s*(?:{prefix}__|{rank}[:=_ ]+)\s*([^;]+)", text, re.I)
    if match and match.group(1).strip():
        return match.group(1).strip()
    parts = [part.strip() for part in text.split(";") if part.strip()]
    if rank == "family":
        for part in parts:
            if re.search(r"viridae(?:$|[ _-])", part, re.I):
                return part
        return "未注释到病毒科"
    if rank == "genus":
        for part in reversed(parts):
            if re.search(r"virus$", part, re.I) and "viridae" not in part.lower():
                return part
        return "未注释到病毒属"
    return "未注释到病毒种"


def quality_label(value: str) -> str:
    normalized = (value or "").strip().lower()
    mapping = {
        "complete": "Complete",
        "high-quality": "High-quality",
        "medium-quality": "Medium-quality",
        "low-quality": "Low-quality",
    }
    return mapping.get(normalized, "Not-determined")


def quality_high(row: dict[str, Any]) -> bool:
    return quality_label(str(row.get("checkv_quality", ""))) in {"Complete", "High-quality", "Medium-quality"}


def detected(row: dict[str, Any]) -> bool:
    return str(row.get("detected", "")).strip().lower() == "yes" or number(row.get("relative_abundance")) > 0


def is_priority(row: dict[str, Any]) -> bool:
    value = str(row.get("importance", ""))
    return "高优先" in value or value.lower() in {"high", "high-priority", "priority"}


def read_diamond_annotations(path: Path) -> dict[str, dict[str, str]]:
    """Read optional DIAMOND/TaxonKit annotations keyed by CheckV sequence ID."""
    annotations: dict[str, dict[str, str]] = {}
    for row in read_tsv(path):
        query_id = row.get("query_id", "").split()[0]
        if query_id and row.get("lca_taxid", "") not in {"", "0"}:
            annotations[query_id] = row
    return annotations


def enrich_annotation(row: dict[str, str], diamond: dict[str, dict[str, str]]) -> dict[str, Any]:
    result: dict[str, Any] = dict(row)
    sequence_id = str(result.get("representative_sequence_id", "")).split()[0]
    annotation = diamond.get(sequence_id)
    if annotation:
        result["display_taxonomy"] = annotation.get("lca_lineage") or annotation.get("lca_name") or result.get("taxonomy", "")
        result["display_taxonomy_source"] = "DIAMOND/TaxonKit LCA"
        result["diamond_lca_name"] = annotation.get("lca_name", "")
        result["diamond_lca_rank"] = annotation.get("lca_rank", "")
    else:
        result["display_taxonomy"] = result.get("taxonomy", "")
        result["display_taxonomy_source"] = "geNomad"
        result["diamond_lca_name"] = ""
        result["diamond_lca_rank"] = ""
    return result


def display_taxonomy_rank(row: dict[str, Any], rank: str) -> str:
    if str(row.get("diamond_lca_rank", "")).lower() == rank and row.get("diamond_lca_name", ""):
        return str(row["diamond_lca_name"])
    result = taxonomy_rank(str(row.get("display_taxonomy") or row.get("taxonomy", "")), rank)
    if result != taxonomy_rank("", rank):
        return result
    return taxonomy_rank(str(row.get("taxonomy", "")), rank)


def serialise_row(row: dict[str, Any]) -> dict[str, Any]:
    """Keep a small, stable dashboard schema independent of upstream TSV headers."""
    quality = quality_label(str(row.get("checkv_quality", "")))
    completeness = number(row.get("completeness"))
    contamination = number(row.get("contamination"))
    representative_length = integer(row.get("representative_length"))
    covered_bases = integer(row.get("covered_bases")) or 0
    covered_fraction = None
    if representative_length and representative_length > 0:
        covered_fraction = round(min(100.0, 100.0 * covered_bases / representative_length), 4)
    return {
        "votu_id": str(row.get("votu_id", "")),
        "representative_sequence_id": str(row.get("representative_sequence_id", "")),
        "representative_length": representative_length,
        "member_count": integer(row.get("member_count")) or 0,
        "detected": detected(row),
        "relative_abundance": round(number(row.get("relative_abundance")), 6),
        "mean_coverage": round(number(row.get("mean_coverage")), 6),
        "covered_bases": covered_bases,
        "covered_fraction": covered_fraction,
        "read_count": integer(row.get("read_count")) or 0,
        "checkv_quality": quality,
        "miuvig_quality": str(row.get("miuvig_quality", "Genome-fragment")),
        "completeness": round(completeness, 4) if str(row.get("completeness", "")).strip() not in {"", "NA", "nan"} else None,
        "contamination": round(contamination, 4) if str(row.get("contamination", "")).strip() not in {"", "NA", "nan"} else None,
        "virus_score": str(row.get("virus_score", "NA")),
        "importance": str(row.get("importance", "常规")),
        "priority": is_priority(row),
        "taxonomy": str(row.get("display_taxonomy") or row.get("taxonomy", "Unclassified virus")),
        "taxonomy_source": str(row.get("display_taxonomy_source", "geNomad")),
        "family": display_taxonomy_rank(row, "family"),
        "genus": display_taxonomy_rank(row, "genus"),
        "species": display_taxonomy_rank(row, "species"),
    }


def sample_payload(
    sample: str,
    status: str,
    rows: list[dict[str, Any]],
    priority_taxa: dict[str, dict[str, str]],
    ictv_reference: dict[str, dict[str, str]],
) -> dict[str, Any]:
    for row in rows:
        genome = reference_for_family(str(row["family"]), ictv_reference)
        row.update(genome)
        row["virus_type"] = genome["genome_label"]  # Backward-compatible presentation field.
        priority = priority_taxa.get(normalise_taxon_name(str(row["family"])))
        row["priority_family"] = bool(priority)
        row["priority_badge"] = priority["display_badge"] if priority else ""
    ordered = sorted(rows, key=lambda item: (not item["priority"], not item["detected"], -item["relative_abundance"], item["votu_id"]))
    quality_counts = Counter(row["checkv_quality"] for row in rows)
    summary = {
        "local_votu_count": len(rows),
        "detected_local_votu_count": sum(row["detected"] for row in rows),
        "quality_supported_count": sum(row["detected"] and row["checkv_quality"] in {"Complete", "High-quality", "Medium-quality"} for row in rows),
        "priority_count": sum(row["priority"] for row in rows),
        "priority_family_votu_count": sum(row["detected"] and row["priority_family"] for row in rows),
        "read_count": sum(row["read_count"] for row in rows),
        "relative_abundance": round(sum(row["relative_abundance"] for row in rows if row["detected"]), 6),
    }
    return {
        "sample_id": sample,
        "status": status,
        "summary": summary,
        "quality_counts": {label: quality_counts.get(label, 0) for label in QUALITY_ORDER},
        "rows": ordered,
    }


def build_batch_payload(samples: list[str], payloads: dict[str, dict[str, Any]], rank: str, top_taxa: int, sample_groups: dict[str, str] | None = None) -> dict[str, Any]:
    sample_groups = sample_groups or {}
    presence: dict[str, dict[str, list[dict[str, Any]]]] = defaultdict(lambda: defaultdict(list))
    long_rows: list[dict[str, Any]] = []
    sample_rows: list[dict[str, Any]] = []
    for sample in samples:
        payload = payloads[sample]
        summary = payload["summary"]
        sample_rows.append({"sample_id": sample, "group": sample_groups.get(sample, ""), "status": payload["status"], **summary})
        grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for row in payload["rows"]:
            if row["detected"]:
                grouped[str(row[rank])].append(row)
        for taxon, items in grouped.items():
            presence[taxon][sample] = items
            long_rows.append({
                "rank": rank,
                "taxon": taxon,
                "sample_id": sample,
                "group": sample_groups.get(sample, ""),
                "detected_local_votu_count": len(items),
                "top_relative_abundance": round(max(item["relative_abundance"] for item in items), 6),
                "read_count": sum(item["read_count"] for item in items),
                "candidate_contig_count": sum(item["member_count"] for item in items),
                "high_quality_count": sum(item["checkv_quality"] in {"Complete", "High-quality", "Medium-quality"} for item in items),
                "report_path": f"samples/{safe_id(sample)}.html",
            })
    taxa = []
    for taxon, sample_map in sorted(presence.items(), key=lambda item: (-len(item[1]), item[0])):
        all_items = [item for items in sample_map.values() for item in items]
        groups = Counter(item["genome_group"] for item in all_items)
        group = groups.most_common(1)[0][0] if len(groups) == 1 else "MIXED"
        taxa.append({
            "taxon": taxon,
            "genome_group": group,
            "virus_type": GROUP_LABELS.get(group, GROUP_LABELS["UNCLASSIFIED"]),
            "sample_count": len(sample_map),
            "total_local_votu_count": len(all_items),
            "candidate_contig_count": sum(item["member_count"] for item in all_items),
            "read_count": sum(item["read_count"] for item in all_items),
            "total_relative_abundance": round(sum(item["relative_abundance"] for item in all_items), 6),
            "samples": [{"sample_id": sample, "count": len(items), "top_relative_abundance": max(item["relative_abundance"] for item in items), "read_count": sum(item["read_count"] for item in items)} for sample, items in sorted(sample_map.items())],
        })
    has_read_counts = any(row["read_count"] > 0 for payload in payloads.values() for row in payload["rows"])
    priority_samples = sorted(
        sample_rows,
        key=lambda row: (-row["priority_family_votu_count"], -(row["read_count"] if has_read_counts else row["relative_abundance"]), -row["quality_supported_count"], -row["detected_local_votu_count"], row["sample_id"]),
    )
    return {
        "summary": {
            "sample_count": len(samples),
            "success_count": sum(item["status"] == "SUCCESS" for item in sample_rows),
            "detected_votu_count": sum(item["detected_local_votu_count"] for item in sample_rows),
            "taxon_count": len(taxa),
        },
        "samples": sample_rows,
        "taxa": taxa,
        "matrix_taxa": taxa if top_taxa == 0 else taxa[:top_taxa],
        "priority_samples": priority_samples,
        "presence_rows": long_rows,
        "has_read_counts": has_read_counts,
    }


def dashboard_html(data: dict[str, Any]) -> str:
    """Return one offline HTML document with batch and per-sample navigation."""
    payload = json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
    return f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Virome screening dashboard</title>
<style>
:root{{--ink:#17212b;--muted:#64748b;--paper:#f7f9fc;--surface:#fff;--line:#dde5ee;--navy:#163a5f;--blue:#2374ab;--teal:#1b8a83;--gold:#d39136;--rose:#b7525a;--pale:#eaf1f7;--shadow:0 8px 24px rgba(21,42,64,.08)}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--paper);color:var(--ink);font:14px/1.55 Inter,"Microsoft YaHei",Arial,sans-serif}}button,select{{font:inherit}}.wrap{{max-width:1540px;margin:auto;padding:28px}}.top{{display:flex;align-items:flex-start;justify-content:space-between;gap:20px;margin-bottom:22px}}h1{{margin:0;color:var(--navy);font-size:30px;letter-spacing:-.04em}}h2{{margin:0;font-size:18px;color:var(--navy)}}h3{{margin:0;font-size:14px;color:var(--navy)}}.subtitle,.muted{{color:var(--muted)}}.subtitle{{margin:6px 0 0}}.nav{{display:flex;gap:8px;flex-wrap:wrap}}.nav button,.sample-select{{border:1px solid var(--line);background:#fff;border-radius:8px;padding:8px 12px;color:var(--navy);cursor:pointer}}.nav button.active{{background:var(--navy);color:#fff;border-color:var(--navy)}}.metrics{{display:grid;grid-template-columns:repeat(4,minmax(145px,1fr));gap:14px;margin-bottom:18px}}.metric,.panel,.notice{{background:var(--surface);border:1px solid var(--line);border-radius:12px;box-shadow:var(--shadow)}}.metric{{padding:16px}}.metric .value{{display:block;color:var(--navy);font-size:30px;font-weight:700;line-height:1.1;margin-top:5px}}.metric .label{{color:var(--muted);font-size:12px}}.grid{{display:grid;grid-template-columns:minmax(0,1.35fr) minmax(300px,.65fr);gap:16px;margin:16px 0}}.panel{{padding:18px;min-width:0}}.panel-head{{display:flex;align-items:flex-start;justify-content:space-between;gap:12px;margin-bottom:12px}}details{{margin-top:12px;border-top:1px solid var(--line);padding-top:10px;color:var(--muted)}}summary{{cursor:pointer;color:var(--blue);font-weight:650}}.chart{{min-height:270px}}.bar-row{{display:grid;grid-template-columns:minmax(112px,1fr) minmax(120px,3fr) 42px;gap:10px;align-items:center;margin:10px 0}}.bar-name{{overflow:hidden;text-overflow:ellipsis;white-space:nowrap}}.bar-track{{height:13px;background:var(--pale);border-radius:99px;overflow:hidden}}.bar-fill{{height:100%;border-radius:99px;background:linear-gradient(90deg,var(--blue),var(--teal))}}.bar-value{{font-variant-numeric:tabular-nums;text-align:right;color:var(--muted)}}.rank-list{{list-style:none;margin:0;padding:0}}.rank-list button{{display:grid;grid-template-columns:auto 1fr auto;gap:9px;align-items:center;width:100%;border:0;border-top:1px solid var(--line);background:transparent;padding:11px 0;text-align:left;cursor:pointer;color:var(--ink)}}.rank-list button:first-child{{border-top:0}}.rank-number{{color:var(--gold);font-weight:750}}.pill{{display:inline-block;border-radius:99px;padding:2px 8px;font-size:11px;background:#e8f4f2;color:#087269}}.matrix-wrap,.table-wrap{{overflow:auto;max-width:100%}}table{{border-collapse:collapse;width:100%;font-size:12px}}th,td{{border-bottom:1px solid var(--line);padding:8px;text-align:right;white-space:nowrap}}th:first-child,td:first-child{{text-align:left}}th{{color:var(--muted);font-weight:650;background:#fbfcfe;position:sticky;top:0}}.matrix td button{{width:15px;height:15px;border:0;border-radius:4px;background:#dce5ed;cursor:pointer}}.matrix td button.present{{background:var(--teal)}}.matrix td button:hover{{outline:2px solid var(--gold)}}.notice{{padding:15px 17px;background:#fffaf0;border-left:4px solid var(--gold);box-shadow:none}}.notice strong{{color:#815513}}.page{{display:none}}.page.active{{display:block}}.sample-title{{display:flex;align-items:center;gap:12px;flex-wrap:wrap;margin-bottom:18px}}.quality-grid{{display:grid;grid-template-columns:repeat(5,1fr);gap:8px;margin-top:12px}}.quality-cell{{border-radius:9px;padding:10px;background:var(--pale);min-height:74px}}.quality-cell .qnum{{display:block;font-size:23px;font-weight:720;color:var(--navy)}}.quality-cell .qlabel{{font-size:11px;color:var(--muted)}}.scatter{{position:relative;height:285px;border-left:1px solid #9cacbd;border-bottom:1px solid #9cacbd;margin:20px 14px 30px 38px;background:linear-gradient(to top,transparent 24%,#edf2f6 25%,transparent 26%,transparent 49%,#edf2f6 50%,transparent 51%,transparent 74%,#edf2f6 75%,transparent 76%)}}.dot{{position:absolute;width:10px;height:10px;border-radius:50%;background:var(--teal);transform:translate(-50%,50%);border:1px solid white;cursor:default}}.dot.priority{{background:var(--gold);width:13px;height:13px}}.axis-x,.axis-y{{position:absolute;color:var(--muted);font-size:11px}}.axis-x{{bottom:-25px;left:40%;transform:translateX(-20%)}}.axis-y{{left:-35px;top:45%;transform:rotate(-90deg)}}.source{{font-size:11px;color:var(--muted)}}.source b{{color:var(--ink)}}.link{{color:var(--blue);cursor:pointer;text-decoration:underline}}@media(max-width:930px){{.grid{{grid-template-columns:1fr}}.metrics{{grid-template-columns:repeat(2,1fr)}}.top{{display:block}}.nav{{margin-top:14px}}}}@media(max-width:520px){{.wrap{{padding:17px}}.metrics{{grid-template-columns:1fr 1fr}}h1{{font-size:24px}}}}
</style></head><body><main class="wrap">
<header class="top"><div><h1>病毒筛查报告中心</h1><p class="subtitle">逐样本病毒候选解释与批次检出导航 · 报告生成于 <span id="created"></span></p></div><nav class="nav"><button id="overview-btn" class="active">批次总览</button><button id="sample-btn">单样本判读</button><select id="sample-select" class="sample-select"></select></nav></header>
<section id="overview" class="page active"></section><section id="sample" class="page"></section>
</main><script>const DATA={payload};
const $=s=>document.querySelector(s), esc=v=>String(v??'').replace(/[&<>"']/g,c=>({{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}}[c]));
const n=v=>new Intl.NumberFormat('zh-CN').format(v||0); const pct=v=>`${{Number(v||0).toFixed(2)}}%`;
const sampleMap=Object.fromEntries(DATA.samples.map(x=>[x.sample_id,x]));let current=location.hash.startsWith('#sample=')?decodeURIComponent(location.hash.slice(8)):DATA.samples[0]?.sample_id;
function explanation(title,body){{return `<details><summary>如何解读与证据边界</summary><p><b>${{title}}</b>：${{body}}</p></details>`}}
function cards(items){{return `<div class="metrics">${{items.map(x=>`<article class="metric"><span class="label">${{esc(x[0])}}</span><strong class="value">${{esc(x[1])}}</strong><span class="muted">${{esc(x[2]||'')}}</span></article>`).join('')}}</div>`}}
function renderOverview(){{const b=DATA.batch, top=b.taxa.slice(0,DATA.top_taxa), max=Math.max(1,...top.map(x=>x.sample_count));
 $('#overview').innerHTML=cards([['纳入样本',n(b.summary.sample_count),'本次报告覆盖的样本'],['已完成样本',n(b.summary.success_count),'生信流程状态为 SUCCESS'],['reads 支持本地 vOTU',n(b.summary.detected_votu_count),'逐样本独立计算'],[DATA.rank_label+'检出单元',n(b.summary.taxon_count),'用于批次导航']])+`
 <div class="grid"><article class="panel"><div class="panel-head"><div><h2>主要${{DATA.rank_label}}</h2><p class="muted">按“检出该类群的样本数”排序，而非跨样本丰度。</p></div><span class="pill">Top ${{DATA.top_taxa}}</span></div><div class="chart">${{top.map(x=>`<div class="bar-row"><span class="bar-name" title="${{esc(x.taxon)}}">${{esc(x.taxon)}}</span><span class="bar-track"><span class="bar-fill" style="width:${{x.sample_count/max*100}}%"></span></span><span class="bar-value">${{x.sample_count}}</span></div>`).join('')||'<p class="muted">暂无 reads 支持且可归入所选分类层级的记录。</p>'}}</div>${{explanation('回答的问题','此图显示哪些分类单元在更多样本中出现。一个出现标记表示该样本至少有一个 reads 支持的本地 vOTU；它不比较样本间的病毒丰度，也不代表跨样本的同一病毒株。')}}</article>
 <article class="panel"><div class="panel-head"><div><h2>优先查看样本</h2><p class="muted">按高优先级候选、质量支持和检出数导航。</p></div></div><ol class="rank-list">${{b.priority_samples.map((x,i)=>`<li><button data-sample="${{esc(x.sample_id)}}"><span class="rank-number">${{i+1}}</span><span>${{esc(x.sample_id)}}<br><small class="muted">高优先级 ${{x.priority_count}} · reads 支持 ${{x.detected_local_votu_count}}</small></span><span class="pill">查看</span></button></li>`).join('')||'<li class="muted">暂无可排序的样本结果。</li>'}}</ol>${{explanation('排序逻辑','该列表是浏览顺序建议：优先显示同时具有高质量候选、reads 支持或较高相对丰度的样本。它不是临床风险排序，也不对样本间生物学效应作出结论。')}}</article></div>
 <article class="panel"><div class="panel-head"><div><h2>${{DATA.rank_label}} × 样本检测矩阵</h2><p class="muted">点击绿色方格进入对应样本的本地证据页。</p></div></div><div class="matrix-wrap"><table class="matrix"><thead><tr><th>${{DATA.rank_label}}</th>${{DATA.samples.map(s=>`<th>${{esc(s.sample_id)}}</th>`).join('')}}</tr></thead><tbody>${{b.matrix_taxa.map(t=>`<tr><th title="${{esc(t.taxon)}}">${{esc(t.taxon)}}</th>${{DATA.samples.map(s=>{{const hit=t.samples.find(x=>x.sample_id===s.sample_id);return `<td>${{hit?`<button class="present" data-sample="${{esc(s.sample_id)}}" title="${{esc(s.sample_id)}}：${{hit.count}} 个本地 vOTU；最高相对丰度 ${{pct(hit.top_relative_abundance)}}"></button>`:'<button disabled title="未检出"></button>'}}</td>`}}).join('')}}</tr>`).join('')}}</tbody></table></div>${{explanation('检测矩阵','绿色方格表示对应样本有 reads 支持的本地 vOTU 被注释到这一分类层级。未显示不等于样本不含任何病毒，而是指没有满足当前筛查、质量、回贴和分类条件的记录。')}}</article>
 <aside class="notice"><strong>科研证据边界</strong><br>“检出”表示 reads 支持的病毒候选或本地 vOTU；分类名称来源于 DIAMOND/TaxonKit LCA 或 geNomad 的序列相似性/模型推断。该报告不单独证明感染、宿主范围、致病性、病毒活性或跨样本同一病毒株。</aside>`;
 document.querySelectorAll('[data-sample]').forEach(x=>x.addEventListener('click',()=>openSample(x.dataset.sample)));}}
function qualityCells(q){{return `<div class="quality-grid">${{Object.entries(q).map(([label,count])=>`<div class="quality-cell"><span class="qnum">${{n(count)}}</span><span class="qlabel">${{esc(label)}}</span></div>`).join('')}}</div>`}}
function renderSample(id){{const s=sampleMap[id];if(!s) return; current=id;location.hash='sample='+encodeURIComponent(id);$('#sample-select').value=id;const rows=s.rows, tops=rows.filter(x=>x.detected).slice(0,12), max=Math.max(1,...tops.map(x=>x.relative_abundance));
 $('#sample').innerHTML=`<div class="sample-title"><div><h2>样本 ${{esc(s.sample_id)}}</h2><p class="muted">状态：${{esc(s.status)}}；每一个 vOTU 均仅在本样本内聚类、回贴和判读。</p></div><span class="pill">本地证据页</span></div>${{cards([['本地 vOTU',n(s.summary.local_votu_count),'CheckV 后样本内去冗余'],['reads 支持',n(s.summary.detected_local_votu_count),'有相对丰度记录'],['中/高质量或完整',n(s.summary.quality_supported_count),'CheckV 质量支持'],['高优先级候选',n(s.summary.priority_count),'建议优先复核']])}}
 <div class="grid"><article class="panel"><div class="panel-head"><div><h2>Top 本地 vOTU 相对丰度</h2><p class="muted">仅在本样本内展示 reads 回贴结果。</p></div></div><div class="chart">${{tops.map(x=>`<div class="bar-row"><span class="bar-name" title="${{esc(x.votu_id)}}">${{esc(x.votu_id)}}</span><span class="bar-track"><span class="bar-fill" style="width:${{x.relative_abundance/max*100}}%"></span></span><span class="bar-value">${{pct(x.relative_abundance)}}</span></div>`).join('')||'<p class="muted">无 reads 支持的本地 vOTU。</p>'}}</div>${{explanation('相对丰度','数值来自本样本 clean reads 对本地 vOTU representative 的回贴汇总。它适合在同一样本内确定优先查看对象，不应用于不同样本之间的丰度比较。')}}</article>
 <article class="panel"><div class="panel-head"><div><h2>CheckV 质量概览</h2><p class="muted">候选序列的质量标签，而非独立实验验证。</p></div></div>${{qualityCells(s.quality_counts)}}${{explanation('质量标签','Complete、High-quality 和 Medium-quality 为 CheckV 的序列质量评估结果。质量较低或未确定的候选仍保留在证据表中，避免因展示规则而被静默忽略。')}}</article></div>
 <div class="grid"><article class="panel"><div class="panel-head"><div><h2>完整度与覆盖证据</h2><p class="muted">金色点表示高优先级候选；悬停可读取代表序列信息。</p></div></div><div class="scatter"><span class="axis-y">平均覆盖度</span><span class="axis-x">CheckV 完整度 (%)</span>${{rows.filter(x=>x.completeness!==null).map(x=>{{const left=Math.max(1,Math.min(99,x.completeness));const bottom=Math.max(1,Math.min(99,Math.log10(x.mean_coverage+1)/Math.log10(Math.max(2,...rows.map(y=>y.mean_coverage))+1)*100));return `<span class="dot ${{x.priority?'priority':''}}" style="left:${{left}}%;bottom:${{bottom}}%" title="${{esc(x.votu_id)}} | completeness: ${{x.completeness}}% | mean coverage: ${{x.mean_coverage.toFixed(3)}} | ${{esc(x.checkv_quality)}}"></span>`}}).join('')}}</div>${{explanation('证据散点','横轴是 CheckV 估计完整度，纵轴是回贴的平均覆盖度（对数视觉缩放）。两者用于辅助挑选复核对象，不能单独证明病毒的生物学活性或感染状态。')}}</article>
 <article class="panel"><div class="panel-head"><div><h2>优先候选的分类证据</h2><p class="muted">优先显示分类来源，不将不同工具的结果强制合并。</p></div></div><div class="table-wrap"><table><thead><tr><th>本地 vOTU</th><th>分类名称</th><th>来源</th><th>质量</th></tr></thead><tbody>${{rows.filter(x=>x.priority||x.detected).slice(0,18).map(x=>`<tr><td>${{esc(x.votu_id)}}</td><td title="${{esc(x.taxonomy)}}">${{esc(x.family)}}</td><td><span class="pill">${{esc(x.taxonomy_source)}}</span></td><td>${{esc(x.checkv_quality)}}</td></tr>`).join('')||'<tr><td colspan="4" class="muted">暂无可展示记录。</td></tr>'}}</tbody></table></div>${{explanation('分类来源','若代表序列具有有效 DIAMOND/TaxonKit LCA，报告优先展示该来源；否则显示 geNomad 注释。分类名称是注释线索，而不是最终物种鉴定。')}}</article></div>
 <article class="panel"><div class="panel-head"><div><h2>完整本地 vOTU 证据表</h2><p class="muted">包含未检出、未注释和未确定质量的记录，便于审计与下载。</p></div><a href="data/samples/${{encodeURIComponent(s.sample_id)}}.json" download>下载 JSON</a></div><div class="table-wrap"><table><thead><tr><th>本地 vOTU</th><th>reads 支持</th><th>相对丰度</th><th>平均覆盖度</th><th>CheckV</th><th>完整度</th><th>分类</th><th>来源</th><th>代表序列</th></tr></thead><tbody>${{rows.map(x=>`<tr><td>${{esc(x.votu_id)}}</td><td>${{x.detected?'是':'否'}}</td><td>${{pct(x.relative_abundance)}}</td><td>${{x.mean_coverage.toFixed(3)}}</td><td>${{esc(x.checkv_quality)}}</td><td>${{x.completeness===null?'NA':x.completeness+'%'}}</td><td title="${{esc(x.taxonomy)}}">${{esc(x.family)}}</td><td>${{esc(x.taxonomy_source)}}</td><td>${{esc(x.representative_sequence_id)}}</td></tr>`).join('')}}</tbody></table></div>${{explanation('审计用途','此表应与 reports/data/samples 中的机器可读 JSON 及原始 votu_summary.tsv 一同保存。报告的视觉筛选不会删除原始候选记录。')}}</article>`;}}
function openSample(id){{renderSample(id);$('#overview').classList.remove('active');$('#sample').classList.add('active');$('#overview-btn').classList.remove('active');$('#sample-btn').classList.add('active');}}
function openOverview(){{$('#sample').classList.remove('active');$('#overview').classList.add('active');$('#sample-btn').classList.remove('active');$('#overview-btn').classList.add('active');history.replaceState(null,'','#overview')}}
// Dashboard v2 enhancements: all taxa, attention-family navigation and numeric axes.
function supportLabelV2(){{return DATA.batch.has_read_counts ? 'mapped reads' : '相对丰度合计 (%)';}}
function supportValueV2(item){{return Number(DATA.batch.has_read_counts ? (item.read_count || 0) : (item.total_relative_abundance || 0));}}
function bubbleSizeV2(value, maximum){{return Math.max(7, Math.min(28, 7 + 21 * Math.sqrt((value || 0) / Math.max(1, maximum))));}}
function explanationV2(title, body){{return '<details><summary>如何解读与证据边界</summary><p><b>'+esc(title)+'</b>：'+esc(body)+'</p></details>';}}
function priorityFamiliesV2(sample){{return [...new Set(sample.rows.filter(function(x){{return x.detected && x.priority_family;}}).map(function(x){{return x.family;}}))];}}
var renderOverviewBase=renderOverview;
renderOverview=function(){{renderOverviewBase(); var batch=DATA.batch, taxa=batch.taxa, panels=document.querySelectorAll('#overview .panel'), maxContigs=Math.max(1,...taxa.map(function(x){{return x.candidate_contig_count||0;}})), maxSupport=Math.max(1,...taxa.map(supportValueV2));
 var familyRows=''; taxa.forEach(function(x){{var support=supportValueV2(x), size=bubbleSizeV2(support,maxSupport); familyRows+='<div class="bar-row" style="grid-template-columns:minmax(150px,1.2fr) minmax(150px,3fr) 80px 64px"><span class="bar-name" title="'+esc(x.taxon)+'">'+esc(x.taxon)+'</span><span class="bar-track" title="候选 contig 数：'+n(x.candidate_contig_count)+'"><span class="bar-fill" style="width:'+(x.candidate_contig_count/maxContigs*100)+'%"></span></span><span title="'+supportLabelV2()+'：'+n(support)+'" style="display:inline-flex;align-items:center;gap:7px;white-space:nowrap"><i style="display:inline-block;width:'+size+'px;height:'+size+'px;border-radius:50%;background:var(--gold);opacity:.82"></i>'+n(support)+'</span><span class="bar-value">'+n(x.candidate_contig_count)+'</span></div>';}}); if(!familyRows) familyRows='<p class="muted">暂无 reads 支持且可归入所选分类层级的记录。</p>';
 panels[0].innerHTML='<div class="panel-head"><div><h2>主要'+esc(DATA.rank_label)+'</h2><p class="muted">横条为候选 contig 数；气泡为 '+supportLabelV2()+'。全部分类单元可在面板内滚动查看。</p></div><span class="pill">全部 '+n(taxa.length)+'</span></div><div class="chart" style="max-height:620px;overflow:auto;padding-right:8px">'+familyRows+'</div>'+explanationV2('构图与边界','横条统计所有样本内本地 vOTU 所包含的候选 contig 成员数；气泡显示真实 mapped reads。若报告来自旧版丰度表且没有 count 列，气泡会明确改为相对丰度合计。不同样本仅作导航汇总，不等同于同一病毒株。');
 var priorityRows=''; batch.priority_samples.forEach(function(x,index){{var sample=sampleMap[x.sample_id], families=priorityFamiliesV2(sample), support=DATA.batch.has_read_counts?x.read_count:x.relative_abundance; priorityRows+='<li><button data-sample="'+esc(x.sample_id)+'"><span class="rank-number">'+(index+1)+'</span><span>'+esc(x.sample_id)+'<br><small class="muted">关注科 vOTU '+n(x.priority_family_votu_count)+' · '+supportLabelV2()+' '+n(support)+' · 质量支持 '+n(x.quality_supported_count)+'</small>'+(families.length?'<br><small class="source">关注科：'+esc(families.join('、'))+'</small>':'')+'</span><span class="pill">查看</span></button></li>';}}); if(!priorityRows) priorityRows='<li class="muted">暂无可排序的样本结果。</li>';
 panels[1].innerHTML='<div class="panel-head"><div><h2>优先查看样本</h2><p class="muted">先按关注病毒科的 reads 支持，再按质量支持与检出数排序；全部样本可滚动查看。</p></div><span class="pill">关注科优先</span></div><ol class="rank-list" style="max-height:620px;overflow:auto">'+priorityRows+'</ol>'+explanationV2('优先复核规则','关注病毒科来自管理员配置，默认包括 Coronaviridae、Paramyxoviridae、Orthomyxoviridae 和 Flaviviridae，仅用于优先打开样本报告。它不代表样本存在人畜共患、感染性、致病性或公共卫生风险。');
 var matrixRows=''; batch.matrix_taxa.forEach(function(t){{matrixRows+='<tr><td><span class="pill">'+esc(t.virus_type)+'</span></td><th title="'+esc(t.taxon)+'">'+esc(t.taxon)+'</th>'; DATA.samples.forEach(function(s){{var hit=t.samples.find(function(x){{return x.sample_id===s.sample_id;}}); matrixRows+='<td>'+(hit?'<button class="present" data-sample="'+esc(s.sample_id)+'" title="'+esc(s.sample_id)+'：'+hit.count+' 个本地 vOTU；'+supportLabelV2()+' '+n(hit.read_count||hit.top_relative_abundance)+'"></button>':'<button disabled title="未检出"></button>')+'</td>';}}); matrixRows+='</tr>';}}); 
 panels[2].innerHTML='<div class="panel-head"><div><h2>'+esc(DATA.rank_label)+' × 样本检测矩阵</h2><p class="muted">第一列为科级参考映射的基因组类型；点击绿色方格进入对应样本的本地证据页。全部分类单元可滚动查看。</p></div><span class="pill">全部 '+n(batch.matrix_taxa.length)+'</span></div><div class="matrix-wrap" style="max-height:720px"><table class="matrix"><thead><tr><th>基因组类型</th><th>'+esc(DATA.rank_label)+'</th>'+DATA.samples.map(function(s){{return '<th>'+esc(s.sample_id)+'</th>';}}).join('')+'</tr></thead><tbody>'+matrixRows+'</tbody></table></div>'+explanationV2('基因组类型与矩阵','DNA、RNA、逆转录病毒标签来自内置且可审计的病毒科级显示映射；未在映射表中的分类单元显示为未归类，不会由 contig 自动猜测。绿色方格仍只表示 reads 支持的本地 vOTU。');
 document.querySelectorAll('[data-sample]').forEach(function(x){{x.addEventListener('click',function(){{openSample(x.dataset.sample);}});}});}};
var renderSampleBase=renderSample;
renderSample=function(id){{renderSampleBase(id); var sample=sampleMap[id], grids=document.querySelectorAll('#sample .grid'), topPanel=grids[0].children[0], scatterPanel=grids[1].children[0], rows=sample.rows, topRows=rows.filter(function(x){{return x.detected;}}).slice(0,15), maxAbundance=Math.max(1,...topRows.map(function(x){{return x.relative_abundance;}})), maxSupport=Math.max(1,...topRows.map(function(x){{return DATA.batch.has_read_counts?x.read_count:x.mean_coverage;}}));
 var lollipops=''; topRows.forEach(function(x){{var support=DATA.batch.has_read_counts?x.read_count:x.mean_coverage, size=bubbleSizeV2(support,maxSupport), position=x.relative_abundance/maxAbundance*100, tooltip=x.votu_id+' | 病毒科：'+x.family+' | 属：'+x.genus+' | 分类来源：'+x.taxonomy_source+' | CheckV：'+x.checkv_quality+' | mapped reads：'+n(x.read_count); lollipops+='<div title="'+esc(tooltip)+'" style="display:grid;grid-template-columns:minmax(150px,1.2fr) minmax(150px,3fr) 62px;gap:10px;align-items:center;margin:14px 0;cursor:help"><span class="bar-name">'+esc(x.votu_id)+'<br><small class="source">'+esc(x.family)+' · '+esc(x.checkv_quality)+'</small></span><span class="bar-track" style="height:3px;overflow:visible;position:relative"><span class="bar-fill" style="height:3px;width:'+position+'%;background:var(--blue)"></span><i style="position:absolute;left:calc('+position+'% - '+(size/2)+'px);top:-'+((size-3)/2)+'px;display:block;width:'+size+'px;height:'+size+'px;border-radius:50%;background:'+(x.priority?'var(--gold)':'var(--teal)')+';border:2px solid #fff;box-shadow:0 1px 4px rgba(0,0,0,.2)"></i></span><span class="bar-value">'+pct(x.relative_abundance)+'</span></div>';}}); if(!lollipops) lollipops='<p class="muted">无 reads 支持的本地 vOTU。</p>';
 topPanel.innerHTML='<div class="panel-head"><div><h2>Top 本地 vOTU：丰度与测序支持</h2><p class="muted">横线末端位置为相对丰度；气泡大小为 '+(DATA.batch.has_read_counts?'mapped reads':'平均覆盖度（旧结果回退）')+'。悬停可查看分类信息。</p></div></div><div class="chart" style="max-height:610px;overflow:auto">'+lollipops+'</div>'+explanationV2('点线图编码','相对丰度以线段终点的横向位置精确排序；圆点面积表示 mapped reads。鼠标悬停于任一 vOTU 行可读取病毒科、属、分类来源、CheckV 质量和测序支持。金色圆点表示高优先级候选。');
 var scatterRows=rows.filter(function(x){{return x.completeness!==null;}}), maxCoverage=Math.max(1,...scatterRows.map(function(x){{return x.mean_coverage;}})), logMaximum=Math.max(1,Math.ceil(Math.log10(maxCoverage+1))), xTicks=[0,25,50,75,100].map(function(v){{return '<span style="position:absolute;bottom:-21px;left:'+v+'%;transform:translateX(-50%);font-size:10px;color:var(--muted)">'+v+'</span>';}}).join(''), yTicks=Array.from({{length:logMaximum+1}},function(_,i){{return '<span style="position:absolute;left:-38px;bottom:'+(i/logMaximum*100)+'%;transform:translateY(50%);font-size:10px;color:var(--muted)">'+(Math.pow(10,i)-1).toFixed(i<2?0:1)+'</span>';}}).join(''), dots=''; scatterRows.forEach(function(x){{var bottom=Math.log10(x.mean_coverage+1)/logMaximum*100, tooltip=x.votu_id+' | 病毒科：'+x.family+' | completeness: '+x.completeness+'% | mean coverage: '+x.mean_coverage.toFixed(3)+' | mapped reads: '+n(x.read_count)+' | '+x.checkv_quality; dots+='<span class="dot '+(x.priority?'priority':'')+'" style="left:'+Math.max(1,Math.min(99,x.completeness))+'%;bottom:'+Math.max(1,Math.min(99,bottom))+'%" title="'+esc(tooltip)+'"></span>';}}); 
 scatterPanel.innerHTML='<div class="panel-head"><div><h2>完整度与覆盖证据</h2><p class="muted">横轴为 CheckV 完整度（%）；纵轴为平均覆盖度的 log10(x+1) 变换，并显示实际刻度。</p></div></div><div class="scatter"><span class="axis-y">平均覆盖度（log10(x+1)）</span><span class="axis-x">CheckV 完整度 (%)</span>'+xTicks+yTicks+dots+'</div>'+explanationV2('坐标与证据','完整度坐标保留 0–100% 的实际百分比。平均覆盖度跨度通常很大，因此采用 log10(x+1) 视觉变换，并在纵轴显示对应实际值；悬停点可读取原始平均覆盖度和 mapped reads。该图用于复核优先级，不证明病毒活性或感染状态。');
}};
// Dashboard v3: versioned ICTV groups, compact attention labels and natural evidence axes.
function genomeColourV3(group){{return {{DNA:'#2374ab',RNA:'#1b8a83',RT:'#7c5ab5',MIXED:'#b7525a',UNCLASSIFIED:'#94a3b8'}}[group]||'#94a3b8';}}
function qualityColourV3(quality){{return {{'Complete':'#163a5f','High-quality':'#1b8a83','Medium-quality':'#d39136','Low-quality':'#b7525a','Not-determined':'#94a3b8'}}[quality]||'#94a3b8';}}
function genomeLabelV3(item){{return item.virus_type||'未归类';}}
function attentionCountV3(sample){{return sample.rows.filter(function(x){{return x.detected&&x.priority_family;}}).length;}}
var renderOverviewV3Base=renderOverview;
renderOverview=function(){{
 renderOverviewV3Base();
 var batch=DATA.batch, panels=document.querySelectorAll('#overview .panel');
 var taxa=batch.taxa, maxContigs=Math.max(1,...taxa.map(function(x){{return x.candidate_contig_count||0;}})), maxSupport=Math.max(1,...taxa.map(supportValueV2)), evidenceRows='';
 taxa.forEach(function(x){{
   var contigs=x.candidate_contig_count||0, width=Math.max(contigs?1:0,contigs/maxContigs*100), support=supportValueV2(x), size=Math.max(15,Math.min(42,15+27*Math.sqrt(support/maxSupport))), bubbleLeft=Math.max(0,Math.min(100,width)), bubbleTitle=DATA.batch.has_read_counts?'mapped reads：'+n(support):'相对丰度合计：'+pct(support);
   evidenceRows+='<div style="display:grid;grid-template-columns:minmax(150px,1.15fr) minmax(170px,3.3fr) 76px;gap:12px;align-items:center;margin:17px 0;padding:3px 3px 3px 0"><span class="bar-name" title="'+esc(x.taxon)+'" style="font-weight:680">'+esc(x.taxon)+'<br><small class="source">'+esc(x.virus_type||'未归类')+'</small></span><span class="bar-track" style="height:12px;overflow:visible;position:relative;background:linear-gradient(90deg,#edf2f7,rgba(237,242,247,.42));box-shadow:inset 0 1px 0 rgba(255,255,255,.95)"><span class="bar-fill" style="width:'+width+'%;background:linear-gradient(90deg,rgba(35,116,171,.96),rgba(27,138,131,.82));box-shadow:0 2px 8px rgba(35,116,171,.20)"></span><i title="'+esc(bubbleTitle)+'" style="position:absolute;left:calc('+bubbleLeft+'% - '+(size/2)+'px);top:50%;transform:translateY(-50%);display:block;width:'+size+'px;height:'+size+'px;border-radius:50%;background:radial-gradient(circle at 31% 27%,rgba(255,255,255,.96) 0 7%,rgba(255,241,203,.76) 19%,rgba(219,167,70,.48) 48%,rgba(211,145,54,.14) 73%,rgba(211,145,54,.05) 100%);border:1px solid rgba(211,145,54,.36);box-shadow:0 0 0 5px rgba(211,145,54,.055),0 9px 22px rgba(98,69,25,.18),inset -5px -6px 10px rgba(112,72,13,.08);backdrop-filter:blur(2px);cursor:help"></i></span><span class="bar-value" style="line-height:1.2"><b style="color:var(--ink)">'+n(contigs)+'</b><br><small>'+n(support)+'</small></span></div>';
 }});
 if(!evidenceRows) evidenceRows='<p class="muted">暂无 reads 支持且可归入所选分类层级的记录。</p>';
 panels[0].innerHTML='<div class="panel-head"><div><h2>主要'+esc(DATA.rank_label)+'：contig 与测序证据</h2><p class="muted">横条长度为候选 contig 数；气泡锚定在横条末端，面积表示 '+supportLabelV2()+'。全部分类单元可在面板内滚动查看。</p></div><span class="pill">全部 '+n(taxa.length)+'</span></div><div class="chart" style="max-height:680px;overflow:auto;padding:5px 14px 10px 0">'+evidenceRows+'</div>'+explanationV2('横条与气泡','每一横条的长度仅由候选 contig 成员数决定；末端气泡的面积独立编码 mapped reads。两类证据不相互替代，且批次汇总只用于定位需要查看的样本，不比较样本间丰度。');
 var priorityRows='';
 batch.priority_samples.forEach(function(x,index){{
   var sample=sampleMap[x.sample_id], attention=attentionCountV3(sample), support=DATA.batch.has_read_counts?x.read_count:x.relative_abundance;
   priorityRows+='<li><button data-sample="'+esc(x.sample_id)+'"><span class="rank-number">'+(index+1)+'</span><span>'+esc(x.sample_id)+'<br><small class="muted">'+(attention?'<span style="color:#815513;font-weight:700">[关注科 ×'+n(attention)+']</span> · ':'')+supportLabelV2()+' '+n(support)+' · 质量支持 '+n(x.quality_supported_count)+'</small></span><span class="pill">查看</span></button></li>';
 }});
 if(!priorityRows) priorityRows='<li class="muted">暂无可排序的样本结果。</li>';
 panels[1].innerHTML='<div class="panel-head"><div><h2>优先查看样本</h2><p class="muted">先显示含“关注科”注释证据的样本；同类样本再按 mapped reads、质量支持和检出数排序。</p></div><span class="pill">全部样本</span></div><ol class="rank-list" style="max-height:620px;overflow:auto">'+priorityRows+'</ol>'+explanationV2('导航规则','方括号中的关注科仅来自本地可维护清单，用于压缩和排序报告入口；它不等同于宿主范围、感染性、致病性或公共卫生风险。');
 var orderedMatrix=[].concat(batch.matrix_taxa).sort(function(a,b){{var order={{DNA:0,RNA:1,RT:2,MIXED:3,UNCLASSIFIED:4}};return (order[a.genome_group]??9)-(order[b.genome_group]??9)||a.taxon.localeCompare(b.taxon);}});
 var matrixRows='';
 orderedMatrix.forEach(function(t){{
   matrixRows+='<tr><td><span class="pill" style="background:'+genomeColourV3(t.genome_group)+'18;color:'+genomeColourV3(t.genome_group)+';border:1px solid '+genomeColourV3(t.genome_group)+'55">'+esc(genomeLabelV3(t))+'</span></td><th title="'+esc(t.taxon)+'">'+esc(t.taxon)+'</th>';
   DATA.samples.forEach(function(s){{var hit=t.samples.find(function(x){{return x.sample_id===s.sample_id;}});matrixRows+='<td>'+(hit?'<button class="present" data-sample="'+esc(s.sample_id)+'" title="'+esc(s.sample_id)+'：'+hit.count+' 个本地 vOTU；'+supportLabelV2()+' '+n(hit.read_count||hit.top_relative_abundance)+'"></button>':'<button disabled title="未检出"></button>')+'</td>';}});matrixRows+='</tr>';
 }});
 panels[2].innerHTML='<div class="panel-head"><div><h2>'+esc(DATA.rank_label)+' × 样本检测矩阵</h2><p class="muted">按 DNA、RNA、逆转录、家族内异质和未归类分组；点击方格进入对应样本证据页。</p></div><span class="pill">全部 '+n(orderedMatrix.length)+'</span></div><div class="matrix-wrap" style="max-height:760px"><table class="matrix"><thead><tr><th>ICTV 基因组组</th><th>'+esc(DATA.rank_label)+'</th>'+DATA.samples.map(function(s){{return '<th>'+esc(s.sample_id)+'</th>';}}).join('')+'</tr></thead><tbody>'+matrixRows+'</tbody></table></div>'+explanationV2('ICTV 基因组组','标签由报告随附的版本化 ICTV MSL 科级参考表提供。家族内出现跨组 Genome 记录时显示为“家族内异质”，未匹配时显示为“未归类”；两者均不由 contig 自动猜测。');
 document.querySelectorAll('[data-sample]').forEach(function(x){{x.addEventListener('click',function(){{openSample(x.dataset.sample);}});}});
}};
var renderSampleV3Base=renderSample;
renderSample=function(id){{
 renderSampleV3Base(id);
 var sample=sampleMap[id], grids=document.querySelectorAll('#sample .grid'), topPanel=grids[0].children[0], scatterPanel=grids[1].children[0];
 var rows=sample.rows, topRows=rows.filter(function(x){{return x.detected;}}).slice().sort(function(a,b){{return b.relative_abundance-a.relative_abundance||b.read_count-a.read_count;}}).slice(0,20);
 var maxAbundance=Math.max(1,...topRows.map(function(x){{return x.relative_abundance;}})), maxSupport=Math.max(1,...topRows.map(function(x){{return DATA.batch.has_read_counts?x.read_count:x.mean_coverage;}})), lollipops='';
 topRows.forEach(function(x){{
   var support=DATA.batch.has_read_counts?x.read_count:x.mean_coverage, size=bubbleSizeV2(support,maxSupport), position=x.relative_abundance/maxAbundance*100, colour=genomeColourV3(x.genome_group), border=x.priority_family?'3px solid #d39136':'2px solid #fff';
   var tooltip=x.votu_id+' | 病毒科：'+x.family+' | 属：'+x.genus+' | ICTV 组：'+genomeLabelV3(x)+' | 分类来源：'+x.taxonomy_source+' | CheckV：'+x.checkv_quality+' | mapped reads：'+n(x.read_count);
   lollipops+='<div title="'+esc(tooltip)+'" style="display:grid;grid-template-columns:minmax(132px,1fr) minmax(140px,1.1fr) minmax(160px,3fr) 75px;gap:10px;align-items:center;margin:16px 0;cursor:help"><span class="bar-name">'+esc(x.votu_id)+'<br><small class="source">'+esc(x.checkv_quality)+'</small></span><span class="source"><b>'+esc(x.family)+'</b><br><span style="color:'+colour+'">'+esc(genomeLabelV3(x))+'</span> · '+esc(x.taxonomy_source)+'</span><span class="bar-track" style="height:4px;overflow:visible;position:relative"><span class="bar-fill" style="height:4px;width:'+position+'%;background:'+colour+'"></span><i style="position:absolute;left:calc('+position+'% - '+(size/2)+'px);top:-'+((size-4)/2)+'px;display:block;width:'+size+'px;height:'+size+'px;border-radius:50%;background:'+colour+';border:'+border+';box-shadow:0 1px 5px rgba(0,0,0,.22)"></i></span><span class="bar-value">'+pct(x.relative_abundance)+'<br><small>'+n(support)+'</small></span></div>';
 }});
 if(!lollipops) lollipops='<p class="muted">无 reads 支持的本地 vOTU。</p>';
 topPanel.innerHTML='<div class="panel-head"><div><h2>Top 本地 vOTU：丰度、分类与测序支持</h2><p class="muted">横向位置精确表示本样本相对丰度；圆点面积表示 '+(DATA.batch.has_read_counts?'mapped reads':'平均覆盖度（旧结果回退）')+'；颜色表示 ICTV 基因组组，金色外圈表示“关注科”。</p></div></div><div class="chart" style="max-height:680px;overflow:auto;padding-right:8px">'+lollipops+'</div>'+explanationV2('点线图编码','每一行都是该样本内独立聚类和回贴的本地 vOTU。悬停可读取病毒科、属、分类来源、CheckV 标签和回贴支持。不同样本间的位置和丰度不用于比较。');
 var scatterRows=rows.filter(function(x){{return x.completeness!==null&&x.covered_fraction!==null;}}), maxReads=Math.max(1,...scatterRows.map(function(x){{return x.read_count||x.mean_coverage||0;}})), xTicks=[0,25,50,75,100].map(function(v){{return '<span style="position:absolute;bottom:-21px;left:'+v+'%;transform:translateX(-50%);font-size:10px;color:var(--muted)">'+v+'%</span>';}}).join(''), yTicks=[0,25,50,75,100].map(function(v){{return '<span style="position:absolute;left:-36px;bottom:'+v+'%;transform:translateY(50%);font-size:10px;color:var(--muted)">'+v+'%</span>';}}).join(''), dots='';
 scatterRows.forEach(function(x){{
   var support=x.read_count||x.mean_coverage||0, size=Math.max(8,Math.min(28,8+20*Math.sqrt(support/maxReads))), colour=qualityColourV3(x.checkv_quality), tooltip=x.votu_id+' | 完整度：'+x.completeness+'% | 覆盖碱基比例：'+x.covered_fraction+'% | 覆盖碱基：'+n(x.covered_bases)+'/'+n(x.representative_length)+' | mapped reads：'+n(x.read_count)+' | '+x.checkv_quality;
   dots+='<span class="dot" style="width:'+size+'px;height:'+size+'px;background:'+colour+';left:'+Math.max(1,Math.min(99,x.completeness))+'%;bottom:'+Math.max(1,Math.min(99,x.covered_fraction))+'%;border:'+(x.priority_family?'3px solid #d39136':'1px solid #fff')+'" title="'+esc(tooltip)+'"></span>';
 }});
 scatterPanel.innerHTML='<div class="panel-head"><div><h2>完整度与覆盖证据</h2><p class="muted">两轴均为 0–100% 的原始百分比：横轴为 CheckV 完整度，纵轴为代表序列被 reads 覆盖的碱基比例；气泡面积为 mapped reads，颜色为 CheckV 质量。</p></div></div><div class="scatter"><span class="axis-y">覆盖碱基比例 (%)</span><span class="axis-x">CheckV 完整度 (%)</span>'+xTicks+yTicks+dots+'</div>'+explanationV2('坐标与证据边界','覆盖碱基比例按 CoverM covered_bases / representative_length 计算，数值截断至 0–100%。它与 CheckV 完整度共同辅助选择复核对象，但不证明病毒活性、感染状态或病毒学功能。');
}};

+// A single delegated handler survives every panel redraw and works inside Streamlit's iframe.
function reportHashSample(){{return location.hash.startsWith('#sample=')?decodeURIComponent(location.hash.slice(8)):'';}}
function navigateToSampleV4(sampleId){{if(sampleMap[sampleId]&&(current!==sampleId||!$('#sample').classList.contains('active'))){{openSample(sampleId);}}}}
document.addEventListener('click',function(event){{
 var target=event.target.closest?event.target.closest('[data-sample]'):null;
 if(!target||target.disabled) return;
 event.preventDefault();
 navigateToSampleV4(target.dataset.sample);
}});
window.addEventListener('hashchange',function(){{
 var sampleId=reportHashSample();
 if(sampleId&&sampleMap[sampleId]&&sampleId!==current) navigateToSampleV4(sampleId);
}});

$('#created').textContent=DATA.metadata.created_at;$('#sample-select').innerHTML=DATA.samples.map(x=>`<option value="${{esc(x.sample_id)}}">${{esc(x.sample_id)}}</option>`).join('');$('#sample-select').value=current;renderOverview();renderSample(current);$('#overview-btn').onclick=openOverview;$('#sample-btn').onclick=()=>openSample(current);$('#sample-select').onchange=e=>openSample(e.target.value);if(location.hash.startsWith('#sample='))openSample(current);
</script></body></html>"""


def redirect_page(title: str, target: str) -> str:
    escaped = html.escape(target, quote=True)
    return f"""<!doctype html><html lang='zh-CN'><head><meta charset='utf-8'><meta http-equiv='refresh' content='0;url={escaped}'><title>{html.escape(title)}</title></head><body><p>正在打开交互式报告；若未自动跳转，请<a href='{escaped}'>点击此处</a>。</p></body></html>"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--groups-file", type=Path)
    parser.add_argument("--overview-rank", choices=["family", "genus", "species"], default="family")
    parser.add_argument("--theme", choices=["quarto-scientific"], default="quarto-scientific")
    parser.add_argument("--top-taxa", type=int, default=0, help="0 means display all classified taxa")
    parser.add_argument("--priority-families", default="", help="Comma-separated family names for navigation priority")
    parser.add_argument(
        "--ictv-reference", type=Path,
        default=Path(__file__).resolve().parents[2] / "config" / "ictv_family_genome_reference.tsv",
        help="Versioned ICTV family-level Genome reference TSV",
    )
    parser.add_argument(
        "--priority-reference", type=Path,
        default=Path(__file__).resolve().parents[2] / "config" / "priority_review_taxa.tsv",
        help="Editable attention-taxon TSV used only for report navigation",
    )
    args = parser.parse_args()
    if args.top_taxa < 0:
        raise SystemExit("--top-taxa must be zero or a positive integer")

    root = args.output_dir
    report_dir = root / "reports"
    data_dir = report_dir / "data"
    sample_data_dir = data_dir / "samples"
    sample_report_dir = report_dir / "samples"
    sample_data_dir.mkdir(parents=True, exist_ok=True)
    sample_report_dir.mkdir(parents=True, exist_ok=True)

    samples = [row["sample_id"] for row in read_tsv(args.manifest) if row.get("sample_id")]
    sample_groups = load_sample_groups(args.groups_file, samples)
    status = {row.get("sample_id", ""): row.get("status", "NOT_RUN") for row in read_tsv(root / "04_sample_votu" / "sample_status.tsv")}
    diamond = read_diamond_annotations(root / "02_diamond_nr_taxonomy" / "contig_taxonomy_lca.tsv")
    priority_families = normalise_family_set(args.priority_families)
    ictv_reference = load_ictv_reference(args.ictv_reference)
    priority_taxa = load_priority_reference(args.priority_reference, priority_families)
    payloads: dict[str, dict[str, Any]] = {}
    for sample in samples:
        rows = [serialise_row(enrich_annotation(row, diamond)) for row in read_tsv(root / "04_sample_votu" / sample / "votu_summary.tsv")]
        payloads[sample] = sample_payload(sample, status.get(sample, "NOT_RUN"), rows, priority_taxa, ictv_reference)

    batch = build_batch_payload(samples, payloads, args.overview_rank, args.top_taxa, sample_groups)
    rank_label = {"family": "病毒科", "genus": "病毒属", "species": "病毒种"}[args.overview_rank]
    metadata = {
        "report_schema_version": SCHEMA_VERSION,
        "created_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "report_theme": args.theme,
        "overview_rank": args.overview_rank,
        "overview_rank_label": rank_label,
        "top_taxa": args.top_taxa,
        "priority_taxa": sorted(record["taxon_name"] for record in priority_taxa.values()),
        "ictv_reference": {
            **source_fingerprint(args.ictv_reference),
            "matched_family_count": len(ictv_reference),
            "description": "版本化 ICTV MSL 科级 Genome 参考；未匹配或跨组分类单元不从 contig 自动推断。",
        },
        "priority_reference": {
            **source_fingerprint(args.priority_reference),
            "description": "可维护的关注分类单元导航清单；仅用于优先复核，不构成风险或宿主范围结论。",
        },
        "source_files": {
            "manifest": str(args.manifest),
            "groups_file": str(args.groups_file) if args.groups_file else "",
            "sample_status": str(root / "04_sample_votu" / "sample_status.tsv"),
            "local_votu_root": str(root / "04_sample_votu"),
            "diamond_taxonomy": str(root / "02_diamond_nr_taxonomy" / "contig_taxonomy_lca.tsv"),
            "ictv_reference": str(args.ictv_reference),
            "priority_reference": str(args.priority_reference),
        },
        "evidence_boundary": "批次页仅汇总每个样本中 reads 支持的本地 vOTU 及其分类线索；不进行跨样本丰度比较，不将不同样本本地 vOTU 视为同一病毒株。关注病毒科仅用于优先复核导航，不等同于人畜共患或致病风险结论。",
    }
    dashboard_data = {"metadata": metadata, "rank": args.overview_rank, "rank_label": rank_label, "top_taxa": args.top_taxa, "batch": batch, "samples": [payloads[sample] for sample in samples]}

    (data_dir / "report_metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (data_dir / "batch_dashboard_data.json").write_text(json.dumps(dashboard_data, ensure_ascii=False) + "\n", encoding="utf-8")
    for sample, sample_data in payloads.items():
        (sample_data_dir / f"{safe_id(sample)}.json").write_text(json.dumps(sample_data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        target = f"../virome_dashboard.html#sample={quote(sample)}"
        (sample_report_dir / f"{safe_id(sample)}.html").write_text(redirect_page(f"样本 {sample} 病毒判读", target), encoding="utf-8")

    write_tsv(report_dir / "batch_taxa_presence.tsv", batch["presence_rows"], ["rank", "taxon", "sample_id", "group", "detected_local_votu_count", "top_relative_abundance", "read_count", "candidate_contig_count", "high_quality_count", "report_path"])
    write_tsv(report_dir / "sample_report_index.tsv", batch["samples"], ["sample_id", "group", "status", "local_votu_count", "detected_local_votu_count", "quality_supported_count", "priority_count", "priority_family_votu_count", "read_count"])
    write_tsv(data_dir / "batch_summary.tsv", batch["samples"], ["sample_id", "group", "status", "local_votu_count", "detected_local_votu_count", "quality_supported_count", "priority_count", "priority_family_votu_count", "read_count"])
    write_tsv(data_dir / "evidence_legend.tsv", [
        {"field": "reads_supported_local_votu", "meaning": "样本内本地 vOTU representative 获得 clean reads 回贴支持", "not_evidence_for": "绝对丰度、病毒活性、感染或致病性"},
        {"field": "checkv_quality", "meaning": "CheckV 的病毒基因组质量评估", "not_evidence_for": "独立实验验证或宿主范围"},
        {"field": "taxonomy", "meaning": "DIAMOND/TaxonKit LCA 或 geNomad 的分类线索", "not_evidence_for": "最终物种鉴定"},
        {"field": "batch_taxon_presence", "meaning": "分类单元在多个样本中出现的导航汇总", "not_evidence_for": "跨样本生态比较或同一病毒株认定"},
    ], ["field", "meaning", "not_evidence_for"])

    dashboard = dashboard_html(dashboard_data)
    (report_dir / "virome_dashboard.html").write_text(dashboard, encoding="utf-8")
    (report_dir / "batch_overview.html").write_text(redirect_page("批次病毒检出总览", "virome_dashboard.html#overview"), encoding="utf-8")
    (report_dir / "index.html").write_text(redirect_page("病毒筛查报告中心", "virome_dashboard.html#overview"), encoding="utf-8")
    print(f"[INFO] Offline dashboard written: {report_dir / 'virome_dashboard.html'}")
    print(f"[INFO] Machine-readable report data: {data_dir}")


if __name__ == "__main__":
    main()
