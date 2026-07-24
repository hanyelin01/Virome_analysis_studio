#!/usr/bin/env python3
"""Create a compact standalone report for a custom fine-annotation run."""
from __future__ import annotations

import argparse
import csv
import html
from collections import Counter
from pathlib import Path


def rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--input-description", required=True)
    parser.add_argument("--taxon-scope", required=True)
    parser.add_argument("--taxonlist", default="all NR")
    args = parser.parse_args()
    root = args.output_dir
    report_dir = root / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    table = rows(root / "02_diamond_nr_taxonomy" / "contig_taxonomy_lca.tsv")
    classified = [row for row in table if row.get("lca_taxid") not in {"", "0"}]
    taxa = Counter((row.get("lca_name") or "Unclassified") for row in classified).most_common(20)
    rma = root / "01_diamond_megan" / "viral_candidates.nr.rma6"
    list_rows = "".join(f"<tr><td>{html.escape(name)}</td><td>{count}</td></tr>" for name, count in taxa) or "<tr><td colspan='2'>No TaxonKit LCA result was generated.</td></tr>"
    rma_text = f"RMA6 is available: <code>{html.escape(str(rma))}</code>" if rma.is_file() else "RMA6 was not requested or was not generated."
    body = f"""<!doctype html><html lang='en'><head><meta charset='utf-8'><title>Fine annotation report</title>
<style>body{{font-family:Arial,'Microsoft YaHei',sans-serif;background:#f5f8fb;color:#152535;margin:0}}main{{max-width:980px;margin:auto;padding:34px}}.hero,.card{{background:#fff;border:1px solid #dce7ef;border-radius:14px;padding:20px;margin:14px 0}}.hero{{border-top:5px solid #10a59b}}.grid{{display:grid;grid-template-columns:repeat(3,1fr);gap:12px}}strong{{display:block;font-size:28px;color:#126c88;margin-top:4px}}small{{color:#607080}}table{{border-collapse:collapse;width:100%}}th,td{{padding:9px;border-bottom:1px solid #e5edf2;text-align:left}}th{{background:#edf6f8}}code{{word-break:break-all}}</style>
</head><body><main><section class='hero'><h1>Custom contig fine-annotation report</h1><p>Input: {html.escape(args.input_description)}</p><p>NR search scope: <b>{html.escape(args.taxon_scope)}</b>; TaxID list: <b>{html.escape(args.taxonlist)}</b></p></section>
<section class='grid'><div class='card'><small>Sequences with DIAMOND hit</small><strong>{len(table)}</strong></div><div class='card'><small>TaxonKit LCA assigned</small><strong>{len(classified)}</strong></div><div class='card'><small>MEGAN</small><strong>{'RMA6 ready' if rma.is_file() else 'Not run'}</strong></div></section>
<section class='card'><h2>Top LCA assignments</h2><table><tr><th>Taxon</th><th>Contigs</th></tr>{list_rows}</table></section><section class='card'><h2>Files</h2><p>{rma_text}</p><p>Tabular annotations: <code>02_diamond_nr_taxonomy/contig_taxonomy_lca.tsv</code></p><p>Raw DIAMOND outfmt 6: <code>02_diamond_nr_taxonomy/nr_virus_hits.outfmt6.tsv</code></p></section>
</main></body></html>"""
    output = report_dir / "fine_annotation_report.html"
    output.write_text(body, encoding="utf-8")
    print(f"[INFO] Fine annotation report: {output}")


if __name__ == "__main__":
    main()
