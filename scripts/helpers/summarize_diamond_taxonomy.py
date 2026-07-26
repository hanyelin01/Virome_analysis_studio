#!/usr/bin/env python3
"""Join DIAMOND best hits with TaxonKit LCA/lineage output."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def tab_rows(path: Path):
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.reader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hits", required=True, type=Path)
    parser.add_argument("--lca", required=True, type=Path)
    parser.add_argument("--lineage", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    best: dict[str, list[str]] = {}
    for row in tab_rows(args.hits):
        if len(row) >= 13 and row[0] not in best:
            best[row[0]] = row
    lca: dict[str, list[str]] = {}
    for row in tab_rows(args.lca):
        if len(row) >= 3:
            lca[row[0]] = row
    lineage: dict[str, list[str]] = {}
    for row in tab_rows(args.lineage):
        if len(row) >= 4:
            lineage[row[0]] = row

    query_ids = list(dict.fromkeys([*best, *lca]))
    fields = ["query_id", "query_taxids", "lca_taxid", "lca_lineage", "lca_name", "lca_rank", "best_subject_id", "best_evalue", "best_bitscore", "best_staxids", "best_scientific_names", "best_lineages"]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for query in query_ids:
            lca_row, lineage_row, hit = lca.get(query, []), lineage.get(query, []), best.get(query, [])
            writer.writerow({
                "query_id": query,
                "query_taxids": lca_row[1] if len(lca_row) > 1 else "",
                "lca_taxid": lca_row[2] if len(lca_row) > 2 else "",
                "lca_lineage": lineage_row[3] if len(lineage_row) > 3 else "",
                "lca_name": lineage_row[4] if len(lineage_row) > 4 else "",
                "lca_rank": lineage_row[5] if len(lineage_row) > 5 else "",
                "best_subject_id": hit[10] if len(hit) > 10 else "",
                "best_evalue": hit[6] if len(hit) > 6 else "",
                "best_bitscore": hit[7] if len(hit) > 7 else "",
                "best_staxids": hit[11] if len(hit) > 11 else "",
                "best_scientific_names": hit[12] if len(hit) > 12 else "",
                # Lineage is derived by TaxonKit from best_staxids; DIAMOND
                # 2.0.x does not provide the later slineages output field.
                "best_lineages": "",
            })
    print(f"[INFO] DIAMOND/TaxonKit annotation table: {args.output}")


if __name__ == "__main__":
    main()
