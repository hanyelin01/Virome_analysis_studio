#!/usr/bin/env python3
"""Reconcile discovery and full-NR evidence without discarding novel viruses."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fasta_records(path: Path):
    header: str | None = None; parts: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line: continue
            if line.startswith(">"):
                if header is not None: yield header, "".join(parts)
                header, parts = line[1:].strip(), []
            else: parts.append(line)
    if header is not None: yield header, "".join(parts)


def write_fasta(handle, ident: str, sequence: str) -> None:
    handle.write(f">{ident}\n")
    for offset in range(0, len(sequence), 80): handle.write(sequence[offset:offset + 80] + "\n")


def is_viral(annotation: dict[str, str]) -> bool:
    text = " ".join(annotation.get(key, "") for key in ("lca_lineage", "best_lineages", "best_scientific_names")).lower()
    return "viruses" in text or "virus;" in text


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalogue-fasta", required=True, type=Path)
    parser.add_argument("--discovery-evidence", required=True, type=Path)
    parser.add_argument("--taxonomy", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    evidence = {row["vc_id"]: row for row in read_tsv(args.discovery_evidence) if row.get("vc_id")}
    taxonomy = {row["query_id"]: row for row in read_tsv(args.taxonomy) if row.get("query_id")}
    out = args.output_dir; out.mkdir(parents=True, exist_ok=True)
    decisions: list[dict[str, str]] = []
    sequences = list(fasta_records(args.catalogue_fasta))
    fields = ["vc_id", "decision", "reason", "supporting_method_count", "discovery_pattern", "nr_lca_name", "nr_lca_rank", "nr_lineage", "best_subject_id", "best_bitscore"]
    with (out / "viral_decision.tsv").open("w", newline="", encoding="utf-8") as table, \
         (out / "confirmed_viral.fna").open("w", encoding="utf-8") as confirmed, \
         (out / "putative_novel_virus.fna").open("w", encoding="utf-8") as novel, \
         (out / "checkv_input.fna").open("w", encoding="utf-8") as checkv_input:
        writer = csv.DictWriter(table, fieldnames=fields, delimiter="\t"); writer.writeheader()
        for header, sequence in sequences:
            vc_id = header.split()[0]
            discovery = evidence.get(vc_id, {})
            annotation = taxonomy.get(vc_id, {})
            support = int(discovery.get("supporting_method_count", "0") or 0)
            viral = is_viral(annotation)
            has_nr = bool(annotation.get("best_subject_id"))
            if viral:
                decision, reason = "confirmed_viral", "full NR taxonomy contains viral support"
                write_fasta(confirmed, vc_id, sequence); write_fasta(checkv_input, vc_id, sequence)
            elif not has_nr and support >= 2:
                decision, reason = "putative_novel_virus", "two or more discovery methods support a sequence without an NR hit"
                write_fasta(novel, vc_id, sequence); write_fasta(checkv_input, vc_id, sequence)
            elif has_nr and support >= 2:
                decision, reason = "ambiguous", "multi-tool discovery conflicts with nonviral full-NR evidence"
            else:
                decision, reason = "nonviral_or_insufficient", "insufficient viral evidence after full-NR review"
            row = {"vc_id": vc_id, "decision": decision, "reason": reason,
                   "supporting_method_count": str(support), "discovery_pattern": discovery.get("discovery_pattern", ""),
                   "nr_lca_name": annotation.get("lca_name", ""), "nr_lca_rank": annotation.get("lca_rank", ""),
                   "nr_lineage": annotation.get("lca_lineage", ""), "best_subject_id": annotation.get("best_subject_id", ""),
                   "best_bitscore": annotation.get("best_bitscore", "")}
            writer.writerow(row); decisions.append(row)
    print(f"[INFO] Viral decision table: {out / 'viral_decision.tsv'}; CheckV input candidates: {sum(row['decision'] in {'confirmed_viral', 'putative_novel_virus'} for row in decisions)}")


if __name__ == "__main__":
    main()
