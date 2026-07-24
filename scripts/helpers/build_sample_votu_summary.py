#!/usr/bin/env python3
"""Join local Vclust, CheckV/geNomad and CoverM results into one sample table."""
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value: str | None) -> float:
    try:
        return float(value or 0)
    except ValueError:
        return 0.0


def coverage_values(rows: list[dict[str, str]]) -> dict[str, dict[str, float]]:
    values: dict[str, dict[str, float]] = {}
    for row in rows:
        genome = row.get("Genome") or row.get("genome") or ""
        if not genome or genome.lower() == "unmapped":
            continue
        votu_id = Path(genome).stem
        relative = next((number(value) for key, value in row.items() if "relative abundance" in key.lower()), 0.0)
        mean = next((number(value) for key, value in row.items() if key.lower().endswith(" mean") or key.lower() == "mean"), 0.0)
        covered_bases = next((number(value) for key, value in row.items() if "covered bases" in key.lower()), 0.0)
        read_count = next((number(value) for key, value in row.items() if key.lower().endswith(" count") or key.lower() == "count"), 0.0)
        values[votu_id] = {"relative_abundance": relative, "mean_coverage": mean, "covered_bases": covered_bases, "read_count": read_count}
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--members", required=True, type=Path)
    parser.add_argument("--representatives", required=True, type=Path)
    parser.add_argument("--coverage", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--importance-abundance", type=float, default=5.0)
    parser.add_argument("--allow-empty", action="store_true")
    args = parser.parse_args()

    fields = [
        "sample_id", "votu_id", "representative_sequence_id", "representative_length", "member_count",
        "checkv_quality", "miuvig_quality", "completeness", "contamination", "taxonomy", "virus_score",
        "relative_abundance", "mean_coverage", "covered_bases", "read_count", "detected", "importance",
    ]
    metadata = {row.get("sequence_id", ""): row for row in read_tsv(args.metadata)}
    members = read_tsv(args.members)
    representatives = {row.get("votu_id", ""): row for row in read_tsv(args.representatives)}
    if not members and not args.allow_empty:
        raise SystemExit("Local vOTU membership table is empty")
    counts = Counter(row.get("votu_id", "") for row in members)
    coverage = coverage_values(read_tsv(args.coverage) if args.coverage else [])
    rows: list[dict[str, object]] = []
    for votu_id, representative in sorted(representatives.items()):
        sequence_id = representative.get("representative_sequence_id", "")
        metadata_row = metadata.get(sequence_id, {})
        abundance = coverage.get(votu_id, {})
        relative = abundance.get("relative_abundance", 0.0)
        quality = metadata_row.get("checkv_quality", "Not-determined")
        high_quality = quality.lower() in {"complete", "high-quality", "medium-quality"}
        if high_quality and relative >= args.importance_abundance:
            importance = "高优先级"
        elif high_quality or relative >= args.importance_abundance:
            importance = "关注"
        else:
            importance = "常规"
        rows.append({
            "sample_id": args.sample_id,
            "votu_id": votu_id,
            "representative_sequence_id": sequence_id,
            "representative_length": representative.get("representative_length", metadata_row.get("length", "")),
            "member_count": counts[votu_id],
            "checkv_quality": quality,
            "miuvig_quality": metadata_row.get("miuvig_quality", "Genome-fragment"),
            "completeness": metadata_row.get("completeness", "NA"),
            "contamination": metadata_row.get("contamination", "NA"),
            "taxonomy": metadata_row.get("taxonomy", "Unclassified virus"),
            "virus_score": metadata_row.get("virus_score", "NA"),
            "relative_abundance": f"{relative:.8f}",
            "mean_coverage": f"{abundance.get('mean_coverage', 0.0):.8f}",
            "covered_bases": f"{abundance.get('covered_bases', 0.0):.0f}",
            "read_count": f"{abundance.get('read_count', 0.0):.0f}",
            "detected": "yes" if relative > 0 else "no",
            "importance": importance,
        })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"[INFO] Wrote {len(rows)} local vOTU summary row(s): {args.output}")


if __name__ == "__main__":
    main()
