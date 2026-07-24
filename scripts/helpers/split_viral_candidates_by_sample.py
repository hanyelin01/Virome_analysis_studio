#!/usr/bin/env python3
"""Restore sample ownership for CheckV candidates while retaining provenance."""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


def fasta_records(path: Path):
    header: str | None = None
    parts: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header, parts = line[1:].strip(), []
            else:
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)


def read_tsv(path: Path, key: str) -> dict[str, dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return {}
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        return {row.get(key, ""): row for row in csv.DictReader(handle, delimiter="\t") if row.get(key, "")}


def lookup(table: dict[str, dict[str, str]], sequence_id: str) -> dict[str, str]:
    return table.get(sequence_id) or table.get(sequence_id.split("|", 1)[0]) or {}


def resolve_provenance(
    table: dict[str, dict[str, str]], sequence_id: str
) -> tuple[str, dict[str, str], str]:
    """Resolve exact IDs first, then CheckV's provirus child suffix (_1, _2, ...)."""
    direct_id = sequence_id.split("|", 1)[0]
    if direct_id in table:
        return direct_id, table[direct_id], "exact"
    match = re.fullmatch(r"(.+)_([1-9][0-9]*)", direct_id)
    if match and match.group(1) in table:
        parent_id = match.group(1)
        return parent_id, table[parent_id], "checkv_provirus_suffix"
    return "", {}, "unresolved"


def write_fasta(handle, sequence_id: str, sequence: str) -> None:
    handle.write(f">{sequence_id}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start:start + 80] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--provenance", required=True, type=Path)
    parser.add_argument("--quality", required=True, type=Path)
    parser.add_argument("--taxonomy", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--min-length", required=True, type=int)
    args = parser.parse_args()

    provenance = read_tsv(args.provenance, "sequence_id")
    quality = read_tsv(args.quality, "contig_id")
    taxonomy = read_tsv(args.taxonomy, "seq_name")
    manifest_rows = read_tsv(args.manifest, "sample_id")
    samples = list(manifest_rows)
    args.output_root.mkdir(parents=True, exist_ok=True)

    fasta_handles: dict[str, object] = {}
    metadata_handles: dict[str, object] = {}
    writers: dict[str, csv.DictWriter] = {}
    fieldnames = [
        "sequence_id", "provenance_sequence_id", "provenance_match", "sample_id",
        "length", "checkv_quality", "miuvig_quality",
        "completeness", "contamination", "taxonomy", "virus_score",
    ]
    try:
        for sample in samples:
            candidate_dir = args.output_root / sample / "01_candidates"
            candidate_dir.mkdir(parents=True, exist_ok=True)
            fasta_handles[sample] = (candidate_dir / "viral_candidates_checkv.fna").open("w", encoding="utf-8")
            metadata_handles[sample] = (candidate_dir / "candidate_metadata.tsv").open("w", newline="", encoding="utf-8")
            writers[sample] = csv.DictWriter(metadata_handles[sample], fieldnames=fieldnames, delimiter="\t")
            writers[sample].writeheader()

        kept = Counter()
        short = Counter()
        unassigned: list[dict[str, str]] = []
        for header, sequence in fasta_records(args.input):
            sequence_id = header.split()[0]
            provenance_id, source, provenance_match = resolve_provenance(
                provenance, sequence_id
            )
            if provenance_match == "checkv_provirus_suffix":
                parent_quality = lookup(quality, provenance_id)
                if parent_quality.get("provirus", "").strip().lower() not in {"yes", "true"}:
                    provenance_id, source, provenance_match = "", {}, "unresolved"
            sample = source.get("sample_id", "")
            if sample not in writers:
                unassigned.append({
                    "sequence_id": sequence_id,
                    "length": str(len(sequence)),
                    "reason": "missing exact/CheckV-parent provenance or sample not in manifest",
                })
                continue
            if len(sequence) < args.min_length:
                short[sample] += 1
                continue
            quality_row = lookup(quality, provenance_id)
            taxonomy_row = lookup(taxonomy, provenance_id)
            write_fasta(fasta_handles[sample], sequence_id, sequence)
            writers[sample].writerow({
                "sequence_id": sequence_id,
                "provenance_sequence_id": provenance_id,
                "provenance_match": provenance_match,
                "sample_id": sample,
                "length": len(sequence),
                "checkv_quality": quality_row.get("checkv_quality", "Not-determined"),
                "miuvig_quality": quality_row.get("miuvig_quality", "Genome-fragment"),
                "completeness": quality_row.get("completeness", "NA"),
                "contamination": quality_row.get("contamination", "NA"),
                "taxonomy": taxonomy_row.get("taxonomy", "Unclassified virus"),
                "virus_score": taxonomy_row.get("virus_score", "NA"),
            })
            kept[sample] += 1
    finally:
        for handle in fasta_handles.values():
            handle.close()
        for handle in metadata_handles.values():
            handle.close()

    with (args.output_root / "split_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "kept_candidates", "excluded_short_candidates"], delimiter="\t")
        writer.writeheader()
        for sample in samples:
            writer.writerow({"sample_id": sample, "kept_candidates": kept[sample], "excluded_short_candidates": short[sample]})
    with (args.output_root / "unassigned_candidates.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sequence_id", "length", "reason"], delimiter="\t")
        writer.writeheader()
        writer.writerows(unassigned)
    if unassigned:
        raise SystemExit(f"{len(unassigned)} CheckV candidates could not be restored to a manifest sample; inspect unassigned_candidates.tsv")
    completion = {
        "schema_version": 1,
        "manifest_sample_count": len(samples),
        "kept_candidate_count": sum(kept.values()),
        "excluded_short_candidate_count": sum(short.values()),
        "unassigned_candidate_count": 0,
    }
    (args.output_root / "split_complete.json").write_text(
        json.dumps(completion, indent=2) + "\n", encoding="utf-8"
    )
    print(f"[INFO] Restored {sum(kept.values())} CheckV candidates to {len(samples)} sample directories")


if __name__ == "__main__":
    main()
