#!/usr/bin/env python3
"""Sample FASTQ reads and report exact occurrences of versioned reference sequences."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path

from adapter_evidence import load_reference


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="ascii", errors="replace") if path.suffix == ".gz" else path.open(
        encoding="ascii", errors="replace"
    )


def scan_fastq(path: Path, references: list[dict[str, str]], max_reads: int) -> tuple[int, dict[str, list[int]]]:
    counts = {row["sequence_id"]: [0, 0, 0] for row in references}
    read_count = 0
    with open_text(path) as handle:
        while read_count < max_reads:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().strip().upper()
            plus = handle.readline()
            quality = handle.readline()
            if not (header.startswith("@") and plus.startswith("+") and quality):
                raise ValueError(f"invalid FASTQ record near read {read_count + 1}: {path}")
            read_count += 1
            for row in references:
                needle = row["sequence"]
                if needle not in sequence:
                    continue
                counts[row["sequence_id"]][0] += 1
                if sequence.startswith(needle):
                    counts[row["sequence_id"]][1] += 1
                if sequence.endswith(needle):
                    counts[row["sequence_id"]][2] += 1
    return read_count, counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--r1", type=Path, required=True)
    parser.add_argument("--r2", type=Path, required=True)
    parser.add_argument("--max-reads", type=int, default=100_000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.max_reads < 1:
        parser.error("--max-reads must be positive")
    references = load_reference(args.reference)
    output_rows = []
    for read_label, path in (("R1", args.r1), ("R2", args.r2)):
        total, counts = scan_fastq(path, references, args.max_reads)
        for row in references:
            anywhere, at_5p, at_3p = counts[row["sequence_id"]]
            if anywhere == 0:
                continue
            output_rows.append({
                "sample_id": args.sample,
                "read": read_label,
                "sequence_id": row["sequence_id"],
                "aliases": row["aliases"],
                "sequence": row["sequence"],
                "category": row["category"],
                "trimming_action": row["trimming_action"],
                "status": row["status"],
                "reads_scanned": total,
                "reads_with_sequence": anywhere,
                "match_fraction": f"{anywhere / total:.8f}" if total else "",
                "matches_at_5prime": at_5p,
                "matches_at_3prime": at_3p,
                "source_title": row["source_title"],
                "source_url": row["source_url"],
                "interpretation": (
                    "禁用：未核验序列" if row["status"] == "rejected"
                    else "协议特异引物；需按实验方案决定是否另行剪切"
                    if row["trimming_action"] == "protocol_specific"
                    else "仅作来源和污染判断，不自动加入剪切"
                ),
            })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_id", "read", "sequence_id", "aliases", "sequence", "category",
        "trimming_action", "status", "reads_scanned", "reads_with_sequence",
        "match_fraction", "matches_at_5prime", "matches_at_3prime",
        "source_title", "source_url", "interpretation",
    ]
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)


if __name__ == "__main__":
    main()
