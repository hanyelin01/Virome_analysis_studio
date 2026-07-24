#!/usr/bin/env python3
"""Create a traceable, length-filtered contig catalogue from MEGAHIT outputs."""
from __future__ import annotations

import argparse
import csv
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--min-length", required=True, type=int)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fasta_out = args.output_dir / "merged_assembled_contigs.fna"
    provenance_out = args.output_dir / "contig_provenance.tsv"
    summary_out = args.output_dir / "contig_preparation_summary.tsv"
    assemblies = sorted(args.assembly_dir.glob("*/final.contigs.fa"))
    if not assemblies:
        raise SystemExit(f"No final.contigs.fa files found under: {args.assembly_dir}")

    used: set[str] = set()
    total_records = kept_records = 0
    summaries: list[dict[str, object]] = []
    with fasta_out.open("w", encoding="utf-8") as fasta, provenance_out.open("w", newline="", encoding="utf-8") as prov:
        writer = csv.DictWriter(prov, fieldnames=["sequence_id", "sample_id", "original_contig_id", "length"] , delimiter="\t")
        writer.writeheader()
        for assembly in assemblies:
            sample = assembly.parent.name
            sample_total = sample_kept = 0
            for original_header, sequence in fasta_records(assembly):
                total_records += 1
                sample_total += 1
                original_id = original_header.split()[0]
                safe_original = "".join(char if char.isalnum() or char in "._-" else "_" for char in original_id)
                seq_id = f"{sample}__{safe_original}"
                if seq_id in used:
                    raise SystemExit(f"Duplicate contig identifier after standardisation: {seq_id}")
                used.add(seq_id)
                if len(sequence) < args.min_length:
                    continue
                kept_records += 1
                sample_kept += 1
                fasta.write(f">{seq_id}\n")
                for start in range(0, len(sequence), 80):
                    fasta.write(sequence[start:start + 80] + "\n")
                writer.writerow({"sequence_id": seq_id, "sample_id": sample, "original_contig_id": original_id, "length": len(sequence)})
            summaries.append({"sample_id": sample, "assembled_contigs": sample_total, "kept_contigs": sample_kept})
    with summary_out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "assembled_contigs", "kept_contigs"], delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)
    if kept_records == 0:
        raise SystemExit(f"No contigs >= {args.min_length} bp were retained")
    print(f"[INFO] Prepared {kept_records}/{total_records} contigs from {len(assemblies)} sample(s): {fasta_out}")


if __name__ == "__main__":
    main()
