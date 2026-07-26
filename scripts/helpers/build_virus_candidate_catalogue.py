#!/usr/bin/env python3
"""Build a provenance-preserving, exact-deduplicated viral candidate catalogue.

The catalogue deliberately does not create vOTUs.  A ``VC`` identifier denotes
one exact nucleotide sequence; every discovery call and every source sample is
retained in accompanying TSV files.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from collections import defaultdict
from pathlib import Path


def configure_csv_field_limit() -> None:
    """Allow DIAMOND's potentially long slineages field without overflow."""
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit //= 10


configure_csv_field_limit()


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
                    yield header, "".join(parts).upper()
                header, parts = line[1:].strip(), []
            else:
                parts.append(line)
    if header is not None:
        yield header, "".join(parts).upper()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parent_id(value: str) -> str:
    """Recover the prepared-contig ID from a VirSorter2-style suffix."""
    return value.split()[0].split("||", 1)[0].split("|", 1)[0]


def write_fasta(handle, ident: str, sequence: str) -> None:
    handle.write(f">{ident}\n")
    for offset in range(0, len(sequence), 80):
        handle.write(sequence[offset : offset + 80] + "\n")


def diamond_discovery_calls(path: Path, all_contigs: dict[str, str]):
    """Yield one best DIAMOND discovery call per contig, without loading hits."""
    seen: set[str] = set()
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 8:
                continue
            raw_id = row[0].split()[0]
            if raw_id in seen or raw_id not in all_contigs:
                continue
            seen.add(raw_id)
            yield {
                "tool": "DIAMOND-NR-virus", "raw_sequence_id": raw_id,
                "source_sequence_id": raw_id, "sequence": all_contigs[raw_id],
                "evidence": f"best_bitscore={row[7]};evalue={row[6]}",
            }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-fasta", required=True, type=Path)
    parser.add_argument("--provenance", required=True, type=Path)
    parser.add_argument("--genomad-fasta", required=True, type=Path)
    parser.add_argument("--virsorter-fasta", required=True, type=Path)
    parser.add_argument("--diamond-virus-hits", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    provenance = {row["sequence_id"]: row for row in read_tsv(args.provenance) if row.get("sequence_id")}
    all_contigs = {header.split()[0]: sequence for header, sequence in fasta_records(args.input_fasta)}
    calls: list[dict[str, str]] = []

    def add_fasta_calls(path: Path, tool: str) -> None:
        for header, sequence in fasta_records(path):
            raw_id = header.split()[0]
            source_id = parent_id(raw_id)
            calls.append({
                "tool": tool, "raw_sequence_id": raw_id, "source_sequence_id": source_id,
                "sequence": sequence, "evidence": "called",
            })

    add_fasta_calls(args.genomad_fasta, "geNomad")
    add_fasta_calls(args.virsorter_fasta, "VirSorter2")
    calls.extend(diamond_discovery_calls(args.diamond_virus_hits, all_contigs))

    if not calls:
        raise SystemExit("No discovery calls were supplied by geNomad, VirSorter2, or DIAMOND-NR-virus")

    by_hash: dict[str, dict[str, object]] = {}
    for call in calls:
        sequence = call["sequence"]
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        entry = by_hash.setdefault(digest, {"sequence": sequence, "calls": []})
        entry["calls"].append(call)

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    fasta = out / "VC_catalogue.fna"
    catalogue = out / "VC_catalogue.tsv"
    source_map = out / "VC_source_mapping.tsv"
    evidence = out / "VC_discovery_evidence.tsv"
    summary = out / "discovery_summary.tsv"
    vc_by_hash = {digest: f"VC_{index:07d}" for index, digest in enumerate(sorted(by_hash), 1)}

    with fasta.open("w", encoding="utf-8") as fasta_handle, \
         catalogue.open("w", newline="", encoding="utf-8") as catalogue_handle, \
         source_map.open("w", newline="", encoding="utf-8") as source_handle, \
         evidence.open("w", newline="", encoding="utf-8") as evidence_handle:
        catalogue_writer = csv.DictWriter(catalogue_handle, fieldnames=["vc_id", "sequence_sha256", "length", "source_call_count"], delimiter="\t")
        source_writer = csv.DictWriter(source_handle, fieldnames=["vc_id", "sample_id", "source_sequence_id", "original_contig_id", "tool", "raw_sequence_id", "evidence"], delimiter="\t")
        evidence_writer = csv.DictWriter(evidence_handle, fieldnames=["vc_id", "length", "geNomad", "VirSorter2", "DIAMOND_NR_virus", "discovery_pattern", "supporting_method_count"], delimiter="\t")
        catalogue_writer.writeheader(); source_writer.writeheader(); evidence_writer.writeheader()
        for digest in sorted(by_hash):
            entry = by_hash[digest]
            vc_id = vc_by_hash[digest]
            sequence = str(entry["sequence"])
            vc_calls = list(entry["calls"])
            write_fasta(fasta_handle, vc_id, sequence)
            catalogue_writer.writerow({"vc_id": vc_id, "sequence_sha256": digest, "length": len(sequence), "source_call_count": len(vc_calls)})
            tools = {str(call["tool"]) for call in vc_calls}
            for call in vc_calls:
                source = provenance.get(str(call["source_sequence_id"]), {})
                source_writer.writerow({
                    "vc_id": vc_id, "sample_id": source.get("sample_id", ""),
                    "source_sequence_id": call["source_sequence_id"],
                    "original_contig_id": source.get("original_contig_id", ""),
                    "tool": call["tool"], "raw_sequence_id": call["raw_sequence_id"],
                    "evidence": call["evidence"],
                })
            evidence_writer.writerow({
                "vc_id": vc_id, "length": len(sequence),
                "geNomad": "yes" if "geNomad" in tools else "no",
                "VirSorter2": "yes" if "VirSorter2" in tools else "no",
                "DIAMOND_NR_virus": "yes" if "DIAMOND-NR-virus" in tools else "no",
                "discovery_pattern": " + ".join(sorted(tools)),
                "supporting_method_count": len(tools),
            })

    counts: defaultdict[str, int] = defaultdict(int)
    for entry in by_hash.values():
        counts[" + ".join(sorted({str(call["tool"]) for call in entry["calls"]}))] += 1
    with summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["discovery_pattern", "vc_count"], delimiter="\t")
        writer.writeheader()
        for pattern, count in sorted(counts.items()):
            writer.writerow({"discovery_pattern": pattern, "vc_count": count})
    print(f"[INFO] Built {len(by_hash)} exact nonredundant VC record(s): {fasta}")


if __name__ == "__main__":
    main()
