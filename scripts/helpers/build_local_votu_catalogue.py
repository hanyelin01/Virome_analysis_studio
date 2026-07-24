#!/usr/bin/env python3
"""Create a reproducible representative FASTA from Vclust cluster membership."""
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
                    yield header.split()[0], "".join(parts)
                header, parts = line[1:].strip(), []
            else:
                parts.append(line)
    if header is not None:
        yield header.split()[0], "".join(parts)


def read_clusters(path: Path, sequence_ids: set[str], singletons: bool) -> dict[str, str]:
    if singletons:
        return {sequence_id: sequence_id for sequence_id in sequence_ids}
    membership: dict[str, str] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            member, representative = fields[0].strip(), fields[1].strip()
            if member.lower() in {"object", "sequence_id"}:
                continue
            if member in sequence_ids and representative in sequence_ids:
                membership[member] = representative
    return {sequence_id: membership.get(sequence_id, sequence_id) for sequence_id in sequence_ids}


def write_fasta(handle, sequence_id: str, sequence: str) -> None:
    handle.write(f">{sequence_id}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start:start + 80] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--clusters", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--singletons", action="store_true")
    args = parser.parse_args()
    sequences = dict(fasta_records(args.input))
    if not sequences:
        raise SystemExit("No sequences available for local vOTU catalogue")
    if not args.singletons and (args.clusters is None or not args.clusters.is_file()):
        raise SystemExit("--clusters is required unless --singletons is used")
    membership = read_clusters(args.clusters, set(sequences), args.singletons)
    representatives = sorted(set(membership.values()), key=lambda item: (-len(sequences[item]), item))
    votu_for_representative = {representative: f"vOTU_{index:06d}" for index, representative in enumerate(representatives, 1)}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    representative_dir = args.output_dir / "representatives"
    representative_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "votu_cluster_members.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["votu_id", "representative_sequence_id", "member_sequence_id", "is_representative", "member_length"], delimiter="\t")
        writer.writeheader()
        for member, representative in sorted(membership.items(), key=lambda item: (votu_for_representative[item[1]], item[0])):
            writer.writerow({
                "votu_id": votu_for_representative[representative],
                "representative_sequence_id": representative,
                "member_sequence_id": member,
                "is_representative": int(member == representative),
                "member_length": len(sequences[member]),
            })
    with (args.output_dir / "votu_representative_map.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["votu_id", "representative_sequence_id", "representative_length", "member_count"], delimiter="\t")
        writer.writeheader()
        for representative in representatives:
            votu_id = votu_for_representative[representative]
            members = [member for member, rep in membership.items() if rep == representative]
            writer.writerow({"votu_id": votu_id, "representative_sequence_id": representative, "representative_length": len(sequences[representative]), "member_count": len(members)})
            with (representative_dir / f"{votu_id}.fna").open("w", encoding="utf-8") as fasta:
                write_fasta(fasta, votu_id, sequences[representative])
    print(f"[INFO] Created {len(representatives)} local vOTU representative sequence(s)")


if __name__ == "__main__":
    main()
