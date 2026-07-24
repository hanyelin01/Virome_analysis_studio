#!/usr/bin/env python3
"""Validate and standardise user-supplied candidate contig FASTA input."""
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


FASTA_SUFFIXES = {".fa", ".fna", ".fasta"}


def records(path: Path):
    header: str | None = None
    sequence: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header, sequence = line[1:].strip(), []
            else:
                sequence.append(line)
    if header is not None:
        yield header, "".join(sequence)


def input_files(source: Path, kind: str) -> list[Path]:
    if kind == "file":
        if not source.is_file():
            raise SystemExit(f"Custom FASTA is missing or not a file: {source}")
        return [source]
    if not source.is_dir():
        raise SystemExit(f"Custom FASTA directory is missing or not a directory: {source}")
    files = sorted(path for path in source.iterdir() if path.is_file() and path.suffix.lower() in FASTA_SUFFIXES)
    if not files:
        raise SystemExit(f"No .fa/.fna/.fasta file found in custom input directory: {source}")
    return files


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--input-type", required=True, choices=["file", "directory"])
    parser.add_argument("--output-fasta", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    args = parser.parse_args()

    files = input_files(args.input, args.input_type)
    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    seen: set[str] = set()
    count = 0
    with args.output_fasta.open("w", encoding="utf-8") as fasta, args.manifest.open("w", newline="", encoding="utf-8") as manifest:
        writer = csv.DictWriter(manifest, fieldnames=["source_file", "sha256", "original_id", "sequence_id", "length"], delimiter="\t")
        writer.writeheader()
        for source in files:
            checksum = digest(source)
            for header, sequence in records(source):
                original_id = header.split()[0]
                if not original_id:
                    raise SystemExit(f"Empty FASTA identifier in: {source}")
                safe_id = "".join(char if char.isalnum() or char in "._-" else "_" for char in original_id)
                if safe_id in seen:
                    raise SystemExit(f"Duplicate FASTA identifier after standardisation: {safe_id}")
                if not sequence:
                    raise SystemExit(f"Empty sequence for {original_id} in: {source}")
                seen.add(safe_id)
                fasta.write(f">{safe_id}\n")
                for start in range(0, len(sequence), 80):
                    fasta.write(sequence[start:start + 80] + "\n")
                writer.writerow({"source_file": str(source.resolve()), "sha256": checksum, "original_id": original_id, "sequence_id": safe_id, "length": len(sequence)})
                count += 1
    if not count:
        raise SystemExit("No FASTA record was found in the custom input")
    print(f"[INFO] Prepared {count} custom candidate contig(s): {args.output_fasta}")


if __name__ == "__main__":
    main()
