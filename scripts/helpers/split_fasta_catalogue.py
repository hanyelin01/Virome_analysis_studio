#!/usr/bin/env python3
"""Split viral candidates into single-record FASTA files for CoverM genome mode."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def records(path: Path):
    header = None
    seq: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:].strip(), []
            else:
                seq.append(line)
    if header is not None:
        yield header, "".join(seq)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--mapping", required=True, type=Path)
    args = parser.parse_args()
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise SystemExit(f"Refusing to overwrite non-empty catalogue directory: {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with args.mapping.open("w", newline="", encoding="utf-8") as map_handle:
        writer = csv.DictWriter(map_handle, fieldnames=["candidate_id", "sequence_id", "length"], delimiter="\t")
        writer.writeheader()
        count = 0
        for count, (header, sequence) in enumerate(records(args.input), 1):
            candidate_id = f"candidate_{count:07d}"
            file = args.output_dir / f"{candidate_id}.fna"
            with file.open("w", encoding="utf-8") as out:
                out.write(f">{header.split()[0]}\n")
                for start in range(0, len(sequence), 80):
                    out.write(sequence[start:start + 80] + "\n")
            writer.writerow({"candidate_id": candidate_id, "sequence_id": header.split()[0], "length": len(sequence)})
    if not count:
        raise SystemExit("No sequences available for catalogue splitting")
    print(f"[INFO] Prepared {count} single-sequence files for vOTU clustering")


if __name__ == "__main__":
    main()
