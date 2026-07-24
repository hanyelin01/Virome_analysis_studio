#!/usr/bin/env python3
"""Merge FASTA files while preserving the first occurrence of each record ID."""
from __future__ import annotations

import argparse
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
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("inputs", nargs="+", type=Path)
    args = parser.parse_args()
    seen: set[str] = set()
    count = 0
    with args.output.open("w", encoding="utf-8") as out:
        for input_file in args.inputs:
            if not input_file.is_file():
                continue
            for header, sequence in records(input_file):
                record_id = header.split()[0]
                if record_id in seen:
                    continue
                seen.add(record_id)
                out.write(f">{record_id}\n")
                for start in range(0, len(sequence), 80):
                    out.write(sequence[start:start + 80] + "\n")
                count += 1
    if not count:
        raise SystemExit("No FASTA records were written")
    print(f"[INFO] Merged {count} unique candidate sequence(s): {args.output}")


if __name__ == "__main__":
    main()
