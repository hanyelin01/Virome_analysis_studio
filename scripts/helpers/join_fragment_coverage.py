#!/usr/bin/env python3
"""Join CoverM output to per-sample final-fragment annotations."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path


def rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0: return []
    with path.open(encoding="utf-8", errors="replace", newline="") as handle: return list(csv.DictReader(handle, delimiter="\t"))


def number(row: dict[str, str], needle: str) -> str:
    for key, value in row.items():
        if needle.lower() in key.lower(): return value
    return "0"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotations", required=True, type=Path)
    parser.add_argument("--coverage", required=True, type=Path)
    parser.add_argument("--counts", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    coverage: dict[str, dict[str, str]] = {}
    for row in rows(args.coverage):
        genome = row.get("Genome") or row.get("genome") or ""
        if genome and genome.lower() != "unmapped": coverage[Path(genome).stem] = row
    counts: dict[str, dict[str, str]] = {}
    for row in rows(args.counts):
        genome = row.get("Genome") or row.get("genome") or ""
        if genome and genome.lower() != "unmapped": counts[Path(genome).stem] = row
    annotations = rows(args.annotations)
    fields = list(annotations[0]) if annotations else ["vf_id"]
    fields += [name for name in ("relative_abundance", "mean_coverage", "covered_bases", "read_count", "detected") if name not in fields]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader()
        for row in annotations:
            value = coverage.get(row.get("vf_id", ""), {})
            count = counts.get(row.get("vf_id", ""), {})
            row.update({"relative_abundance": number(value, "relative abundance"), "mean_coverage": number(value, "mean"), "covered_bases": number(value, "covered bases"), "read_count": number(count, "count")})
            try: row["detected"] = "yes" if float(row["read_count"]) > 0 else "no"
            except ValueError: row["detected"] = "no"
            writer.writerow(row)


if __name__ == "__main__": main()
