#!/usr/bin/env python3
"""Convert DIAMOND outfmt-6 hits into query-wise TaxonKit LCA input."""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import OrderedDict
from pathlib import Path


def configure_csv_field_limit() -> None:
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit //= 10


configure_csv_field_limit()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hits", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    values: OrderedDict[str, set[str]] = OrderedDict()
    with args.hits.open(encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 12:
                continue
            query = row[0]
            taxids = re.findall(r"\d+", row[11])
            if taxids:
                values.setdefault(query, set()).update(taxids)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["query_id", "taxids"])
        for query, taxids in values.items():
            writer.writerow([query, ",".join(sorted(taxids, key=int))])
    print(f"[INFO] TaxonKit LCA input: {len(values)} query sequence(s): {args.output}")


if __name__ == "__main__":
    main()
