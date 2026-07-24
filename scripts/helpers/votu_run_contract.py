#!/usr/bin/env python3
"""Create or verify the effective parameter contract for the vOTU stage."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--min-length", type=int, required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--ani", required=True)
    parser.add_argument("--aligned-fraction", required=True)
    parser.add_argument("--read-identity", required=True)
    parser.add_argument("--read-aligned-percent", required=True)
    parser.add_argument("--covered-fraction", required=True)
    parser.add_argument("--importance-abundance", required=True)
    args = parser.parse_args()
    expected = {
        "schema_version": 1,
        "manifest_sha256": sha256(args.manifest),
        "checkv_candidates_sha256": sha256(args.input),
        "min_length": args.min_length,
        "threads": args.threads,
        "votu_ani": str(args.ani),
        "votu_aligned_fraction": str(args.aligned_fraction),
        "coverm_min_read_percent_identity": str(args.read_identity),
        "coverm_min_read_aligned_percent": str(args.read_aligned_percent),
        "coverm_min_covered_fraction": str(args.covered_fraction),
        "importance_relative_abundance": str(args.importance_abundance),
    }
    if args.contract.is_file():
        existing = json.loads(args.contract.read_text(encoding="utf-8"))
        if existing != expected:
            changed = [
                key for key in expected if existing.get(key) != expected.get(key)
            ]
            raise SystemExit(
                "vOTU resume contract differs for: "
                + ", ".join(changed)
                + ". Use a new report output directory."
            )
        print("[INFO] vOTU parameter contract matches the existing run")
        return
    args.contract.write_text(
        json.dumps(expected, indent=2) + "\n", encoding="utf-8"
    )
    print(f"[INFO] vOTU parameter contract written: {args.contract}")


if __name__ == "__main__":
    main()
