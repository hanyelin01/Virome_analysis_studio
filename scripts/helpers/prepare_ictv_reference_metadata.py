#!/usr/bin/env python3
"""Validate a pinned ICTV-VMR-derived local reference metadata table.

Sequence retrieval is intentionally a separate, reviewable administrator step:
the VMR accession list and the downloaded proteins must be frozen together
before a production run.  This script produces the metadata contract consumed
by the ICTV refinement stage.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path


REQUIRED = {"reference_id", "family", "genus", "species", "baltimore_group"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True, type=Path, help="Reviewed TSV derived from one pinned ICTV VMR release")
    parser.add_argument("--protein-fasta", required=True, type=Path, help="Protein FASTA whose IDs equal metadata reference_id")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--version", required=True, help="For example VMR_MSL41.v1.20260721")
    args = parser.parse_args()
    with args.metadata.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        missing = REQUIRED - fields
        if missing: raise SystemExit(f"ICTV metadata is missing required columns: {', '.join(sorted(missing))}")
        records = list(reader)
    if not records: raise SystemExit("ICTV metadata has no reference records")
    identifiers = [row["reference_id"].strip() for row in records]
    if any(not value for value in identifiers) or len(identifiers) != len(set(identifiers)):
        raise SystemExit("ICTV reference_id values must be non-empty and unique")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / "ictv_reference_metadata.tsv"
    output.write_text(args.metadata.read_text(encoding="utf-8"), encoding="utf-8")
    protein_hash = hashlib.sha256(args.protein_fasta.read_bytes()).hexdigest()
    metadata_hash = hashlib.sha256(output.read_bytes()).hexdigest()
    (args.output_dir / "ictv_reference_manifest.json").write_text(json.dumps({
        "schema_version": 1, "ictv_reference_version": args.version,
        "metadata_file": output.name, "metadata_sha256": metadata_hash,
        "protein_fasta": str(args.protein_fasta), "protein_fasta_sha256": protein_hash,
        "reference_count": len(records),
    }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"[INFO] Validated {len(records)} ICTV reference records: {output}")


if __name__ == "__main__": main()
