#!/usr/bin/env python3
"""Create a traceable, length-filtered contig catalogue from MEGAHIT outputs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
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


def manifest_assemblies(manifest: Path, assembly_root: Path) -> list[tuple[str, Path]]:
    root = assembly_root.resolve()
    with manifest.open(encoding="utf-8", errors="replace", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or "sample_id" not in rows[0] or "assembly_dir" not in rows[0]:
        raise SystemExit(f"Manifest lacks sample_id/assembly_dir rows: {manifest}")
    result: list[tuple[str, Path]] = []
    seen: set[str] = set()
    for row in rows:
        sample = row["sample_id"].strip()
        assembly_dir = Path(row["assembly_dir"]).resolve()
        if not sample or sample in seen:
            raise SystemExit(f"Manifest contains an empty or duplicate sample ID: {sample!r}")
        seen.add(sample)
        if assembly_dir != (root / sample).resolve():
            raise SystemExit(
                f"{sample}: manifest assembly path is inconsistent with the selected assembly root: {assembly_dir}"
            )
        fasta = assembly_dir / "final.contigs.fa"
        if not fasta.is_file():
            raise SystemExit(f"{sample}: final.contigs.fa is missing: {fasta}")
        result.append((sample, fasta))
    return result


def input_fingerprint(assemblies: list[tuple[str, Path]], min_length: int) -> dict[str, object]:
    inputs = []
    for sample, path in assemblies:
        stat = path.stat()
        inputs.append({
            "sample_id": sample,
            "assembly_fasta": str(path),
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
        })
    payload: dict[str, object] = {
        "schema_version": 1,
        "min_length": min_length,
        "inputs": inputs,
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["fingerprint_sha256"] = hashlib.sha256(canonical).hexdigest()
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly-dir", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--min-length", required=True, type=int)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fasta_out = args.output_dir / "merged_assembled_contigs.fna"
    provenance_out = args.output_dir / "contig_provenance.tsv"
    summary_out = args.output_dir / "contig_preparation_summary.tsv"
    fingerprint_out = args.output_dir / "preparation_inputs.json"
    assemblies = manifest_assemblies(args.manifest, args.assembly_dir)
    expected_fingerprint = input_fingerprint(assemblies, args.min_length)
    required_outputs = (fasta_out, provenance_out, summary_out, fingerprint_out)
    if args.resume and all(path.is_file() and path.stat().st_size > 0 for path in required_outputs):
        existing = json.loads(fingerprint_out.read_text(encoding="utf-8"))
        if existing == expected_fingerprint:
            print("[INFO] Prepared contigs match the current manifest and parameters; skipped")
            return
        raise SystemExit(
            "Existing prepared contigs were created from different samples, assembly files, or "
            "minimum-length settings. Use a new report output directory."
        )
    if args.resume and any(path.exists() for path in required_outputs):
        raise SystemExit(
            "Existing prepared-contig output predates manifest fingerprinting or is incomplete; "
            "it cannot be safely resumed. Use a new report output directory."
        )
    if any(path.exists() for path in required_outputs):
        raise SystemExit(f"Prepared-contig output already exists: {args.output_dir}")

    used: set[str] = set()
    total_records = kept_records = 0
    summaries: list[dict[str, object]] = []
    with fasta_out.open("w", encoding="utf-8") as fasta, provenance_out.open("w", newline="", encoding="utf-8") as prov:
        writer = csv.DictWriter(prov, fieldnames=["sequence_id", "sample_id", "original_contig_id", "length"] , delimiter="\t")
        writer.writeheader()
        for sample, assembly in assemblies:
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
    fingerprint_out.write_text(
        json.dumps(expected_fingerprint, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"[INFO] Prepared {kept_records}/{total_records} contigs from {len(assemblies)} sample(s): {fasta_out}")


if __name__ == "__main__":
    main()
