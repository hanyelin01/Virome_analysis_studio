#!/usr/bin/env python3
"""Join CheckV, NR and ICTV evidence, then distribute final fragments by sample."""
from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0: return []
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fasta_records(path: Path):
    header: str | None = None; parts: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line: continue
            if line.startswith(">"):
                if header is not None: yield header, "".join(parts)
                header, parts = line[1:].strip(), []
            else: parts.append(line)
    if header is not None: yield header, "".join(parts)


def write_fasta(handle, ident: str, sequence: str) -> None:
    handle.write(f">{ident}\n")
    for offset in range(0, len(sequence), 80): handle.write(sequence[offset:offset + 80] + "\n")


def ictv_hit_rows(path: Path):
    """Read both historical headerless and current headered DIAMOND outfmt-6."""
    fields = ["qseqid", "qlen", "qstart", "qend", "pident", "length", "evalue", "bitscore", "sseqid"]
    if not path.is_file() or path.stat().st_size == 0:
        return
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        for values in csv.reader(handle, delimiter="\t"):
            if not values or not values[0].strip() or values[0].strip() == "qseqid":
                continue
            if len(values) < len(fields):
                continue
            yield dict(zip(fields, values))


def parent_vc(ident: str, known: set[str]) -> str:
    ident = ident.split()[0].split("|", 1)[0]
    if ident in known: return ident
    match = re.fullmatch(r"(.+)_([1-9][0-9]*)", ident)
    return match.group(1) if match and match.group(1) in known else ""


def taxon_at_rank(annotation: dict[str, str], rank: str) -> str:
    if annotation.get("lca_rank", "").lower() == rank and annotation.get("lca_name"): return annotation["lca_name"]
    lineage = annotation.get("lca_lineage", "")
    prefix = {"order": "o", "family": "f", "genus": "g", "species": "s"}[rank]
    match = re.search(rf"(?:^|;)\s*(?:{prefix}__|{rank}[:=_ ]+)\s*([^;]+)", lineage, re.I)
    return match.group(1).strip() if match else ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkv-fasta", required=True, type=Path)
    parser.add_argument("--checkv-quality", required=True, type=Path)
    parser.add_argument("--decision", required=True, type=Path)
    parser.add_argument("--nr-taxonomy", required=True, type=Path)
    parser.add_argument("--source-mapping", required=True, type=Path)
    parser.add_argument("--ictv-hits", required=True, type=Path)
    parser.add_argument("--ictv-metadata", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--catalogue-dir", required=True, type=Path)
    parser.add_argument("--sample-dir", required=True, type=Path)
    args = parser.parse_args()

    decision = {row["vc_id"]: row for row in read_tsv(args.decision) if row.get("vc_id")}
    taxonomy = {row["query_id"]: row for row in read_tsv(args.nr_taxonomy) if row.get("query_id")}
    quality = {row["contig_id"]: row for row in read_tsv(args.checkv_quality) if row.get("contig_id")}
    by_vc_sources: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(args.source_mapping):
        if row.get("vc_id") and row.get("sample_id"): by_vc_sources[row["vc_id"]].append(row)
    ictv_meta = {row.get("reference_id", row.get("ref_id", "")): row for row in read_tsv(args.ictv_metadata)}
    best_ictv: dict[str, list[str]] = {}
    for row in ictv_hit_rows(args.ictv_hits):
        raw_query = row.get("qseqid", "").strip()
        if not raw_query:
            continue
        query = raw_query.split()[0]
        if query and query not in best_ictv: best_ictv[query] = [row.get("sseqid", ""), row.get("pident", ""), row.get("length", ""), row.get("evalue", ""), row.get("bitscore", "")]
    samples = [row.get("sample_id", "") for row in read_tsv(args.manifest) if row.get("sample_id")]
    known = set(decision)
    args.catalogue_dir.mkdir(parents=True, exist_ok=True); args.sample_dir.mkdir(parents=True, exist_ok=True)
    refs = args.catalogue_dir / "references"; refs.mkdir(exist_ok=True)
    fields = ["vf_id", "checkv_sequence_id", "parent_vc_id", "length", "decision", "checkv_quality", "completeness", "contamination", "nr_order", "nr_family", "nr_genus", "nr_species", "nr_lca_name", "ictv_reference_id", "ictv_species", "ictv_genus", "ictv_family", "baltimore_group", "ictv_pident", "ictv_alignment_length", "ictv_evalue", "ictv_bitscore"]
    all_rows: list[dict[str, str]] = []
    sequences_by_vf: dict[str, str] = {}
    with (args.catalogue_dir / "VF_catalogue.fna").open("w", encoding="utf-8") as final_fasta, \
         (args.catalogue_dir / "VF_catalogue.tsv").open("w", newline="", encoding="utf-8") as table:
        writer = csv.DictWriter(table, fieldnames=fields, delimiter="\t"); writer.writeheader()
        for index, (header, sequence) in enumerate(fasta_records(args.checkv_fasta), 1):
            checkv_id = header.split()[0]; vf_id = f"VF_{index:07d}"; vc_id = parent_vc(checkv_id, known)
            if not vc_id: raise SystemExit(f"Cannot restore CheckV fragment to a VC parent: {checkv_id}")
            q = quality.get(checkv_id) or quality.get(vc_id, {})
            nr = taxonomy.get(vc_id, {}); hit = best_ictv.get(checkv_id, ["", "", "", "", ""]); meta = ictv_meta.get(hit[0], {})
            row = {"vf_id": vf_id, "checkv_sequence_id": checkv_id, "parent_vc_id": vc_id, "length": str(len(sequence)),
                   "decision": decision[vc_id].get("decision", ""), "checkv_quality": q.get("checkv_quality", "Not-determined"),
                   "completeness": q.get("completeness", "NA"), "contamination": q.get("contamination", "NA"),
                   "nr_order": taxon_at_rank(nr, "order"), "nr_family": taxon_at_rank(nr, "family"),
                   "nr_genus": taxon_at_rank(nr, "genus"), "nr_species": taxon_at_rank(nr, "species"), "nr_lca_name": nr.get("lca_name", ""),
                   "ictv_reference_id": hit[0], "ictv_species": meta.get("species", ""), "ictv_genus": meta.get("genus", ""),
                   "ictv_family": meta.get("family", ""), "baltimore_group": meta.get("baltimore_group", "UNCLASSIFIED"),
                   "ictv_pident": hit[1], "ictv_alignment_length": hit[2], "ictv_evalue": hit[3], "ictv_bitscore": hit[4]}
            writer.writerow(row); all_rows.append(row); sequences_by_vf[vf_id] = sequence; write_fasta(final_fasta, vf_id, sequence)
            with (refs / f"{vf_id}.fna").open("w", encoding="utf-8") as ref: write_fasta(ref, vf_id, sequence)
    sample_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    by_parent = {row["vf_id"]: row for row in all_rows}
    for vf_id, row in by_parent.items():
        for source in by_vc_sources.get(row["parent_vc_id"], []): sample_rows[source["sample_id"]].append(row)
    for sample in samples:
        root = args.sample_dir / sample; root.mkdir(parents=True, exist_ok=True)
        rows = {row["vf_id"]: row for row in sample_rows[sample]}
        with (root / "viral_fragments.tsv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader(); writer.writerows(rows.values())
        sample_refs = root / "references"; sample_refs.mkdir(exist_ok=True)
        for vf_id in rows:
            with (sample_refs / f"{vf_id}.fna").open("w", encoding="utf-8") as ref:
                write_fasta(ref, vf_id, sequences_by_vf[vf_id])
    with (args.sample_dir / "sample_fragment_presence.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "vf_id", "parent_vc_id"], delimiter="\t"); writer.writeheader()
        for sample in samples:
            for row in {row["vf_id"]: row for row in sample_rows[sample]}.values(): writer.writerow({"sample_id": sample, "vf_id": row["vf_id"], "parent_vc_id": row["parent_vc_id"]})
    print(f"[INFO] Final global catalogue: {len(all_rows)} CheckV-refined fragment(s); distributed to {len(samples)} sample(s)")


if __name__ == "__main__": main()
