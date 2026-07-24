#!/usr/bin/env python3
"""Select only CheckV fragments with a family-level NR classification for ICTV."""
from __future__ import annotations
import argparse
import csv
import re
from pathlib import Path


def rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0: return []
    with path.open(encoding="utf-8", errors="replace", newline="") as handle: return list(csv.DictReader(handle, delimiter="\t"))


def records(path: Path):
    header=None; parts=[]
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line=raw.strip()
            if not line: continue
            if line.startswith(">"):
                if header is not None: yield header, "".join(parts)
                header, parts=line[1:].strip(), []
            else: parts.append(line)
    if header is not None: yield header, "".join(parts)


def family(row: dict[str, str]) -> str:
    if row.get("lca_rank", "").lower() == "family": return row.get("lca_name", "")
    match = re.search(r"(?:^|;)\s*(?:f__|family[:=_ ]+)\s*([^;]+)", row.get("lca_lineage", ""), re.I)
    return match.group(1).strip() if match else ""


def parent(ident: str, known: set[str]) -> str:
    ident=ident.split()[0].split("|",1)[0]
    if ident in known: return ident
    match=re.fullmatch(r"(.+)_([1-9][0-9]*)", ident)
    return match.group(1) if match and match.group(1) in known else ""


def write(handle, ident: str, sequence: str) -> None:
    handle.write(f">{ident}\n")
    for i in range(0,len(sequence),80): handle.write(sequence[i:i+80]+"\n")


def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument("--checkv-fasta", required=True, type=Path)
    parser.add_argument("--decision", required=True, type=Path)
    parser.add_argument("--taxonomy", required=True, type=Path)
    parser.add_argument("--output-fasta", required=True, type=Path)
    parser.add_argument("--output-table", required=True, type=Path)
    args=parser.parse_args()
    decisions={row["vc_id"]: row for row in rows(args.decision) if row.get("vc_id")}
    taxonomy={row["query_id"]: row for row in rows(args.taxonomy) if row.get("query_id")}
    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    selected=[]
    with args.output_fasta.open("w", encoding="utf-8") as fasta:
        for header, sequence in records(args.checkv_fasta):
            vf_id=header.split()[0]; vc_id=parent(vf_id,set(decisions)); annotation=taxonomy.get(vc_id,{})
            taxon=family(annotation)
            if vc_id and taxon and decisions[vc_id].get("decision") in {"confirmed_viral","putative_novel_virus"}:
                write(fasta,vf_id,sequence); selected.append({"vf_id":vf_id,"parent_vc_id":vc_id,"nr_family":taxon})
    with args.output_table.open("w",newline="",encoding="utf-8") as handle:
        writer=csv.DictWriter(handle,fieldnames=["vf_id","parent_vc_id","nr_family"],delimiter="\t"); writer.writeheader(); writer.writerows(selected)
    print(f"[INFO] Selected {len(selected)} family-classified CheckV fragment(s) for ICTV refinement")


if __name__ == "__main__": main()
