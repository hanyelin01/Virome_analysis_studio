#!/usr/bin/env python3
"""Validate the adapter catalogue and turn fastp JSON into auditable evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

DNA = re.compile(r"^[ACGT]+$")
REQUIRED = {
    "profile_id", "display_name", "platform", "library_kits", "r1_sequence",
    "r2_sequence", "source_level", "source_title", "source_url",
    "source_version", "verified_on", "status", "notes",
}
REFERENCE_REQUIRED = {
    "sequence_id", "aliases", "sequence", "reverse_complement_id", "category",
    "trimming_action", "source_level", "source_title", "source_url",
    "verified_on", "status", "notes",
}


def load_catalog(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or not REQUIRED.issubset(rows[0]):
        raise ValueError(f"catalogue columns are incomplete: {path}")
    seen: set[str] = set()
    for line_no, row in enumerate(rows, 2):
        profile_id = row["profile_id"].strip()
        if not re.fullmatch(r"[a-z0-9][a-z0-9_-]*", profile_id):
            raise ValueError(f"line {line_no}: invalid profile_id")
        if profile_id in seen:
            raise ValueError(f"line {line_no}: duplicate profile_id {profile_id}")
        seen.add(profile_id)
        for key in ("r1_sequence", "r2_sequence"):
            sequence = row[key].strip().upper()
            if sequence and (len(sequence) < 6 or not DNA.fullmatch(sequence)):
                raise ValueError(f"line {line_no}: invalid {key}")
            row[key] = sequence
        if profile_id != "auto" and not (row["r1_sequence"] and row["r2_sequence"]):
            raise ValueError(f"line {line_no}: a manual profile requires R1 and R2")
        if row["source_level"] not in {"vendor", "software", "community"}:
            raise ValueError(f"line {line_no}: invalid source_level")
        if row["status"] not in {"active", "review", "retired"}:
            raise ValueError(f"line {line_no}: invalid status")
        if not row["source_url"].startswith("https://"):
            raise ValueError(f"line {line_no}: source_url must use HTTPS")
    return rows


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def load_reference(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or not REFERENCE_REQUIRED.issubset(rows[0]):
        raise ValueError(f"reference columns are incomplete: {path}")
    by_id: dict[str, dict[str, str]] = {}
    for line_no, row in enumerate(rows, 2):
        sequence_id = row["sequence_id"].strip()
        sequence = row["sequence"].strip().upper()
        if not re.fullmatch(r"[a-z0-9][a-z0-9_-]*", sequence_id):
            raise ValueError(f"reference line {line_no}: invalid sequence_id")
        if sequence_id in by_id:
            raise ValueError(f"reference line {line_no}: duplicate sequence_id")
        if len(sequence) < 6 or not DNA.fullmatch(sequence):
            raise ValueError(f"reference line {line_no}: invalid sequence")
        if row["status"] not in {"active", "review", "rejected", "retired"}:
            raise ValueError(f"reference line {line_no}: invalid status")
        row["sequence"] = sequence
        by_id[sequence_id] = row
    sequences: dict[str, str] = {}
    for row in rows:
        sequence = row["sequence"]
        if sequence in sequences:
            raise ValueError(
                f"duplicate reference sequence: {row['sequence_id']} and {sequences[sequence]}"
            )
        sequences[sequence] = row["sequence_id"]
        rc_id = row["reverse_complement_id"]
        if rc_id:
            if rc_id not in by_id:
                raise ValueError(f"{row['sequence_id']}: missing reverse complement {rc_id}")
            if reverse_complement(sequence) != by_id[rc_id]["sequence"]:
                raise ValueError(f"{row['sequence_id']}: reverse complement mismatch")
    return rows


def profile(rows: list[dict[str, str]], profile_id: str) -> dict[str, str]:
    match = next((row for row in rows if row["profile_id"] == profile_id), None)
    if match is None:
        raise ValueError(f"unknown adapter profile: {profile_id}")
    if match["status"] == "retired":
        raise ValueError(f"adapter profile is retired: {profile_id}")
    return match


def sequence_match(sequence: str, rows: list[dict[str, str]]) -> tuple[str, str]:
    sequence = sequence.upper()
    if not sequence:
        return "", "not_reported"
    exact: list[str] = []
    family: list[str] = []
    for row in rows:
        for field in ("r1_sequence", "r2_sequence"):
            known = row[field]
            if not known:
                continue
            if sequence == known:
                exact.append(row["profile_id"])
            elif min(len(sequence), len(known)) >= 12 and (
                sequence.startswith(known) or known.startswith(sequence)
            ):
                family.append(row["profile_id"])
    if exact:
        return ";".join(sorted(set(exact))), "exact_catalogue_match"
    if family:
        return ";".join(sorted(set(family))), "prefix_family_match"
    return "", "not_in_catalogue"


def reference_match(sequence: str, rows: list[dict[str, str]]) -> str:
    if not sequence:
        return ""
    matches = []
    for row in rows:
        known = row["sequence"]
        if sequence == known or (
            min(len(sequence), len(known)) >= 12
            and (sequence.startswith(known) or known.startswith(sequence))
        ):
            matches.append(f"{row['sequence_id']}[{row['status']}]")
    return ";".join(matches)


def summarize(args: argparse.Namespace) -> None:
    rows = load_catalog(args.catalog)
    reference_rows = load_reference(args.reference)
    selected = profile(rows, args.profile)
    data = json.loads(args.fastp_json.read_text(encoding="utf-8"))
    cutting = data.get("adapter_cutting") or {}
    detected = {
        "R1": str(cutting.get("read1_adapter_sequence") or ""),
        "R2": str(cutting.get("read2_adapter_sequence") or ""),
    }
    total_before = int(
        ((data.get("summary") or {}).get("before_filtering") or {}).get("total_reads") or 0
    )
    trimmed = int(cutting.get("adapter_trimmed_reads") or 0)
    evidence_rows = []
    for read, sequence in detected.items():
        matched, judgement = sequence_match(sequence, rows)
        method = "fastp_auto_detection" if args.profile == "auto" else "configured_fallback_plus_fastp_auto"
        source_row = selected
        if args.profile == "auto" and matched:
            first_match = matched.split(";", 1)[0]
            source_row = profile(rows, first_match)
        if args.profile != "auto":
            expected = selected["r1_sequence" if read == "R1" else "r2_sequence"]
            if sequence == expected:
                judgement = "reported_sequence_consistent_with_selected_profile"
                matched = args.profile
            elif sequence:
                judgement = "reported_sequence_differs_from_selected_profile"
            else:
                judgement = "selected_profile_used_as_fallback_sequence"
        evidence_rows.append({
            "sample_id": args.sample,
            "read": read,
            "selected_profile": args.profile,
            "selected_profile_name": selected["display_name"],
            "fastp_reported_sequence": sequence or "未报告",
            "matched_catalogue_profiles": matched or "未匹配",
            "matched_reference_sequences": reference_match(sequence, reference_rows) or "未匹配",
            "source_judgement": judgement,
            "evidence_method": method,
            "adapter_trimmed_reads": trimmed,
            "input_reads": total_before,
            "trimmed_read_fraction": f"{trimmed / total_before:.6f}" if total_before else "",
            "source_level": source_row["source_level"],
            "source_title": source_row["source_title"],
            "source_url": source_row["source_url"],
            "catalog_sha256": hashlib.sha256(args.catalog.read_bytes()).hexdigest(),
        })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=evidence_rows[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(evidence_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    validate = sub.add_parser("validate")
    validate.add_argument("--catalog", type=Path, required=True)
    validate.add_argument("--reference", type=Path)
    get = sub.add_parser("get")
    get.add_argument("--catalog", type=Path, required=True)
    get.add_argument("--profile", required=True)
    get.add_argument("--field", choices=sorted(REQUIRED), required=True)
    report = sub.add_parser("summarize")
    report.add_argument("--catalog", type=Path, required=True)
    report.add_argument("--reference", type=Path, required=True)
    report.add_argument("--profile", required=True)
    report.add_argument("--fastp-json", type=Path, required=True)
    report.add_argument("--sample", required=True)
    report.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = load_catalog(args.catalog)
    if args.command == "validate":
        message = f"validated {len(rows)} adapter profiles"
        if args.reference:
            message += f" and {len(load_reference(args.reference))} reference sequences"
        print(message)
    elif args.command == "get":
        print(profile(rows, args.profile)[args.field])
    else:
        summarize(args)


if __name__ == "__main__":
    main()
