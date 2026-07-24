#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"
MANIFEST='' OUTPUT=''
while (($#)); do
  case "$1" in --manifest|--output) (($#>=2)) || die 2 "Missing value for $1"; [[ $1 == --manifest ]] && MANIFEST=$2 || OUTPUT=$2; shift 2;; *) die 2 "Unknown option: $1";; esac
done
[[ -f $MANIFEST && -n $OUTPUT ]] || die 2 '--manifest and --output are required'
mkdir -p "$(dirname "$OUTPUT")"; FAILED="$(dirname "$OUTPUT")/failed_samples.tsv"
printf 'sample_id\tstatus\tcontig_file\tcontig_count\tfile_size_bytes\n' > "$OUTPUT"; printf 'sample_id\tstatus\tcontig_file\n' > "$FAILED"; status=0
while IFS= read -r manifest_row; do
  IFS=$'\x1f' read -r sample type raw1 raw2 clean1 clean2 clean_single assembly <<< "${manifest_row//$'\t'/$'\x1f'}"
  contig="$assembly/final.contigs.fa"; count=0; size=0
  if [[ -s $contig ]]; then count=$(grep -c '^>' "$contig" || true); size=$(stat -c '%s' "$contig"); fi
  if [[ -s $contig && $count -gt 0 ]]; then result=SUCCESS; else result=FAILED; printf '%s\t%s\t%s\n' "$sample" "$result" "$contig" >> "$FAILED"; status=5; fi
  printf '%s\t%s\t%s\t%s\t%s\n' "$sample" "$result" "$contig" "$count" "$size" | tee -a "$OUTPUT"
done < <(manifest_rows "$MANIFEST")
exit "$status"
