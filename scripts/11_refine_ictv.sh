#!/usr/bin/env bash
# Compare CheckV-refined fragments to the locally versioned ICTV reference DB.
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"; load_pipeline_config

INPUT='' OUTPUT_DIR='' THREADS='' BLOCK_SIZE="$DIAMOND_BLOCK_SIZE" INDEX_CHUNKS="$DIAMOND_INDEX_CHUNKS" TMPDIR="$DIAMOND_TMPDIR" RESUME=0
while (($#)); do
  case "$1" in
    --input|--output-dir|--threads|--block-size|--index-chunks|--tmpdir)
      (($# >= 2)) || die 2 "Missing value for $1"
      case "$1" in --input) INPUT=$2;; --output-dir) OUTPUT_DIR=$2;; --threads) THREADS=$2;; --block-size) BLOCK_SIZE=$2;; --index-chunks) INDEX_CHUNKS=$2;; --tmpdir) TMPDIR=$2;; esac
      shift 2;;
    --resume) RESUME=1; shift;;
    -h|--help) echo "Usage: $0 --input CHECKV_FASTA --output-dir PATH --threads N [--block-size GB] [--index-chunks N] [--tmpdir PATH] [--resume]"; exit 0;;
    *) die 2 "Unknown option: $1";;
  esac
done
[[ -s $INPUT && -n $OUTPUT_DIR && -n $THREADS ]] || die 2 'Input FASTA, output directory and threads are required'
positive_int "$THREADS" && positive_int "$ICTV_REFERENCE_MAX_TARGET_SEQS" && (( THREADS <= MAX_THREADS_PER_VIRAL_TOOL && THREADS <= MAX_TOTAL_THREADS )) || die 2 'Invalid ICTV DIAMOND thread request or max target count'
validate_diamond_performance "$THREADS" "$BLOCK_SIZE" "$INDEX_CHUNKS" "$TMPDIR"
[[ -n $ICTV_REFERENCE_DMND && -f $ICTV_REFERENCE_DMND ]] || die 3 'ICTV_REFERENCE_DMND is missing or not a readable .dmnd file'
[[ -n $ICTV_REFERENCE_METADATA && -f $ICTV_REFERENCE_METADATA ]] || die 3 'ICTV_REFERENCE_METADATA is missing or not readable'
require_command diamond
OUT="$OUTPUT_DIR/06_ictv_refinement"; HITS="$OUT/ictv_hits.tsv"
if [[ -f $HITS ]]; then
  (( RESUME )) && { echo "[INFO] ICTV refinement hits already exist; skipped"; exit 0; }
  die 4 "ICTV refinement output already exists; use --resume or a new output directory"
fi
if [[ -e $OUT ]] && ! dir_is_empty_or_missing "$OUT"; then die 4 "ICTV refinement output is incomplete or conflicting: $OUT"; fi
mkdir -p "$OUT"
fields=(qseqid qlen qstart qend pident length evalue bitscore sseqid)
args=(blastx --db "$ICTV_REFERENCE_DMND" --query "$INPUT" --out "$HITS" --outfmt 6 "${fields[@]}" --threads "$THREADS" --block-size "$BLOCK_SIZE" --index-chunks "$INDEX_CHUNKS" -t "$TMPDIR" --evalue "$DIAMOND_EVALUE" --max-target-seqs "$ICTV_REFERENCE_MAX_TARGET_SEQS")
[[ $DIAMOND_SENSITIVITY == more-sensitive ]] && args+=(--more-sensitive)
args+=(--header)
printf '%q ' diamond "${args[@]}" > "$OUT/diamond_command.sh"; printf '\n' >> "$OUT/diamond_command.sh"
diamond "${args[@]}" >"$OUT/diamond.stdout.log" 2>"$OUT/diamond.stderr.log"
[[ -f $HITS ]] || die 1 "ICTV DIAMOND did not create a hit table; inspect: $OUT/diamond.stderr.log"
cp -- "$ICTV_REFERENCE_METADATA" "$OUT/ictv_reference_metadata.tsv"
printf 'input=%s\nreference_dmnd=%s\nreference_metadata=%s\nreference_version=%s\nthreads=%s\nblock_size=%s\nindex_chunks=%s\ntmpdir=%s\n' "$INPUT" "$ICTV_REFERENCE_DMND" "$ICTV_REFERENCE_METADATA" "$ICTV_REFERENCE_VERSION" "$THREADS" "$BLOCK_SIZE" "$INDEX_CHUNKS" "$TMPDIR" > "$OUT/parameters.env"
echo "[INFO] ICTV refinement hits: $HITS"
