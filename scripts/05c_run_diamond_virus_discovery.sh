#!/usr/bin/env bash
# Discovery-only DIAMOND search against the viral subset of a complete NR DB.
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
    -h|--help) echo "Usage: $0 --input FASTA --output-dir PATH --threads N [--block-size GB] [--index-chunks N] [--tmpdir PATH] [--resume]"; exit 0;;
    *) die 2 "Unknown option: $1";;
  esac
done
[[ -s $INPUT && -n $OUTPUT_DIR && -n $THREADS ]] || die 2 'Input FASTA, output directory and threads are required'
positive_int "$THREADS" && (( THREADS <= MAX_THREADS_PER_VIRAL_TOOL && THREADS <= MAX_TOTAL_THREADS )) || die 2 'Invalid DIAMOND thread request'
validate_diamond_performance "$THREADS" "$BLOCK_SIZE" "$INDEX_CHUNKS" "$TMPDIR"
[[ -n $DIAMOND_NR_DB && -f $DIAMOND_NR_DB ]] || die 3 'DIAMOND_NR_DB is missing or not a readable .dmnd file'
require_command diamond
require_diamond_version
OUT="$OUTPUT_DIR/02c_diamond_virus"; HITS="$OUT/nr_virus_hits.tsv"; COMPLETE="$OUT/parameters.env"
# A zero-row hit table is a valid completed result.  Its companion parameter
# file is the completion marker; a bare/empty TSV can also be left behind when
# DIAMOND is interrupted before it has written any results.
if [[ -f $HITS && -f $COMPLETE ]]; then
  (( RESUME )) && { echo "[INFO] DIAMOND virus-discovery hits already exist; skipped"; exit 0; }
  die 4 "DIAMOND virus-discovery output already exists; use --resume or a new output directory"
fi
if [[ -e $OUT ]] && ! dir_is_empty_or_missing "$OUT"; then
  (( RESUME )) || die 4 "DIAMOND virus-discovery output is incomplete or conflicting: $OUT"
  echo "[INFO] Restarting incomplete DIAMOND virus-discovery step"
fi
mkdir -p "$OUT"
# DIAMOND 2.2.4+ supplies the source-lineage fields retained for evidence
# review alongside the TaxonKit-derived classification.
fields=(qseqid qlen qstart qend pident length evalue bitscore sseqid staxids sscinames slineages)
args=(blastx --db "$DIAMOND_NR_DB" --query "$INPUT" --out "$HITS" --outfmt 6 "${fields[@]}" --threads "$THREADS" --block-size "$BLOCK_SIZE" --index-chunks "$INDEX_CHUNKS" -t "$TMPDIR" --evalue "$DIAMOND_EVALUE" --max-target-seqs "$DIAMOND_NR_MAX_TARGET_SEQS" --taxonlist "$DIAMOND_DEFAULT_TAXONLIST")
[[ $DIAMOND_SENSITIVITY == more-sensitive ]] && args+=(--more-sensitive)
printf '%q ' diamond "${args[@]}" > "$OUT/diamond_command.sh"; printf '\n' >> "$OUT/diamond_command.sh"
diamond "${args[@]}" >"$OUT/diamond.stdout.log" 2>"$OUT/diamond.stderr.log"
[[ -f $HITS ]] || die 1 "DIAMOND did not create a virus-discovery hit table; inspect: $OUT/diamond.stderr.log"
printf 'input=%s\ntaxonlist=%s\nmax_target_seqs=%s\nthreads=%s\nblock_size=%s\nindex_chunks=%s\ntmpdir=%s\n' "$INPUT" "$DIAMOND_DEFAULT_TAXONLIST" "$DIAMOND_NR_MAX_TARGET_SEQS" "$THREADS" "$BLOCK_SIZE" "$INDEX_CHUNKS" "$TMPDIR" > "$OUT/parameters.env"
echo "[INFO] DIAMOND virus-discovery hits: $HITS"
