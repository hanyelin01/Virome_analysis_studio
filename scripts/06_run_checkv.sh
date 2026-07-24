#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"; load_pipeline_config

INPUT='' OUTPUT_DIR='' DATABASE="$CHECKV_DB" THREADS='' STAGE_DIR='03_checkv' RESUME=0
while (($#)); do
  case "$1" in
    --input|--output-dir|--database|--threads|--stage-dir)
      (($# >= 2)) || die 2 "Missing value for $1"
      case "$1" in --input) INPUT=$2;; --output-dir) OUTPUT_DIR=$2;; --database) DATABASE=$2;; --threads) THREADS=$2;; --stage-dir) STAGE_DIR=$2;; esac
      shift 2;;
    --resume) RESUME=1; shift;;
    -h|--help) echo "Usage: $0 --input FASTA --output-dir PATH --threads N [--database PATH] [--stage-dir NAME] [--resume]"; exit 0;;
    *) die 2 "Unknown option: $1";;
  esac
done
[[ -s $INPUT && -n $OUTPUT_DIR && -n $THREADS ]] || die 2 'Input, output and threads are required'
positive_int "$THREADS" && (( THREADS <= MAX_THREADS_PER_VIRAL_TOOL )) || die 2 'Invalid CheckV thread request'
[[ $STAGE_DIR =~ ^[A-Za-z0-9._-]+$ ]] || die 2 '--stage-dir must be a simple directory name'
require_command checkv
if [[ -n $DATABASE ]]; then [[ -d $DATABASE ]] || die 3 "CheckV database is missing or not a directory: $DATABASE"; export CHECKVDB="$DATABASE"; fi
OUT="$OUTPUT_DIR/$STAGE_DIR"; QUALITY="$OUT/quality_summary.tsv"; FILTERED="$OUT/viral_candidates_checkv.fna"
if [[ -s $QUALITY && -s $FILTERED ]]; then
  (( RESUME )) && { echo "[INFO] CheckV output already exists; skipped"; exit 0; }
  die 4 "CheckV output already exists; use --resume or choose a new report output directory"
fi
if [[ -e $OUT ]] && ! dir_is_empty_or_missing "$OUT"; then die 4 "CheckV output is incomplete or conflicting: $OUT"; fi
mkdir -p "$OUT"
checkv end_to_end "$INPUT" "$OUT/run" -t "$THREADS" >"$OUT/checkv.stdout.log" 2>"$OUT/checkv.stderr.log"
[[ -s "$OUT/run/quality_summary.tsv" ]] || die 1 "CheckV finished without quality_summary.tsv; inspect: $OUT/checkv.stderr.log"
cp -- "$OUT/run/quality_summary.tsv" "$QUALITY"
FASTA_FILES=()
[[ -s "$OUT/run/viruses.fna" ]] && FASTA_FILES+=("$OUT/run/viruses.fna")
[[ -s "$OUT/run/proviruses.fna" ]] && FASTA_FILES+=("$OUT/run/proviruses.fna")
(( ${#FASTA_FILES[@]} > 0 )) || die 1 "CheckV finished without viruses.fna/proviruses.fna"
python3 "$SCRIPT_DIR/helpers/merge_fasta_by_id.py" --output "$FILTERED" "${FASTA_FILES[@]}"
echo "[INFO] CheckV quality summary: $QUALITY"
