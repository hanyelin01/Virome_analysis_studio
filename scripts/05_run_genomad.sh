#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"; load_pipeline_config

INPUT='' OUTPUT_DIR='' DATABASE="$GENOMAD_DB" THREADS='' RESUME=0
while (($#)); do
  case "$1" in
    --input|--output-dir|--database|--threads)
      (($# >= 2)) || die 2 "Missing value for $1"
      case "$1" in --input) INPUT=$2;; --output-dir) OUTPUT_DIR=$2;; --database) DATABASE=$2;; --threads) THREADS=$2;; esac
      shift 2;;
    --resume) RESUME=1; shift;;
    -h|--help) echo "Usage: $0 --input FASTA --output-dir PATH --threads N [--database PATH] [--resume]"; exit 0;;
    *) die 2 "Unknown option: $1";;
  esac
done
[[ -s $INPUT && -n $OUTPUT_DIR && -n $DATABASE && -n $THREADS ]] || die 2 'Input, output, database and threads are required'
positive_int "$THREADS" && (( THREADS <= MAX_THREADS_PER_VIRAL_TOOL )) || die 2 'Invalid geNomad thread request'
[[ -d $DATABASE ]] || die 3 "geNomad database is missing or not a directory: $DATABASE"
require_command genomad
OUT="$OUTPUT_DIR/02_genomad"; CANDIDATES="$OUT/viral_candidates_genomad.fna"
if [[ -s $CANDIDATES ]]; then
  (( RESUME )) && { echo "[INFO] geNomad viral candidates already exist; skipped"; exit 0; }
  die 4 "geNomad output already exists; use --resume or choose a new report output directory"
fi
if [[ -e $OUT ]] && ! dir_is_empty_or_missing "$OUT"; then die 4 "geNomad output is incomplete or conflicting: $OUT"; fi
mkdir -p "$OUT"
# Splits reduce peak memory consumption on shared servers.  The thread count
# remains bounded by run_viral_report.sh and the central configuration.
genomad end-to-end --cleanup --splits "$THREADS" "$INPUT" "$OUT/run" "$DATABASE" >"$OUT/genomad.stdout.log" 2>"$OUT/genomad.stderr.log"
VIRUS_FASTA=$(find "$OUT/run" -type f -name '*_virus.fna' -print -quit)
[[ -n $VIRUS_FASTA && -s $VIRUS_FASTA ]] || die 1 "geNomad finished without a virus FASTA; inspect: $OUT/genomad.stderr.log"
cp -- "$VIRUS_FASTA" "$CANDIDATES"
SUMMARY=$(find "$OUT/run" -type f -name '*_virus_summary.tsv' -print -quit || true)
[[ -z $SUMMARY ]] || cp -- "$SUMMARY" "$OUT/virus_summary.tsv"
echo "[INFO] geNomad candidates: $CANDIDATES"
