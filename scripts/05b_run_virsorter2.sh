#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"; load_pipeline_config

INPUT='' OUTPUT_DIR='' THREADS='' RESUME=0
while (($#)); do
  case "$1" in
    --input|--output-dir|--threads)
      (($# >= 2)) || die 2 "Missing value for $1"
      case "$1" in --input) INPUT=$2;; --output-dir) OUTPUT_DIR=$2;; --threads) THREADS=$2;; esac
      shift 2;;
    --resume) RESUME=1; shift;;
    -h|--help) echo "Usage: $0 --input FASTA --output-dir PATH --threads N [--resume]"; exit 0;;
    *) die 2 "Unknown option: $1";;
  esac
done
[[ -s $INPUT && -n $OUTPUT_DIR && -n $THREADS ]] || die 2 'Input, output and threads are required'
positive_int "$THREADS" && (( THREADS <= MAX_THREADS_PER_VIRAL_TOOL )) || die 2 'Invalid VirSorter2 thread request'
[[ $VIRSORTER_USE_CONDA_OFF == 0 || $VIRSORTER_USE_CONDA_OFF == 1 ]] || die 2 'VIRSORTER_USE_CONDA_OFF must be 0 or 1 in pipeline.env'
require_executable "$VIRSORTER_COMMAND"
OUT="$OUTPUT_DIR/02b_virsorter2"; CANDIDATES="$OUT/viral_candidates_virsorter2.fna"
if [[ -f $CANDIDATES ]]; then
  (( RESUME )) && { echo "[INFO] VirSorter2 candidates already exist; skipped"; exit 0; }
  die 4 "VirSorter2 output already exists; use --resume or choose a new report output directory"
fi
if [[ -e $OUT ]] && ! dir_is_empty_or_missing "$OUT"; then die 4 "VirSorter2 output is incomplete or conflicting: $OUT"; fi
mkdir -p "$OUT"
# The installed VirSorter2 database is managed by VirSorter2 itself.  We keep
# its default groups so projects can select the database policy centrally.
VIRSORTER_ARGS=(run -i "$INPUT" -w "$OUT/run" -j "$THREADS")
(( VIRSORTER_USE_CONDA_OFF == 0 )) || VIRSORTER_ARGS+=(--use-conda-off)
VIRSORTER_ARGS+=(all)
"$VIRSORTER_COMMAND" "${VIRSORTER_ARGS[@]}" >"$OUT/virsorter2.stdout.log" 2>"$OUT/virsorter2.stderr.log"
VIRUS_FASTA=$(find "$OUT/run" -type f -name 'final-viral-combined.fa' -print -quit)
if [[ -n $VIRUS_FASTA && -s $VIRUS_FASTA ]]; then
  cp -- "$VIRUS_FASTA" "$CANDIDATES"
else
  : > "$CANDIDATES"
  echo "[WARN] VirSorter2 found no viral candidate FASTA; continuing with the other discovery methods" >&2
fi
echo "[INFO] VirSorter2 candidates: $CANDIDATES"
