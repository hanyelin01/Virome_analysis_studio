#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"; load_pipeline_config

ASSEMBLY_DIR='' OUTPUT_DIR='' MIN_LENGTH="$VIRAL_MIN_CONTIG_LEN" RESUME=0
while (($#)); do
  case "$1" in
    --assembly-dir|--output-dir|--min-contig-length)
      (($# >= 2)) || die 2 "Missing value for $1"
      case "$1" in --assembly-dir) ASSEMBLY_DIR=$2;; --output-dir) OUTPUT_DIR=$2;; --min-contig-length) MIN_LENGTH=$2;; esac
      shift 2;;
    --resume) RESUME=1; shift;;
    -h|--help) echo "Usage: $0 --assembly-dir PATH --output-dir PATH [--min-contig-length N] [--resume]"; exit 0;;
    *) die 2 "Unknown option: $1";;
  esac
done
[[ -n $ASSEMBLY_DIR && -n $OUTPUT_DIR ]] || die 2 '--assembly-dir and --output-dir are required'
positive_int "$MIN_LENGTH" || die 2 '--min-contig-length must be a positive integer'
assert_existing_dir 'assembly directory' "$ASSEMBLY_DIR"
require_allowed_path 'assembly directory' "$ASSEMBLY_DIR"; require_allowed_path 'viral report output directory' "$OUTPUT_DIR"
OUT="$OUTPUT_DIR/01_prepared_contigs"; FASTA="$OUT/merged_assembled_contigs.fna"
if [[ -s $FASTA ]]; then
  (( RESUME )) && { echo "[INFO] Prepared contigs already exist; skipped"; exit 0; }
  die 4 "Prepared-contig output already exists; use --resume or choose a new report output directory"
fi
if [[ -e $OUT ]] && ! dir_is_empty_or_missing "$OUT"; then die 4 "Prepared-contig output is incomplete or conflicting: $OUT"; fi
mkdir -p "$OUT"
python3 "$SCRIPT_DIR/helpers/prepare_viral_contigs.py" --assembly-dir "$ASSEMBLY_DIR" --output-dir "$OUT" --min-length "$MIN_LENGTH"
