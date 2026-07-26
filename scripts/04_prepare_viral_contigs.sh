#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"; load_pipeline_config

ASSEMBLY_DIR='' MANIFEST='' OUTPUT_DIR='' MIN_LENGTH="$VIRAL_MIN_CONTIG_LEN" RESUME=0
while (($#)); do
  case "$1" in
    --assembly-dir|--manifest|--output-dir|--min-contig-length)
      (($# >= 2)) || die 2 "Missing value for $1"
      case "$1" in --assembly-dir) ASSEMBLY_DIR=$2;; --manifest) MANIFEST=$2;; --output-dir) OUTPUT_DIR=$2;; --min-contig-length) MIN_LENGTH=$2;; esac
      shift 2;;
    --resume) RESUME=1; shift;;
    -h|--help) echo "Usage: $0 --assembly-dir PATH --manifest PATH --output-dir PATH [--min-contig-length N] [--resume]"; exit 0;;
    *) die 2 "Unknown option: $1";;
  esac
done
[[ -n $ASSEMBLY_DIR && -f $MANIFEST && -n $OUTPUT_DIR ]] || die 2 '--assembly-dir, --manifest and --output-dir are required'
positive_int "$MIN_LENGTH" || die 2 '--min-contig-length must be a positive integer'
assert_existing_dir 'assembly directory' "$ASSEMBLY_DIR"
require_allowed_path 'assembly directory' "$ASSEMBLY_DIR"; require_allowed_path 'viral report output directory' "$OUTPUT_DIR"
OUT="$OUTPUT_DIR/01_prepared_contigs"; FASTA="$OUT/merged_assembled_contigs.fna"; FINGERPRINT="$OUT/preparation_inputs.json"
if [[ -s $FASTA ]]; then
  (( RESUME )) || die 4 "Prepared-contig output already exists; use --resume or choose a new report output directory"
fi
if [[ -e $OUT ]] && ! dir_is_empty_or_missing "$OUT" && (( ! RESUME )); then die 4 "Prepared-contig output is incomplete or conflicting: $OUT"; fi
# Versions created before manifest fingerprinting cannot be resumed safely.
# Preserve them under the report's audit directory and rebuild from the current
# manifest instead of silently reusing or deleting historical contigs.
if (( RESUME )) && [[ -e $OUT ]] && ! dir_is_empty_or_missing "$OUT" && [[ ! -f $FINGERPRINT ]]; then
  ARCHIVE_ROOT="$OUTPUT_DIR/.contig_pipeline/legacy_prepared_contigs"
  ARCHIVE="$ARCHIVE_ROOT/$(date '+%Y%m%d_%H%M%S')_$$_01_prepared_contigs"
  mkdir -p "$ARCHIVE_ROOT"
  mv -- "$OUT" "$ARCHIVE"
  printf 'archived_at\tarchive_path\treason\n' > "$ARCHIVE_ROOT/latest_archive.tsv"
  printf '%s\t%s\t%s\n' "$(date -Is)" "$ARCHIVE" 'missing_manifest_fingerprint' >> "$ARCHIVE_ROOT/latest_archive.tsv"
  echo "[INFO] Archived legacy prepared contigs without a manifest fingerprint: $ARCHIVE"
  echo "[INFO] Rebuilding prepared contigs from the current manifest"
fi
mkdir -p "$OUT"
command=(python3 "$SCRIPT_DIR/helpers/prepare_viral_contigs.py" --assembly-dir "$ASSEMBLY_DIR" --manifest "$MANIFEST" --output-dir "$OUT" --min-length "$MIN_LENGTH")
(( RESUME )) && command+=(--resume)
"${command[@]}"
