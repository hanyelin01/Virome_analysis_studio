#!/usr/bin/env bash
# Build a version-pinned DIAMOND database from reviewed ICTV-VMR-derived data.
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"; load_pipeline_config

METADATA='' PROTEINS='' OUTPUT_DIR='' VERSION=''
while (($#)); do
  case "$1" in
    --metadata|--protein-fasta|--output-dir|--version)
      (($# >= 2)) || die 2 "Missing value for $1"
      case "$1" in --metadata) METADATA=$2;; --protein-fasta) PROTEINS=$2;; --output-dir) OUTPUT_DIR=$2;; --version) VERSION=$2;; esac
      shift 2;;
    -h|--help) echo "Usage: $0 --metadata ICTV.tsv --protein-fasta ICTV.faa --output-dir PATH --version VMR_RELEASE"; exit 0;;
    *) die 2 "Unknown option: $1";;
  esac
done
[[ -f $METADATA && -s $PROTEINS && -n $OUTPUT_DIR && -n $VERSION ]] || die 2 'Metadata, protein FASTA, output directory and version are required'
if [[ -e $OUTPUT_DIR ]] && ! dir_is_empty_or_missing "$OUTPUT_DIR"; then
  die 4 "ICTV reference output already exists: $OUTPUT_DIR"
fi
require_command diamond
mkdir -p "$OUTPUT_DIR"
python3 "$SCRIPT_DIR/helpers/prepare_ictv_reference_metadata.py" --metadata "$METADATA" --protein-fasta "$PROTEINS" --output-dir "$OUTPUT_DIR" --version "$VERSION"
diamond makedb --in "$PROTEINS" --db "$OUTPUT_DIR/ictv_reference" >"$OUTPUT_DIR/diamond_makedb.stdout.log" 2>"$OUTPUT_DIR/diamond_makedb.stderr.log"
[[ -s "$OUTPUT_DIR/ictv_reference.dmnd" ]] || die 1 "DIAMOND did not create ictv_reference.dmnd"
printf 'ICTV_REFERENCE_DMND=%s\nICTV_REFERENCE_METADATA=%s\nICTV_REFERENCE_VERSION=%s\n' "$OUTPUT_DIR/ictv_reference.dmnd" "$OUTPUT_DIR/ictv_reference_metadata.tsv" "$VERSION"
