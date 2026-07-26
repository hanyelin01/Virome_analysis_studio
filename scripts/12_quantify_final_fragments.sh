#!/usr/bin/env bash
# Quantify only the final, sample-distributed CheckV fragments; no vOTU stage.
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"; load_pipeline_config

MANIFEST='' SAMPLE_DIR='' OUTPUT_DIR='' THREADS='' RESUME=0
while (($#)); do
  case "$1" in
    --manifest|--sample-dir|--output-dir|--threads)
      (($# >= 2)) || die 2 "Missing value for $1"
      case "$1" in --manifest) MANIFEST=$2;; --sample-dir) SAMPLE_DIR=$2;; --output-dir) OUTPUT_DIR=$2;; --threads) THREADS=$2;; esac
      shift 2;;
    --resume) RESUME=1; shift;;
    -h|--help) echo "Usage: $0 --manifest PATH --sample-dir PATH --output-dir PATH --threads N [--resume]"; exit 0;;
    *) die 2 "Unknown option: $1";;
  esac
done
[[ -f $MANIFEST && -d $SAMPLE_DIR && -n $OUTPUT_DIR && -n $THREADS ]] || die 2 'Manifest, sample directory, output directory and threads are required'
positive_int "$THREADS" && (( THREADS <= MAX_THREADS_PER_VIRAL_TOOL && THREADS <= MAX_TOTAL_THREADS )) || die 2 'Invalid abundance thread request'
require_command coverm
mkdir -p "$OUTPUT_DIR"
printf 'sample_id\tstatus\tmessage\n' > "$OUTPUT_DIR/abundance_status.tsv"
while IFS= read -r row; do
  IFS=$'\x1f' read -r sample read_type raw1 raw2 clean1 clean2 clean_single assembly <<< "${row//$'\t'/$'\x1f'}"
  root="$SAMPLE_DIR/$sample"; refs="$root/references"; annotations="$root/viral_fragments.tsv"; sample_out="$OUTPUT_DIR/$sample"; coverage="$sample_out/fragment_coverage.tsv"; counts="$sample_out/fragment_counts.tsv"; final="$root/viral_fragments_quantified.tsv"
  if [[ -s $final && $RESUME -eq 1 ]]; then printf '%s\tSUCCESS\tresumed\n' "$sample" >> "$OUTPUT_DIR/abundance_status.tsv"; continue; fi
  mkdir -p "$sample_out"
  if [[ ! -d $refs || -z $(find "$refs" -maxdepth 1 -name '*.fna' -size +0c -print -quit) ]]; then
    python3 "$SCRIPT_DIR/helpers/join_fragment_coverage.py" --annotations "$annotations" --coverage /dev/null --counts /dev/null --output "$final"
    printf '%s\tSUCCESS\tno source fragments\n' "$sample" >> "$OUTPUT_DIR/abundance_status.tsv"; continue
  fi
  # CoverM 0.8 cannot calculate `count` together with a positive covered-
  # fraction threshold.  Keep the biologically filtered coverage metrics and
  # generate raw mapped-read counts in a second, explicitly unfiltered call.
  coverage_command=(coverm genome --genome-fasta-directory "$refs" --genome-fasta-extension fna --min-read-percent-identity "$COVERM_MIN_READ_PERCENT_IDENTITY" --min-read-aligned-percent "$COVERM_MIN_READ_ALIGNED_PERCENT" --min-covered-fraction "$COVERM_MIN_COVERED_FRACTION" -m mean relative_abundance covered_bases -t "$THREADS" -o "$coverage")
  count_command=(coverm genome --genome-fasta-directory "$refs" --genome-fasta-extension fna --min-read-percent-identity "$COVERM_MIN_READ_PERCENT_IDENTITY" --min-read-aligned-percent "$COVERM_MIN_READ_ALIGNED_PERCENT" --min-covered-fraction 0 -m count -t "$THREADS" -o "$counts")
  if [[ $read_type == pe ]]; then
    coverage_command+=(--coupled "$clean1" "$clean2"); count_command+=(--coupled "$clean1" "$clean2")
  else
    coverage_command+=(--single "$clean_single"); count_command+=(--single "$clean_single")
  fi
  "${coverage_command[@]}" >"$sample_out/coverm_coverage.stdout.log" 2>"$sample_out/coverm_coverage.stderr.log"
  "${count_command[@]}" >"$sample_out/coverm_count.stdout.log" 2>"$sample_out/coverm_count.stderr.log"
  python3 "$SCRIPT_DIR/helpers/join_fragment_coverage.py" --annotations "$annotations" --coverage "$coverage" --counts "$counts" --output "$final"
  printf '%s\tSUCCESS\tquantified\n' "$sample" >> "$OUTPUT_DIR/abundance_status.tsv"
done < <(manifest_rows "$MANIFEST")
echo "[INFO] Final-fragment abundance tables: $OUTPUT_DIR"
