#!/usr/bin/env bash
# Build independent, sample-local vOTU catalogues and abundance tables.
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"; load_pipeline_config

MANIFEST='' INPUT='' OUTPUT_DIR='' THREADS='' MIN_LENGTH="$VOTU_POST_CHECKV_MIN_LEN" RESUME=0
while (($#)); do
  case "$1" in
    --manifest|--input|--output-dir|--threads|--min-length)
      (($# >= 2)) || die 2 "Missing value for $1"
      case "$1" in --manifest) MANIFEST=$2;; --input) INPUT=$2;; --output-dir) OUTPUT_DIR=$2;; --threads) THREADS=$2;; --min-length) MIN_LENGTH=$2;; esac
      shift 2;;
    --resume) RESUME=1; shift;;
    -h|--help) echo "Usage: $0 --manifest PATH --input CHECKV_FASTA --output-dir PATH --threads N [--min-length N] [--resume]"; exit 0;;
    *) die 2 "Unknown option: $1";;
  esac
done
[[ -f $MANIFEST && -s $INPUT && -n $OUTPUT_DIR && -n $THREADS ]] || die 2 'Manifest, input, output and threads are required'
positive_int "$THREADS" && (( THREADS <= MAX_THREADS_PER_VIRAL_TOOL )) || die 2 'Invalid CoverM thread request'
positive_int "$MIN_LENGTH" || die 2 '--min-length must be a positive integer'
require_command vclust; require_command coverm; require_command minimap2; require_command samtools

OUT="$OUTPUT_DIR/04_sample_votu"
SPLIT_SUMMARY="$OUT/split_summary.tsv"
SPLIT_COMPLETE="$OUT/split_complete.json"
STATUS="$OUT/sample_status.tsv"
CONTRACT="$OUT/run_contract.json"
PYTHON_EXEC="${CONTIG_PIPELINE_PYTHON:-$PIPELINE_HOME/.venv/bin/python}"
[[ -x $PYTHON_EXEC ]] || die 127 "Pipeline Python interpreter is unavailable: $PYTHON_EXEC"
if [[ -e $OUT ]] && ! dir_is_empty_or_missing "$OUT" && (( ! RESUME )); then
  die 4 "Sample-local vOTU output conflicts with an existing directory: $OUT"
fi
mkdir -p "$OUT"

FORCE_SPLIT=0
if [[ ! -s $CONTRACT ]]; then
  if find "$OUT" -name votu_summary.tsv -type f -size +0c -print -quit | grep -q .; then
    die 4 "Existing vOTU results predate parameter contracts and cannot be safely resumed; use a new report output directory"
  fi
  FORCE_SPLIT=1
fi
"$PYTHON_EXEC" "$SCRIPT_DIR/helpers/votu_run_contract.py" \
  --contract "$CONTRACT" --manifest "$MANIFEST" --input "$INPUT" \
  --min-length "$MIN_LENGTH" --threads "$THREADS" --ani "$VOTU_ANI" \
  --aligned-fraction "$VOTU_ALIGNED_FRACTION" \
  --read-identity "$COVERM_MIN_READ_PERCENT_IDENTITY" \
  --read-aligned-percent "$COVERM_MIN_READ_ALIGNED_PERCENT" \
  --covered-fraction "$COVERM_MIN_COVERED_FRACTION" \
  --importance-abundance "$VOTU_IMPORTANCE_RELATIVE_ABUNDANCE"

if [[ ! -s $SPLIT_COMPLETE || $FORCE_SPLIT -eq 1 ]]; then
  [[ ! -e $SPLIT_SUMMARY ]] || echo "[INFO] Previous candidate split did not complete; rebuilding it"
  "$PYTHON_EXEC" "$SCRIPT_DIR/helpers/split_viral_candidates_by_sample.py" \
    --input "$INPUT" \
    --provenance "$OUTPUT_DIR/01_prepared_contigs/contig_provenance.tsv" \
    --quality "$OUTPUT_DIR/03_checkv/quality_summary.tsv" \
    --taxonomy "$OUTPUT_DIR/02_genomad/virus_summary.tsv" \
    --manifest "$MANIFEST" \
    --output-root "$OUT" \
    --min-length "$MIN_LENGTH"
else
  echo "[INFO] Completed sample candidate split already exists; skipped"
fi

printf 'sample_id\tstatus\tmessage\n' > "$STATUS"

run_one_sample() {
  local sample=$1 read_type=$2 clean1=$3 clean2=$4 clean_single=$5
  local sample_root="$OUT/$sample" candidate_dir="$OUT/$sample/01_candidates"
  local candidates="$candidate_dir/viral_candidates_checkv.fna" metadata="$candidate_dir/candidate_metadata.tsv"
  local work="$sample_root/02_vclust" catalogue="$sample_root/03_votu" abundance="$sample_root/04_abundance"
  local summary="$sample_root/votu_summary.tsv" count=0
  if [[ -s $summary && $RESUME -eq 1 ]]; then
    echo "[INFO] $sample: local vOTU summary already exists; skipped"
    return 0
  fi
  mkdir -p "$work" "$catalogue" "$abundance"
  [[ -f $candidates ]] && count=$(grep -c '^>' "$candidates" || true)
  if (( count == 0 )); then
    "$PYTHON_EXEC" "$SCRIPT_DIR/helpers/build_sample_votu_summary.py" \
      --sample-id "$sample" --metadata "$metadata" --members "$catalogue/votu_cluster_members.tsv" \
      --representatives "$catalogue/votu_representative_map.tsv" --output "$summary" \
      --importance-abundance "$VOTU_IMPORTANCE_RELATIVE_ABUNDANCE" --allow-empty
    echo "[INFO] $sample: no CheckV candidates after the post-CheckV length filter"
    return 0
  fi
  echo "[INFO] $sample: clustering $count candidate sequence(s) locally"
  if (( count == 1 )); then
    "$PYTHON_EXEC" "$SCRIPT_DIR/helpers/build_local_votu_catalogue.py" --input "$candidates" --output-dir "$catalogue" --singletons
  else
    vclust prefilter -i "$candidates" -o "$work/prefilter.tsv" --min-ident "$(awk "BEGIN {print $VOTU_ANI / 100}")" >"$work/prefilter.stdout.log" 2>"$work/prefilter.stderr.log"
    vclust align -i "$candidates" -o "$work/ani.tsv" --filter "$work/prefilter.tsv" --outfmt lite --out-ani "$(awk "BEGIN {print $VOTU_ANI / 100}")" --out-qcov "$(awk "BEGIN {print $VOTU_ALIGNED_FRACTION / 100}")" >"$work/align.stdout.log" 2>"$work/align.stderr.log"
    vclust cluster -i "$work/ani.tsv" -o "$work/vclust_members.tsv" --ids "$work/ani.ids.tsv" --algorithm cd-hit --metric ani --ani "$(awk "BEGIN {print $VOTU_ANI / 100}")" --qcov "$(awk "BEGIN {print $VOTU_ALIGNED_FRACTION / 100}")" --out-repr >"$work/cluster.stdout.log" 2>"$work/cluster.stderr.log"
    "$PYTHON_EXEC" "$SCRIPT_DIR/helpers/build_local_votu_catalogue.py" --input "$candidates" --clusters "$work/vclust_members.tsv" --output-dir "$catalogue"
  fi
  local reps="$catalogue/representatives" coverage="$abundance/votu_coverage.tsv"
  # `count` is the number of mapped reads for each local vOTU representative.
  # It is retained separately from relative abundance, which is length-aware.
  local -a command=(coverm genome --genome-fasta-directory "$reps" --genome-fasta-extension fna --min-read-percent-identity "$COVERM_MIN_READ_PERCENT_IDENTITY" --min-read-aligned-percent "$COVERM_MIN_READ_ALIGNED_PERCENT" --min-covered-fraction "$COVERM_MIN_COVERED_FRACTION" -m mean relative_abundance covered_bases count -t "$THREADS" -o "$coverage")
  if [[ $read_type == pe ]]; then
    [[ -s $clean1 && -s $clean2 ]] || die 3 "$sample: clean PE files are missing"
    command+=(--coupled "$clean1" "$clean2")
  else
    [[ -s $clean_single ]] || die 3 "$sample: clean SE file is missing"
    command+=(--single "$clean_single")
  fi
  if ! "${command[@]}" >"$abundance/coverm.stdout.log" 2>"$abundance/coverm.stderr.log"; then
    echo "[ERROR] $sample: CoverM failed; last 80 error lines follow" >&2
    tail -n 80 "$abundance/coverm.stderr.log" >&2 || true
    return 1
  fi
  "$PYTHON_EXEC" "$SCRIPT_DIR/helpers/build_sample_votu_summary.py" \
    --sample-id "$sample" --metadata "$metadata" --members "$catalogue/votu_cluster_members.tsv" \
    --representatives "$catalogue/votu_representative_map.tsv" --coverage "$coverage" --output "$summary" \
    --importance-abundance "$VOTU_IMPORTANCE_RELATIVE_ABUNDANCE"
  echo "[INFO] $sample: local viral report data are ready"
}

failed=0
while IFS= read -r row; do
  IFS=$'\x1f' read -r sample read_type raw1 raw2 clean1 clean2 clean_single assembly <<< "${row//$'\t'/$'\x1f'}"
  if (run_one_sample "$sample" "$read_type" "$clean1" "$clean2" "$clean_single"); then
    printf '%s\tSUCCESS\tcompleted or resumed\n' "$sample" >> "$STATUS"
  else
    printf '%s\tFAILED\tinspect sample-specific logs under 04_sample_votu/%s\n' "$sample" "$sample" >> "$STATUS"
    echo "[ERROR] $sample failed; continuing with other samples" >&2
    failed=$((failed + 1))
  fi
done < <(manifest_rows "$MANIFEST")

if (( failed )); then
  echo "[WARN] $failed sample(s) failed in the sample-local vOTU stage; completed samples will still receive reports"
else
  echo "[INFO] All sample-local vOTU and abundance steps completed"
fi
