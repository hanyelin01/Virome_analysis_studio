#!/usr/bin/env bash
# Independently add DIAMOND/MEGAN and/or DIAMOND/TaxonKit annotations.
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"; load_pipeline_config

usage() { cat <<'EOF'
Usage:
  run_fine_annotation.sh --source checkv --viral-report-dir PATH --mode megan|taxonomy|both --threads N
                         [--taxon-scope virus|none|custom] [--taxonlist IDs] [--max-target-seqs N]
                         [--block-size GB] [--index-chunks N] [--tmpdir PATH] [--resume]
  run_fine_annotation.sh --source custom --custom-input PATH --custom-input-type file|directory --output-dir PATH
                         --mode megan|taxonomy|both --threads N
                         [--taxon-scope virus|none|custom] [--taxonlist IDs] [--max-target-seqs N]
                         [--block-size GB] [--index-chunks N] [--tmpdir PATH] [--resume]
EOF
}

SOURCE='' VIRAL_REPORT_DIR='' CUSTOM_INPUT='' CUSTOM_INPUT_TYPE='file' OUTPUT_DIR='' MODE='both' THREADS='' TAXON_SCOPE='virus' TAXONLIST='' MAX_TARGETS="$DIAMOND_NR_MAX_TARGET_SEQS" BLOCK_SIZE="$DIAMOND_BLOCK_SIZE" INDEX_CHUNKS="$DIAMOND_INDEX_CHUNKS" TMPDIR="$DIAMOND_TMPDIR" RESUME=0
while (($#)); do
  case "$1" in
    --source|--viral-report-dir|--custom-input|--custom-input-type|--output-dir|--mode|--threads|--taxon-scope|--taxonlist|--max-target-seqs|--block-size|--index-chunks|--tmpdir)
      (($# >= 2)) || die 2 "Missing value for $1"
      case "$1" in
        --source) SOURCE=$2;; --viral-report-dir) VIRAL_REPORT_DIR=$2;; --custom-input) CUSTOM_INPUT=$2;; --custom-input-type) CUSTOM_INPUT_TYPE=$2;; --output-dir) OUTPUT_DIR=$2;; --mode) MODE=$2;; --threads) THREADS=$2;; --taxon-scope) TAXON_SCOPE=$2;; --taxonlist) TAXONLIST=$2;; --max-target-seqs) MAX_TARGETS=$2;; --block-size) BLOCK_SIZE=$2;; --index-chunks) INDEX_CHUNKS=$2;; --tmpdir) TMPDIR=$2;;
      esac
      shift 2;;
    --resume) RESUME=1; shift;;
    -h|--help) usage; exit 0;;
    *) die 2 "Unknown option: $1";;
  esac
done
[[ $SOURCE == checkv || $SOURCE == custom ]] || die 2 '--source must be checkv or custom'
[[ $MODE == megan || $MODE == taxonomy || $MODE == both ]] || die 2 '--mode must be megan, taxonomy or both'
[[ $TAXON_SCOPE == virus || $TAXON_SCOPE == none || $TAXON_SCOPE == custom ]] || die 2 'Invalid --taxon-scope'
positive_int "$THREADS" && positive_int "$MAX_TARGETS" || die 2 'Threads and max-target-seqs must be positive integers'
(( THREADS <= MAX_THREADS_PER_VIRAL_TOOL && THREADS <= MAX_TOTAL_THREADS )) || die 2 'Fine-annotation thread request exceeds configured limits'
validate_diamond_performance "$THREADS" "$BLOCK_SIZE" "$INDEX_CHUNKS" "$TMPDIR"
case "$TAXON_SCOPE" in
  virus) TAXONLIST="${TAXONLIST:-$DIAMOND_DEFAULT_TAXONLIST}";;
  none) TAXONLIST='';;
  custom) [[ -n $TAXONLIST ]] || die 2 'Custom taxon scope requires --taxonlist';;
esac
[[ -z $TAXONLIST ]] || valid_taxonlist "$TAXONLIST" || die 2 'Taxon list must contain positive NCBI TaxIDs separated by commas'

CORE_MANIFEST='' CORE_GROUPS_FILE='' CORE_OVERVIEW_RANK=''
if [[ $SOURCE == checkv ]]; then
  [[ -n $VIRAL_REPORT_DIR ]] || die 2 '--viral-report-dir is required for CheckV source'
  assert_existing_dir 'viral report directory' "$VIRAL_REPORT_DIR"; require_allowed_path 'viral report directory' "$VIRAL_REPORT_DIR"
  OUTPUT_DIR="$VIRAL_REPORT_DIR"; INPUT="$VIRAL_REPORT_DIR/03_checkv/viral_candidates_checkv.fna"
  [[ -s $INPUT ]] || die 3 "CheckV candidate FASTA is missing or empty: $INPUT"
  CORE_MANIFEST=$(find "$VIRAL_REPORT_DIR/.contig_pipeline/runs" -mindepth 2 -maxdepth 2 -name sample_manifest.tsv -type f -printf '%T@\t%p\n' 2>/dev/null | sort -n | tail -n 1 | cut -f2- || true)
  if [[ -n $CORE_MANIFEST && -f "$(dirname "$CORE_MANIFEST")/parameters.env" ]]; then
    CORE_GROUPS_FILE=$(sed -n 's/^GROUPS_FILE=//p' "$(dirname "$CORE_MANIFEST")/parameters.env" | head -n 1)
    [[ -z $CORE_GROUPS_FILE || -f $CORE_GROUPS_FILE ]] || CORE_GROUPS_FILE=''
    CORE_OVERVIEW_RANK=$(sed -n 's/^OVERVIEW_RANK=//p' "$(dirname "$CORE_MANIFEST")/parameters.env" | head -n 1)
    [[ $CORE_OVERVIEW_RANK == family || $CORE_OVERVIEW_RANK == genus || $CORE_OVERVIEW_RANK == species ]] || CORE_OVERVIEW_RANK=''
  fi
else
  [[ -n $CUSTOM_INPUT && -n $OUTPUT_DIR ]] || die 2 '--custom-input and --output-dir are required for custom source'
  [[ $CUSTOM_INPUT_TYPE == file || $CUSTOM_INPUT_TYPE == directory ]] || die 2 '--custom-input-type must be file or directory'
  [[ $CUSTOM_INPUT_TYPE == file && -f $CUSTOM_INPUT || $CUSTOM_INPUT_TYPE == directory && -d $CUSTOM_INPUT ]] || die 3 "Custom input is unavailable: $CUSTOM_INPUT"
  require_allowed_path 'custom input' "$CUSTOM_INPUT"; require_allowed_path 'custom output directory' "$OUTPUT_DIR"; assert_output_parent 'custom output directory' "$OUTPUT_DIR"
  [[ ! -e $OUTPUT_DIR ]] || (( RESUME )) || die 4 "Custom output already exists; use --resume or choose a new output directory"
fi

mkdir -p "$OUTPUT_DIR/.contig_pipeline/runs" "$OUTPUT_DIR/reports"
RUN_ID="$(date '+%Y%m%d_%H%M%S')_$$"; RUN_DIR="$OUTPUT_DIR/.contig_pipeline/runs/$RUN_ID"; LOG="$RUN_DIR/pipeline.log"
TASK_REGISTRY_ID=${CONTIG_PIPELINE_TASK_ID:-}
mkdir -p "$RUN_DIR"; exec 9>"$OUTPUT_DIR/.contig_pipeline/.pipeline.lock"; flock -n 9 || die 75 "Another pipeline task is already running for: $OUTPUT_DIR"
exec > >(tee -a "$LOG") 2>&1; printf 'RUNNING\n' > "$RUN_DIR/status"
on_error() { local rc=$?; printf 'FAILED\n' > "$RUN_DIR/status"; echo "[ERROR] Fine annotation stopped; log: $LOG"; exit "$rc"; }
trap on_error ERR
run_step() { local label=$1; shift; echo "[STEP] $label"; "$@" 2>&1 | tee -a "$RUN_DIR/${label}.log"; }

if [[ $SOURCE == custom ]]; then
  INPUT="$OUTPUT_DIR/00_input/custom_candidates.fna"; INPUT_MANIFEST="$OUTPUT_DIR/00_input/input_manifest.tsv"
  if [[ -s $INPUT ]]; then
    (( RESUME )) || die 4 "Custom input preparation already exists: $INPUT"
  else
    run_step prepare_custom_input python3 "$SCRIPT_DIR/helpers/prepare_annotation_input.py" --input "$CUSTOM_INPUT" --input-type "$CUSTOM_INPUT_TYPE" --output-fasta "$INPUT" --manifest "$INPUT_MANIFEST"
  fi
else
  INPUT_MANIFEST="$RUN_DIR/checkv_input_manifest.tsv"
  printf 'source_file\tsha256\n%s\t%s\n' "$INPUT" "$(sha256sum "$INPUT" | awk '{print $1}')" > "$INPUT_MANIFEST"
fi

cat > "$RUN_DIR/parameters.env" <<EOF
TASK=fine_annotation
TASK_REGISTRY_ID=$TASK_REGISTRY_ID
SOURCE=$SOURCE
INPUT=$INPUT
OUTPUT_DIR=$OUTPUT_DIR
MODE=$MODE
THREADS=$THREADS
TAXON_SCOPE=$TAXON_SCOPE
TAXONLIST=${TAXONLIST:-all_nr}
MAX_TARGET_SEQS=$MAX_TARGETS
BLOCK_SIZE=$BLOCK_SIZE
INDEX_CHUNKS=$INDEX_CHUNKS
TMPDIR=$TMPDIR
RESUME=$RESUME
EOF
common=(--input "$INPUT" --output-dir "$OUTPUT_DIR" --threads "$THREADS" --taxon-scope "$TAXON_SCOPE" --max-target-seqs "$MAX_TARGETS" --block-size "$BLOCK_SIZE" --index-chunks "$INDEX_CHUNKS" --tmpdir "$TMPDIR")
[[ -z $TAXONLIST ]] || common+=(--taxonlist "$TAXONLIST")
if [[ $MODE == megan || $MODE == both ]]; then
  step=(bash "$SCRIPT_DIR/09_run_diamond_megan.sh" "${common[@]}"); (( RESUME )) && step+=(--resume); run_step diamond_megan "${step[@]}"
fi
if [[ $MODE == taxonomy || $MODE == both ]]; then
  step=(bash "$SCRIPT_DIR/10_run_diamond_taxonomy.sh" "${common[@]}"); (( RESUME )) && step+=(--resume); run_step diamond_taxonomy "${step[@]}"
fi
if [[ $SOURCE == checkv && -n $CORE_MANIFEST && -f $CORE_MANIFEST ]]; then
  report_step=(bash "$SCRIPT_DIR/08_generate_viral_report.sh" --output-dir "$OUTPUT_DIR" --manifest "$CORE_MANIFEST" --refresh)
  [[ -z $CORE_GROUPS_FILE ]] || report_step+=(--groups-file "$CORE_GROUPS_FILE")
  [[ -z $CORE_OVERVIEW_RANK ]] || report_step+=(--overview-rank "$CORE_OVERVIEW_RANK")
  run_step refresh_main_report "${report_step[@]}"
else
  description="custom input: $CUSTOM_INPUT"
  run_step custom_annotation_report python3 "$SCRIPT_DIR/helpers/build_fine_annotation_report.py" --output-dir "$OUTPUT_DIR" --input-description "$description" --taxon-scope "$TAXON_SCOPE" --taxonlist "${TAXONLIST:-all NR}"
fi
printf 'SUCCESS\n' > "$RUN_DIR/status"
echo "[INFO] Fine annotation completed: $OUTPUT_DIR"
