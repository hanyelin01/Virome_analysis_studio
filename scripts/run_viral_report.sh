#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"; load_pipeline_config

usage() { cat <<'EOF'
Usage: run_viral_report.sh --assembly-dir PATH --cleandata-dir PATH
       --clean-layout sample_subdirs|flat --read-type pe|se --output-dir PATH
       --threads N [--groups-file PATH] [--min-contig-length N]
       [--post-checkv-min-length N]
       [--overview-rank family|genus|species] [--enable-virsorter2] [--resume]
EOF
}
ASSEMBLY_DIR='' CLEAN_DIR='' CLEAN_LAYOUT='' READ_TYPE='' OUTPUT_DIR='' THREADS='' GROUPS_FILE='' MIN_LENGTH="$VIRAL_MIN_CONTIG_LEN" POST_CHECKV_MIN_LENGTH='' OVERVIEW_RANK="$BATCH_OVERVIEW_RANK" ENABLE_VIRSORTER2=0 RESUME=0
while (($#)); do
  case "$1" in
    --assembly-dir|--cleandata-dir|--clean-layout|--read-type|--output-dir|--threads|--groups-file|--min-contig-length|--post-checkv-min-length|--overview-rank)
      (($# >= 2)) || die 2 "Missing value for $1"
      case "$1" in
        --assembly-dir) ASSEMBLY_DIR=$2;; --cleandata-dir) CLEAN_DIR=$2;; --clean-layout) CLEAN_LAYOUT=$2;; --read-type) READ_TYPE=$2;; --output-dir) OUTPUT_DIR=$2;; --threads) THREADS=$2;; --groups-file) GROUPS_FILE=$2;; --min-contig-length) MIN_LENGTH=$2;; --post-checkv-min-length) POST_CHECKV_MIN_LENGTH=$2;; --overview-rank) OVERVIEW_RANK=$2;;
      esac
      shift 2;;
    --enable-virsorter2) ENABLE_VIRSORTER2=1; shift;;
    --resume) RESUME=1; shift;;
    -h|--help) usage; exit 0;;
    *) die 2 "Unknown option: $1";;
  esac
done
[[ -n $ASSEMBLY_DIR && -n $CLEAN_DIR && -n $OUTPUT_DIR && -n $THREADS ]] || die 2 'Assembly, clean data, output and threads are required'
[[ $CLEAN_LAYOUT == sample_subdirs || $CLEAN_LAYOUT == flat ]] || die 2 'Invalid --clean-layout'
[[ $READ_TYPE == pe || $READ_TYPE == se ]] || die 2 'Invalid --read-type'
[[ $OVERVIEW_RANK == family || $OVERVIEW_RANK == genus || $OVERVIEW_RANK == species ]] || die 2 'overview rank must be family, genus or species'
POST_CHECKV_MIN_LENGTH=${POST_CHECKV_MIN_LENGTH:-$MIN_LENGTH}
positive_int "$THREADS" && positive_int "$MIN_LENGTH" && positive_int "$POST_CHECKV_MIN_LENGTH" || die 2 'Thread count and minimum contig lengths must be positive integers'
(( THREADS <= MAX_THREADS_PER_VIRAL_TOOL && THREADS <= MAX_TOTAL_THREADS )) || die 2 "Viral-tool thread request exceeds configured limits: $THREADS"
[[ -n $GENOMAD_DB ]] || die 3 'GENOMAD_DB is not configured; set it in config/pipeline.env'
[[ -n $CHECKV_DB ]] || die 3 'CHECKV_DB is not configured; set it in config/pipeline.env'
assert_existing_dir 'assembly directory' "$ASSEMBLY_DIR"; assert_existing_dir 'cleandata directory' "$CLEAN_DIR"
require_allowed_path 'assembly directory' "$ASSEMBLY_DIR"; require_allowed_path 'cleandata directory' "$CLEAN_DIR"; require_allowed_path 'viral report output directory' "$OUTPUT_DIR"
[[ -z $GROUPS_FILE || -f $GROUPS_FILE ]] || die 3 "Group metadata file does not exist: $GROUPS_FILE"
[[ -z $GROUPS_FILE ]] || require_allowed_path 'group metadata file' "$GROUPS_FILE"
# VirSorter2 is optional, but it must be available before any expensive
# preparation or geNomad step is started when the user selected it.
(( ENABLE_VIRSORTER2 == 0 )) || require_executable "$VIRSORTER_COMMAND"

mkdir -p "$OUTPUT_DIR/.contig_pipeline/runs" "$OUTPUT_DIR/reports"
RUN_ID="$(date '+%Y%m%d_%H%M%S')_$$"; RUN_DIR="$OUTPUT_DIR/.contig_pipeline/runs/$RUN_ID"; MANIFEST="$RUN_DIR/sample_manifest.tsv"; LOG="$RUN_DIR/pipeline.log"
mkdir -p "$RUN_DIR"; exec 9>"$OUTPUT_DIR/.contig_pipeline/.pipeline.lock"; flock -n 9 || die 75 "Another viral report is already running for this output location: $OUTPUT_DIR"
exec > >(tee -a "$LOG") 2>&1; printf 'RUNNING\n' > "$RUN_DIR/status"
cat > "$RUN_DIR/parameters.env" <<EOF
TASK=viral_report
ASSEMBLY_DIR=$ASSEMBLY_DIR
CLEAN_DIR=$CLEAN_DIR
CLEAN_LAYOUT=$CLEAN_LAYOUT
READ_TYPE=$READ_TYPE
OUTPUT_DIR=$OUTPUT_DIR
THREADS=$THREADS
MIN_CONTIG_LENGTH=$MIN_LENGTH
POST_CHECKV_MIN_LENGTH=$POST_CHECKV_MIN_LENGTH
ENABLE_VIRSORTER2=$ENABLE_VIRSORTER2
GROUPS_FILE=$GROUPS_FILE
OVERVIEW_RANK=$OVERVIEW_RANK
VOTU_ANI=$VOTU_ANI
VOTU_ALIGNED_FRACTION=$VOTU_ALIGNED_FRACTION
COVERM_MIN_READ_PERCENT_IDENTITY=$COVERM_MIN_READ_PERCENT_IDENTITY
COVERM_MIN_READ_ALIGNED_PERCENT=$COVERM_MIN_READ_ALIGNED_PERCENT
COVERM_MIN_COVERED_FRACTION=$COVERM_MIN_COVERED_FRACTION
VOTU_IMPORTANCE_RELATIVE_ABUNDANCE=$VOTU_IMPORTANCE_RELATIVE_ABUNDANCE
RESUME=$RESUME
EOF
run_step() { local label=$1; shift; echo "[STEP] $label"; "$@" 2>&1 | tee -a "$RUN_DIR/${label}.log"; }
on_error() { local rc=$?; printf 'FAILED\n' > "$RUN_DIR/status"; echo "[ERROR] Viral report stopped; log: $LOG"; exit "$rc"; }
trap on_error ERR

run_step preflight bash "$SCRIPT_DIR/00_preflight.sh" --task assembly_only --cleandata-dir "$CLEAN_DIR" --clean-layout "$CLEAN_LAYOUT" --read-type "$READ_TYPE" --assembly-dir "$ASSEMBLY_DIR" --manifest "$MANIFEST"
prepare=(bash "$SCRIPT_DIR/04_prepare_viral_contigs.sh" --assembly-dir "$ASSEMBLY_DIR" --manifest "$MANIFEST" --output-dir "$OUTPUT_DIR" --min-contig-length "$MIN_LENGTH"); (( RESUME )) && prepare+=(--resume); run_step prepare_contigs "${prepare[@]}"
genomad=(bash "$SCRIPT_DIR/05_run_genomad.sh" --input "$OUTPUT_DIR/01_prepared_contigs/merged_assembled_contigs.fna" --output-dir "$OUTPUT_DIR" --database "$GENOMAD_DB" --threads "$THREADS"); (( RESUME )) && genomad+=(--resume); run_step genomad "${genomad[@]}"
CANDIDATE_INPUT="$OUTPUT_DIR/02_genomad/viral_candidates_genomad.fna"
if (( ENABLE_VIRSORTER2 )); then
  virsorter=(bash "$SCRIPT_DIR/05b_run_virsorter2.sh" --input "$OUTPUT_DIR/01_prepared_contigs/merged_assembled_contigs.fna" --output-dir "$OUTPUT_DIR" --threads "$THREADS"); (( RESUME )) && virsorter+=(--resume); run_step virsorter2 "${virsorter[@]}"
  mkdir -p "$OUTPUT_DIR/02_candidates"
  MERGED="$OUTPUT_DIR/02_candidates/viral_candidates_union.fna"
  if [[ ! -s $MERGED ]]; then python3 "$SCRIPT_DIR/helpers/merge_fasta_by_id.py" --output "$MERGED" "$CANDIDATE_INPUT" "$OUTPUT_DIR/02b_virsorter2/viral_candidates_virsorter2.fna"; fi
  CANDIDATE_INPUT="$MERGED"
fi
checkv=(bash "$SCRIPT_DIR/06_run_checkv.sh" --input "$CANDIDATE_INPUT" --output-dir "$OUTPUT_DIR" --database "$CHECKV_DB" --threads "$THREADS"); (( RESUME )) && checkv+=(--resume); run_step checkv "${checkv[@]}"
abundance=(bash "$SCRIPT_DIR/07_votu_abundance.sh" --manifest "$MANIFEST" --input "$OUTPUT_DIR/03_checkv/viral_candidates_checkv.fna" --output-dir "$OUTPUT_DIR" --threads "$THREADS" --min-length "$POST_CHECKV_MIN_LENGTH"); (( RESUME )) && abundance+=(--resume); run_step votu_abundance "${abundance[@]}"
report=(bash "$SCRIPT_DIR/08_generate_viral_report.sh" --output-dir "$OUTPUT_DIR" --manifest "$MANIFEST" --overview-rank "$OVERVIEW_RANK"); [[ -z $GROUPS_FILE ]] || report+=(--groups-file "$GROUPS_FILE"); (( RESUME )) && report+=(--resume); run_step report "${report[@]}"
printf 'SUCCESS\n' > "$RUN_DIR/status"; echo "[INFO] Viral reports completed: $OUTPUT_DIR/reports/virome_dashboard.html"
