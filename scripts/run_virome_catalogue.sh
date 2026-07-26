#!/usr/bin/env bash
# Global, provenance-preserving virome catalogue workflow (v2).
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"; load_pipeline_config

usage() { cat <<'EOF'
Usage: run_virome_catalogue.sh --assembly-dir PATH --cleandata-dir PATH
       --clean-layout sample_subdirs|flat --read-type pe|se --output-dir PATH
       --threads N [--diamond-threads N] [--diamond-block-size GB]
       [--diamond-index-chunks N] [--diamond-tmpdir PATH]
       [--groups-file PATH] [--min-contig-length N] [--resume]
EOF
}
ASSEMBLY_DIR='' CLEAN_DIR='' CLEAN_LAYOUT='' READ_TYPE='' OUTPUT_DIR='' THREADS='' DIAMOND_THREADS="$DIAMOND_THREADS_PER_JOB" DIAMOND_BLOCK_SIZE_RUN="$DIAMOND_BLOCK_SIZE" DIAMOND_INDEX_CHUNKS_RUN="$DIAMOND_INDEX_CHUNKS" DIAMOND_TMPDIR_RUN="$DIAMOND_TMPDIR" GROUPS_FILE='' MIN_LENGTH="$VIRAL_MIN_CONTIG_LEN" RESUME=0
while (($#)); do
  case "$1" in
    --assembly-dir|--cleandata-dir|--clean-layout|--read-type|--output-dir|--threads|--diamond-threads|--diamond-block-size|--diamond-index-chunks|--diamond-tmpdir|--groups-file|--min-contig-length)
      (($# >= 2)) || die 2 "Missing value for $1"
      case "$1" in
        --assembly-dir) ASSEMBLY_DIR=$2;; --cleandata-dir) CLEAN_DIR=$2;; --clean-layout) CLEAN_LAYOUT=$2;; --read-type) READ_TYPE=$2;; --output-dir) OUTPUT_DIR=$2;; --threads) THREADS=$2;; --diamond-threads) DIAMOND_THREADS=$2;; --diamond-block-size) DIAMOND_BLOCK_SIZE_RUN=$2;; --diamond-index-chunks) DIAMOND_INDEX_CHUNKS_RUN=$2;; --diamond-tmpdir) DIAMOND_TMPDIR_RUN=$2;; --groups-file) GROUPS_FILE=$2;; --min-contig-length) MIN_LENGTH=$2;;
      esac
      shift 2;;
    --resume) RESUME=1; shift;;
    -h|--help) usage; exit 0;;
    *) die 2 "Unknown option: $1";;
  esac
done
[[ -n $ASSEMBLY_DIR && -n $CLEAN_DIR && -n $OUTPUT_DIR && -n $THREADS ]] || die 2 'Assembly, clean data, output and threads are required'
[[ $CLEAN_LAYOUT == sample_subdirs || $CLEAN_LAYOUT == flat ]] || die 2 'Invalid --clean-layout'
[[ $READ_TYPE == pe || $READ_TYPE == se ]] || die 2 'Invalid --read-type'
positive_int "$THREADS" && positive_int "$MIN_LENGTH" || die 2 'Thread count and minimum contig length must be positive integers'
(( THREADS <= MAX_THREADS_PER_VIRAL_TOOL && THREADS <= MAX_TOTAL_THREADS )) || die 2 'Virome-catalogue thread request exceeds configured limits'
validate_diamond_performance "$DIAMOND_THREADS" "$DIAMOND_BLOCK_SIZE_RUN" "$DIAMOND_INDEX_CHUNKS_RUN" "$DIAMOND_TMPDIR_RUN"
(( DIAMOND_THREADS <= MAX_THREADS_PER_VIRAL_TOOL && DIAMOND_THREADS <= MAX_TOTAL_THREADS )) || die 2 'DIAMOND thread request exceeds configured limits'
assert_existing_dir 'assembly directory' "$ASSEMBLY_DIR"; assert_existing_dir 'cleandata directory' "$CLEAN_DIR"
require_allowed_path 'assembly directory' "$ASSEMBLY_DIR"; require_allowed_path 'cleandata directory' "$CLEAN_DIR"; require_allowed_path 'virome output directory' "$OUTPUT_DIR"
[[ -z $GROUPS_FILE || -f $GROUPS_FILE ]] || die 3 "Group metadata file does not exist: $GROUPS_FILE"
[[ -n $GENOMAD_DB && -d $GENOMAD_DB ]] || die 3 'GENOMAD_DB is missing or not a directory'
[[ -n $CHECKV_DB && -d $CHECKV_DB ]] || die 3 'CHECKV_DB is missing or not a directory'
[[ -n $DIAMOND_NR_DB && -f $DIAMOND_NR_DB ]] || die 3 'DIAMOND_NR_DB is missing or not a readable .dmnd file'
[[ -n $TAXONKIT_DB && -d $TAXONKIT_DB ]] || die 3 'TAXONKIT_DB is missing or not a directory'
[[ -n $MEGAN_DAA2RMA && -n $MEGAN_MAP_DB ]] || die 3 'MEGAN_DAA2RMA and MEGAN_MAP_DB are required for the v2 DAA/RMA6 deliverable'
[[ -n $ICTV_REFERENCE_DMND && -f $ICTV_REFERENCE_DMND ]] || die 3 'ICTV_REFERENCE_DMND is missing or not a readable .dmnd file'
[[ -n $ICTV_REFERENCE_METADATA && -f $ICTV_REFERENCE_METADATA ]] || die 3 'ICTV_REFERENCE_METADATA is missing or not readable'
require_command genomad; require_command checkv; require_command diamond; require_command taxonkit; require_command coverm; require_executable "$VIRSORTER_COMMAND"

mkdir -p "$OUTPUT_DIR/.contig_pipeline/runs" "$OUTPUT_DIR/reports"
RUN_ID="$(date '+%Y%m%d_%H%M%S')_$$"; RUN_DIR="$OUTPUT_DIR/.contig_pipeline/runs/$RUN_ID"; MANIFEST="$RUN_DIR/sample_manifest.tsv"; LOG="$RUN_DIR/pipeline.log"
mkdir -p "$RUN_DIR"; exec 9>"$OUTPUT_DIR/.contig_pipeline/.pipeline.lock"; flock -n 9 || die 75 "Another virome catalogue is already running for this output location: $OUTPUT_DIR"
exec > >(tee -a "$LOG") 2>&1; printf 'RUNNING\n' > "$RUN_DIR/status"
cat > "$RUN_DIR/parameters.env" <<EOF
TASK=virome_catalogue_v2
ASSEMBLY_DIR=$ASSEMBLY_DIR
CLEAN_DIR=$CLEAN_DIR
CLEAN_LAYOUT=$CLEAN_LAYOUT
READ_TYPE=$READ_TYPE
OUTPUT_DIR=$OUTPUT_DIR
THREADS=$THREADS
MIN_CONTIG_LENGTH=$MIN_LENGTH
GROUPS_FILE=$GROUPS_FILE
GENOMAD_DB=$GENOMAD_DB
CHECKV_DB=$CHECKV_DB
DIAMOND_NR_DB=$DIAMOND_NR_DB
DIAMOND_DEFAULT_TAXONLIST=$DIAMOND_DEFAULT_TAXONLIST
DIAMOND_EVALUE=$DIAMOND_EVALUE
DIAMOND_NR_MAX_TARGET_SEQS=$DIAMOND_NR_MAX_TARGET_SEQS
DIAMOND_THREADS=$DIAMOND_THREADS
DIAMOND_BLOCK_SIZE=$DIAMOND_BLOCK_SIZE_RUN
DIAMOND_INDEX_CHUNKS=$DIAMOND_INDEX_CHUNKS_RUN
DIAMOND_TMPDIR=$DIAMOND_TMPDIR_RUN
ICTV_REFERENCE_DMND=$ICTV_REFERENCE_DMND
ICTV_REFERENCE_METADATA=$ICTV_REFERENCE_METADATA
ICTV_REFERENCE_VERSION=$ICTV_REFERENCE_VERSION
RESUME=$RESUME
EOF
CONTRACT="$OUTPUT_DIR/.contig_pipeline/virome_catalogue_contract.env"
if [[ -f $CONTRACT ]]; then
  if ! diff -u <(grep -Ev '^(THREADS|RESUME|DIAMOND_THREADS|DIAMOND_BLOCK_SIZE|DIAMOND_INDEX_CHUNKS|DIAMOND_TMPDIR)=' "$CONTRACT") <(grep -Ev '^(THREADS|RESUME|DIAMOND_THREADS|DIAMOND_BLOCK_SIZE|DIAMOND_INDEX_CHUNKS|DIAMOND_TMPDIR)=' "$RUN_DIR/parameters.env") >/dev/null; then
    die 4 'Existing v2 output was created with different inputs, databases, or result parameters; use a new output directory'
  fi
else
  cp -- "$RUN_DIR/parameters.env" "$CONTRACT"
fi
run_step() { local label=$1; shift; echo "[STEP] $label"; "$@" 2>&1 | tee -a "$RUN_DIR/${label}.log"; }
on_error() { local rc=$?; printf 'FAILED\n' > "$RUN_DIR/status"; echo "[ERROR] Virome catalogue stopped; log: $LOG"; exit "$rc"; }
trap on_error ERR

MEGAN_DEFERRED=0
launch_megan_background() {
  local after_foreground=${1:-0} stage="$OUTPUT_DIR/04_nr_megan" pid_file="$stage/background.pid" status_file="$stage/background.status" log_file="$stage/background.log"
  local auxiliary_threads=$DIAMOND_THREADS existing_pid
  mkdir -p "$stage"
  if [[ -s "$stage/viral_candidates.nr.daa" && -s "$stage/viral_candidates.nr.rma6" ]]; then
    printf 'state=completed\nfinished_at=%s\n' "$(date -Is)" > "$status_file"
    echo "[INFO] MEGAN auxiliary files already exist; skipped"
    return 0
  fi
  if [[ -s $pid_file ]]; then
    existing_pid=$(<"$pid_file")
    if [[ $existing_pid =~ ^[1-9][0-9]*$ ]] && kill -0 "$existing_pid" 2>/dev/null; then
      echo "[INFO] MEGAN auxiliary task is already running in background (PID $existing_pid)"
      return 0
    fi
    rm -f "$pid_file"
  fi
  if (( ! after_foreground && DIAMOND_THREADS + auxiliary_threads > MAX_TOTAL_THREADS )); then
    auxiliary_threads=$((MAX_TOTAL_THREADS - DIAMOND_THREADS))
  fi
  if (( auxiliary_threads < 1 )); then
    MEGAN_DEFERRED=1
    echo "[INFO] MEGAN auxiliary task deferred until the foreground DIAMOND task has released thread capacity"
    return 0
  fi
  require_command setsid
  local -a step=(bash "$SCRIPT_DIR/09_run_diamond_megan.sh" --input "$OUTPUT_DIR/03_candidate_catalogue/VC_catalogue.fna" --output-dir "$OUTPUT_DIR" --stage-dir 04_nr_megan --threads "$auxiliary_threads" --taxon-scope none --max-target-seqs "$DIAMOND_NR_MAX_TARGET_SEQS" --block-size "$DIAMOND_BLOCK_SIZE_RUN" --index-chunks "$DIAMOND_INDEX_CHUNKS_RUN" --tmpdir "$DIAMOND_TMPDIR_RUN" --resume)
  export AUXILIARY_STATUS_FILE="$status_file" AUXILIARY_PID_FILE="$pid_file"
  setsid bash -c '
    set +e
    # Do not retain the parent workflow lock after the main analysis exits.
    exec 9>&-
    "$@"
    rc=$?
    if (( rc == 0 )); then
      printf "state=completed\\nfinished_at=%s\\n" "$(date -Is)" > "$AUXILIARY_STATUS_FILE"
    else
      printf "state=failed\\nexit_code=%s\\nfinished_at=%s\\n" "$rc" "$(date -Is)" > "$AUXILIARY_STATUS_FILE"
    fi
    rm -f "$AUXILIARY_PID_FILE"
    exit "$rc"
  ' _ "${step[@]}" >>"$log_file" 2>&1 &
  local auxiliary_pid=$!
  printf '%s\n' "$auxiliary_pid" > "$pid_file"
  printf 'state=running\npid=%s\nstarted_at=%s\nthreads=%s\n' "$auxiliary_pid" "$(date -Is)" "$auxiliary_threads" > "$status_file"
  echo "[INFO] MEGAN DAA/RMA6 auxiliary task started in background (PID $auxiliary_pid; threads $auxiliary_threads)"
}

run_step preflight bash "$SCRIPT_DIR/00_preflight.sh" --task assembly_only --cleandata-dir "$CLEAN_DIR" --clean-layout "$CLEAN_LAYOUT" --read-type "$READ_TYPE" --assembly-dir "$ASSEMBLY_DIR" --manifest "$MANIFEST"
step=(bash "$SCRIPT_DIR/04_prepare_viral_contigs.sh" --assembly-dir "$ASSEMBLY_DIR" --manifest "$MANIFEST" --output-dir "$OUTPUT_DIR" --min-contig-length "$MIN_LENGTH"); (( RESUME )) && step+=(--resume); run_step prepare_contigs "${step[@]}"
step=(bash "$SCRIPT_DIR/05_run_genomad.sh" --input "$OUTPUT_DIR/01_prepared_contigs/merged_assembled_contigs.fna" --output-dir "$OUTPUT_DIR" --database "$GENOMAD_DB" --threads "$THREADS"); (( RESUME )) && step+=(--resume); run_step genomad "${step[@]}"
step=(bash "$SCRIPT_DIR/05b_run_virsorter2.sh" --input "$OUTPUT_DIR/01_prepared_contigs/merged_assembled_contigs.fna" --output-dir "$OUTPUT_DIR" --threads "$THREADS"); (( RESUME )) && step+=(--resume); run_step virsorter2 "${step[@]}"
step=(bash "$SCRIPT_DIR/05c_run_diamond_virus_discovery.sh" --input "$OUTPUT_DIR/01_prepared_contigs/merged_assembled_contigs.fna" --output-dir "$OUTPUT_DIR" --threads "$DIAMOND_THREADS" --block-size "$DIAMOND_BLOCK_SIZE_RUN" --index-chunks "$DIAMOND_INDEX_CHUNKS_RUN" --tmpdir "$DIAMOND_TMPDIR_RUN"); (( RESUME )) && step+=(--resume); run_step diamond_virus_discovery "${step[@]}"
if [[ ! -s "$OUTPUT_DIR/03_candidate_catalogue/VC_catalogue.fna" ]]; then
  run_step build_candidate_catalogue python3 "$SCRIPT_DIR/helpers/build_virus_candidate_catalogue.py" --input-fasta "$OUTPUT_DIR/01_prepared_contigs/merged_assembled_contigs.fna" --provenance "$OUTPUT_DIR/01_prepared_contigs/contig_provenance.tsv" --genomad-fasta "$OUTPUT_DIR/02_genomad/viral_candidates_genomad.fna" --virsorter-fasta "$OUTPUT_DIR/02b_virsorter2/viral_candidates_virsorter2.fna" --diamond-virus-hits "$OUTPUT_DIR/02c_diamond_virus/nr_virus_hits.tsv" --output-dir "$OUTPUT_DIR/03_candidate_catalogue"
elif (( ! RESUME )); then die 4 'Candidate catalogue already exists; use --resume or a new output directory'; fi
launch_megan_background
step=(bash "$SCRIPT_DIR/10_run_diamond_taxonomy.sh" --input "$OUTPUT_DIR/03_candidate_catalogue/VC_catalogue.fna" --output-dir "$OUTPUT_DIR" --stage-dir 04_nr_annotation --threads "$DIAMOND_THREADS" --taxon-scope none --max-target-seqs "$DIAMOND_NR_MAX_TARGET_SEQS" --block-size "$DIAMOND_BLOCK_SIZE_RUN" --index-chunks "$DIAMOND_INDEX_CHUNKS_RUN" --tmpdir "$DIAMOND_TMPDIR_RUN"); (( RESUME )) && step+=(--resume); run_step diamond_nr_taxonomy "${step[@]}"
if [[ ! -s "$OUTPUT_DIR/04_nr_annotation/viral_decision.tsv" ]]; then
  run_step resolve_viral_evidence python3 "$SCRIPT_DIR/helpers/resolve_viral_decision.py" --catalogue-fasta "$OUTPUT_DIR/03_candidate_catalogue/VC_catalogue.fna" --discovery-evidence "$OUTPUT_DIR/03_candidate_catalogue/VC_discovery_evidence.tsv" --taxonomy "$OUTPUT_DIR/04_nr_annotation/contig_taxonomy_lca.tsv" --output-dir "$OUTPUT_DIR/04_nr_annotation"
elif (( ! RESUME )); then die 4 'Viral decision table already exists; use --resume or a new output directory'; fi
step=(bash "$SCRIPT_DIR/06_run_checkv.sh" --input "$OUTPUT_DIR/04_nr_annotation/checkv_input.fna" --output-dir "$OUTPUT_DIR" --stage-dir 05_checkv --database "$CHECKV_DB" --threads "$THREADS"); (( RESUME )) && step+=(--resume); run_step checkv "${step[@]}"
if [[ ! -f "$OUTPUT_DIR/05_checkv/ictv_selection.tsv" ]]; then
  run_step select_ictv_candidates python3 "$SCRIPT_DIR/helpers/select_ictv_candidates.py" --checkv-fasta "$OUTPUT_DIR/05_checkv/viral_candidates_checkv.fna" --decision "$OUTPUT_DIR/04_nr_annotation/viral_decision.tsv" --taxonomy "$OUTPUT_DIR/04_nr_annotation/contig_taxonomy_lca.tsv" --output-fasta "$OUTPUT_DIR/05_checkv/ictv_input.fna" --output-table "$OUTPUT_DIR/05_checkv/ictv_selection.tsv"
elif (( ! RESUME )); then die 4 'ICTV selection already exists; use --resume or a new output directory'; fi
if [[ -s "$OUTPUT_DIR/05_checkv/ictv_input.fna" ]]; then
  step=(bash "$SCRIPT_DIR/11_refine_ictv.sh" --input "$OUTPUT_DIR/05_checkv/ictv_input.fna" --output-dir "$OUTPUT_DIR" --threads "$DIAMOND_THREADS" --block-size "$DIAMOND_BLOCK_SIZE_RUN" --index-chunks "$DIAMOND_INDEX_CHUNKS_RUN" --tmpdir "$DIAMOND_TMPDIR_RUN"); (( RESUME )) && step+=(--resume); run_step ictv_refinement "${step[@]}"
else
  mkdir -p "$OUTPUT_DIR/06_ictv_refinement"
  printf 'qseqid\tqlen\tqstart\tqend\tpident\tlength\tevalue\tbitscore\tsseqid\n' > "$OUTPUT_DIR/06_ictv_refinement/ictv_hits.tsv"
  cp -- "$ICTV_REFERENCE_METADATA" "$OUTPUT_DIR/06_ictv_refinement/ictv_reference_metadata.tsv"
  echo '[INFO] No family-classified fragments; ICTV refinement is not applicable'
fi
if [[ ! -s "$OUTPUT_DIR/07_final_catalogue/VF_catalogue.fna" ]]; then
  run_step build_final_catalogue python3 "$SCRIPT_DIR/helpers/build_final_virome_catalogue.py" --checkv-fasta "$OUTPUT_DIR/05_checkv/viral_candidates_checkv.fna" --checkv-quality "$OUTPUT_DIR/05_checkv/quality_summary.tsv" --decision "$OUTPUT_DIR/04_nr_annotation/viral_decision.tsv" --nr-taxonomy "$OUTPUT_DIR/04_nr_annotation/contig_taxonomy_lca.tsv" --source-mapping "$OUTPUT_DIR/03_candidate_catalogue/VC_source_mapping.tsv" --ictv-hits "$OUTPUT_DIR/06_ictv_refinement/ictv_hits.tsv" --ictv-metadata "$OUTPUT_DIR/06_ictv_refinement/ictv_reference_metadata.tsv" --manifest "$MANIFEST" --catalogue-dir "$OUTPUT_DIR/07_final_catalogue" --sample-dir "$OUTPUT_DIR/08_sample_results"
elif (( ! RESUME )); then die 4 'Final catalogue already exists; use --resume or a new output directory'; fi
step=(bash "$SCRIPT_DIR/12_quantify_final_fragments.sh" --manifest "$MANIFEST" --sample-dir "$OUTPUT_DIR/08_sample_results" --output-dir "$OUTPUT_DIR/09_abundance" --threads "$THREADS"); (( RESUME )) && step+=(--resume); run_step quantify_fragments "${step[@]}"
run_step report python3 "$SCRIPT_DIR/helpers/build_virome_catalogue_report.py" --output-dir "$OUTPUT_DIR"
if (( MEGAN_DEFERRED )); then
  echo "[INFO] Foreground analysis completed; starting deferred MEGAN DAA/RMA6 auxiliary task"
  launch_megan_background 1
fi
printf 'SUCCESS\n' > "$RUN_DIR/status"; echo "[INFO] Virome catalogue completed: $OUTPUT_DIR/reports/virome_catalogue_dashboard.html"
