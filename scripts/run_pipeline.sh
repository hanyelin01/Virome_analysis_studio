#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"; load_pipeline_config
usage() { cat <<'EOF'
Usage: run_pipeline.sh --task qc_only|assembly_only|full --cleandata-dir PATH
  [--rawdata-dir PATH --raw-layout sample_subdirs|flat]
  [--clean-layout sample_subdirs|flat --read-type pe|se --assembly-dir PATH]
  [--qc-parallel N --qc-threads N --assembly-parallel N --assembly-threads N]
  [--min-contig-len N] [--resume]
EOF
}
TASK='' RAW_DIR='' RAW_LAYOUT='' CLEAN_DIR='' CLEAN_LAYOUT='' READ_TYPE='' ASSEMBLY_DIR='' QC_PARALLEL=4 QC_THREADS=8 ASM_PARALLEL=2 ASM_THREADS=30 MINLEN=40 RESUME=0
while (($#)); do
  case "$1" in
    --task|--rawdata-dir|--raw-layout|--cleandata-dir|--clean-layout|--read-type|--assembly-dir|--qc-parallel|--qc-threads|--assembly-parallel|--assembly-threads|--min-contig-len)
      (($#>=2)) || die 2 "Missing value for $1"
      case "$1" in --task) TASK=$2;; --rawdata-dir) RAW_DIR=$2;; --raw-layout) RAW_LAYOUT=$2;; --cleandata-dir) CLEAN_DIR=$2;; --clean-layout) CLEAN_LAYOUT=$2;; --read-type) READ_TYPE=$2;; --assembly-dir) ASSEMBLY_DIR=$2;; --qc-parallel) QC_PARALLEL=$2;; --qc-threads) QC_THREADS=$2;; --assembly-parallel) ASM_PARALLEL=$2;; --assembly-threads) ASM_THREADS=$2;; --min-contig-len) MINLEN=$2;; esac
      shift 2;;
    --resume) RESUME=1; shift;; -h|--help) usage; exit 0;; *) die 2 "Unknown option: $1";;
  esac
done
[[ $TASK == qc_only || $TASK == assembly_only || $TASK == full ]] || die 2 'Invalid --task'
positive_int "$QC_PARALLEL" && positive_int "$QC_THREADS" && positive_int "$ASM_PARALLEL" && positive_int "$ASM_THREADS" && positive_int "$MINLEN" || die 2 'Numeric options must be positive integers'
[[ $TASK != full ]] || { CLEAN_LAYOUT=sample_subdirs; READ_TYPE=pe; }
STATE_BASE=${ASSEMBLY_DIR:-$CLEAN_DIR}; [[ -n $STATE_BASE ]] || die 2 '--cleandata-dir is required'
require_allowed_path 'pipeline state directory' "$STATE_BASE"
mkdir -p "$STATE_BASE/.contig_pipeline/runs" "$STATE_BASE/.contig_pipeline/reports"
RUN_ID="$(date '+%Y%m%d_%H%M%S')_$$"; RUN_DIR="$STATE_BASE/.contig_pipeline/runs/$RUN_ID"; REPORT_DIR="$STATE_BASE/.contig_pipeline/reports"; MANIFEST="$RUN_DIR/sample_manifest.tsv"
mkdir -p "$RUN_DIR"; exec 9>"$STATE_BASE/.contig_pipeline/.pipeline.lock"; flock -n 9 || die 75 "Another pipeline is already running for this output location: $STATE_BASE"
LOG="$RUN_DIR/pipeline.log"; exec > >(tee -a "$LOG") 2>&1; printf 'RUNNING\n' > "$RUN_DIR/status"
cat > "$RUN_DIR/parameters.env" <<EOF
TASK=$TASK
RAW_DIR=$RAW_DIR
RAW_LAYOUT=$RAW_LAYOUT
CLEAN_DIR=$CLEAN_DIR
CLEAN_LAYOUT=$CLEAN_LAYOUT
READ_TYPE=$READ_TYPE
ASSEMBLY_DIR=$ASSEMBLY_DIR
QC_PARALLEL=$QC_PARALLEL
QC_THREADS=$QC_THREADS
ASSEMBLY_PARALLEL=$ASM_PARALLEL
ASSEMBLY_THREADS=$ASM_THREADS
MIN_CONTIG_LEN=$MINLEN
RESUME=$RESUME
EOF
run_step() { local label=$1; shift; echo "[STEP] $label"; "$@" 2>&1 | tee -a "$RUN_DIR/${label}.log"; }
on_error() { local rc=$?; printf 'FAILED\n' > "$RUN_DIR/status"; echo "[ERROR] Pipeline stopped; log: $LOG"; exit "$rc"; }
trap on_error ERR
preflight=(bash "$SCRIPT_DIR/00_preflight.sh" --task "$TASK" --cleandata-dir "$CLEAN_DIR" --manifest "$MANIFEST")
[[ $TASK == qc_only || $TASK == full ]] && preflight+=(--rawdata-dir "$RAW_DIR" --raw-layout "$RAW_LAYOUT")
[[ $TASK == assembly_only || $TASK == full ]] && preflight+=(--clean-layout "$CLEAN_LAYOUT" --read-type "$READ_TYPE" --assembly-dir "$ASSEMBLY_DIR")
run_step preflight "${preflight[@]}"
if [[ $TASK == qc_only || $TASK == full ]]; then
  cmd=(bash "$SCRIPT_DIR/01_fastp_pe.sh" --manifest "$MANIFEST" --parallel-samples "$QC_PARALLEL" --threads-per-sample "$QC_THREADS"); (( RESUME )) && cmd+=(--resume); run_step fastp "${cmd[@]}"
fi
if [[ $TASK == assembly_only || $TASK == full ]]; then
  cmd=(bash "$SCRIPT_DIR/02_megahit.sh" --manifest "$MANIFEST" --parallel-samples "$ASM_PARALLEL" --threads-per-sample "$ASM_THREADS" --min-contig-len "$MINLEN"); (( RESUME )) && cmd+=(--resume); run_step megahit "${cmd[@]}"
  run_step check_contigs bash "$SCRIPT_DIR/03_check_contigs.sh" --manifest "$MANIFEST" --output "$REPORT_DIR/assembly_summary_${RUN_ID}.tsv"
fi
printf 'SUCCESS\n' > "$RUN_DIR/status"; echo "[INFO] Pipeline completed. Run directory: $RUN_DIR"
