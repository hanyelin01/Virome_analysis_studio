#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"; load_pipeline_config

OUTPUT_DIR='' MANIFEST='' GROUPS_FILE='' OVERVIEW_RANK="$BATCH_OVERVIEW_RANK" THEME="$REPORT_THEME" TOP_TAXA="$REPORT_TOP_TAXA" PRIORITY_FAMILIES="$REPORT_PRIORITY_FAMILIES" ICTV_REFERENCE="$ICTV_FAMILY_GENOME_REFERENCE" PRIORITY_REFERENCE="$PRIORITY_REVIEW_TAXA_REFERENCE" RESUME=0 REFRESH=0
while (($#)); do
  case "$1" in
    --output-dir|--manifest|--groups-file|--overview-rank|--theme|--top-taxa|--priority-families|--ictv-reference|--priority-reference)
      (($# >= 2)) || die 2 "Missing value for $1"
      case "$1" in --output-dir) OUTPUT_DIR=$2;; --manifest) MANIFEST=$2;; --groups-file) GROUPS_FILE=$2;; --overview-rank) OVERVIEW_RANK=$2;; --theme) THEME=$2;; --top-taxa) TOP_TAXA=$2;; --priority-families) PRIORITY_FAMILIES=$2;; --ictv-reference) ICTV_REFERENCE=$2;; --priority-reference) PRIORITY_REFERENCE=$2;; esac
      shift 2;;
    --resume) RESUME=1; shift;;
    --refresh) REFRESH=1; shift;;
    -h|--help) echo "Usage: $0 --output-dir PATH --manifest PATH [--overview-rank family|genus|species] [--theme quarto-scientific] [--top-taxa 0|N] [--ictv-reference TSV] [--priority-reference TSV] [--priority-families CSV] [--resume] [--refresh]"; exit 0;;
    *) die 2 "Unknown option: $1";;
  esac
done
[[ -n $OUTPUT_DIR && -f $MANIFEST ]] || die 2 '--output-dir and --manifest are required'
[[ $OVERVIEW_RANK == family || $OVERVIEW_RANK == genus || $OVERVIEW_RANK == species ]] || die 2 'overview rank must be family, genus or species'
[[ $THEME == quarto-scientific ]] || die 2 'report theme must be quarto-scientific'
[[ $TOP_TAXA =~ ^[0-9]+$ ]] || die 2 'top taxa must be zero or a positive integer'
REPORT="$OUTPUT_DIR/reports/batch_overview.html"
if [[ -s $REPORT && $REFRESH -eq 0 ]]; then
  (( RESUME )) && { echo "[INFO] HTML report already exists; skipped"; exit 0; }
  die 4 "Report already exists; use --resume or choose a new report output directory"
fi
PYTHON_EXEC="${CONTIG_PIPELINE_PYTHON:-$PIPELINE_HOME/.venv/bin/python}"
[[ -x $PYTHON_EXEC ]] || die 127 "Pipeline Python interpreter is unavailable: $PYTHON_EXEC"
args=(--output-dir "$OUTPUT_DIR" --manifest "$MANIFEST" --overview-rank "$OVERVIEW_RANK" --theme "$THEME" --top-taxa "$TOP_TAXA" --priority-families "$PRIORITY_FAMILIES")
[[ -z $GROUPS_FILE ]] || args+=(--groups-file "$GROUPS_FILE")
[[ -f $ICTV_REFERENCE ]] && args+=(--ictv-reference "$ICTV_REFERENCE") || echo "[WARN] ICTV reference is unavailable; report will mark genome groups as unclassified: $ICTV_REFERENCE" >&2
[[ -f $PRIORITY_REFERENCE ]] && args+=(--priority-reference "$PRIORITY_REFERENCE") || echo "[WARN] Attention-taxon reference is unavailable; falling back to REPORT_PRIORITY_FAMILIES" >&2
"$PYTHON_EXEC" "$SCRIPT_DIR/helpers/build_viral_report.py" "${args[@]}"
