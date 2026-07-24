#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"; load_pipeline_config
MANIFEST='' PARALLEL='' THREADS='' ADAPTER_PROFILE=$FASTP_DEFAULT_ADAPTER_PROFILE ADAPTER_CATALOG=$ADAPTER_CATALOG RESUME=0
while (($#)); do
  case "$1" in
    --manifest|--parallel-samples|--threads-per-sample|--adapter-profile|--adapter-catalog)
      (($#>=2)) || die 2 "Missing value for $1"; case "$1" in --manifest) MANIFEST=$2;; --parallel-samples) PARALLEL=$2;; --threads-per-sample) THREADS=$2;; --adapter-profile) ADAPTER_PROFILE=$2;; --adapter-catalog) ADAPTER_CATALOG=$2;; esac; shift 2;;
    --resume) RESUME=1; shift;; *) die 2 "Unknown option: $1";;
  esac
done
[[ -f $MANIFEST ]] || die 2 '--manifest is required'
positive_int "$PARALLEL" && positive_int "$THREADS" || die 2 'Parallel and threads must be positive integers'
(( PARALLEL <= MAX_QC_PARALLEL && THREADS <= MAX_THREADS_PER_FASTP && PARALLEL * THREADS <= MAX_TOTAL_THREADS )) || die 2 'fastp resource request exceeds configured limits'
require_command fastp
require_command python3
python3 "$SCRIPT_DIR/helpers/adapter_evidence.py" validate --catalog "$ADAPTER_CATALOG" >/dev/null
R1_ADAPTER=$(python3 "$SCRIPT_DIR/helpers/adapter_evidence.py" get --catalog "$ADAPTER_CATALOG" --profile "$ADAPTER_PROFILE" --field r1_sequence)
R2_ADAPTER=$(python3 "$SCRIPT_DIR/helpers/adapter_evidence.py" get --catalog "$ADAPTER_CATALOG" --profile "$ADAPTER_PROFILE" --field r2_sequence)
run_sample() {
  local sample=$1 r1=$2 r2=$3 out1=$4 out2=$5 outdir report evidence log; local -a fastp_cmd
  outdir=$(dirname "$out1"); report="$outdir/fastp_report/${sample}.fastp.json"; log="$outdir/log"
  evidence="$outdir/fastp_report/${sample}.adapter_evidence.tsv"
  if [[ -s $out1 || -s $out2 || -s $report ]]; then
    if (( RESUME )) && [[ -s $out1 && -s $out2 && -s $report ]]; then
      [[ -s $evidence ]] || python3 "$SCRIPT_DIR/helpers/adapter_evidence.py" summarize --catalog "$ADAPTER_CATALOG" --profile "$ADAPTER_PROFILE" --fastp-json "$report" --sample "$sample" --output "$evidence"
      echo "[INFO] $sample: fastp already complete; skipped"; return 0
    fi
    echo "[ERROR] $sample: incomplete or conflicting fastp output: $outdir" >&2; return 4
  fi
  if [[ -e $outdir ]] && ! dir_is_empty_or_missing "$outdir"; then echo "[ERROR] $sample: output directory is not empty: $outdir" >&2; return 4; fi
  mkdir -p "$outdir/fastp_report" "$log"; echo "[INFO] $sample: fastp started"
  fastp_cmd=(fastp --in1 "$r1" --in2 "$r2" --out1 "$out1" --out2 "$out2" --html "$outdir/fastp_report/${sample}.fastp.html" --json "$report" --report_title "${sample} fastp QC" --thread "$THREADS" --detect_adapter_for_pe --qualified_quality_phred "$FASTP_QUALIFIED_QUALITY_PHRED" --unqualified_percent_limit "$FASTP_UNQUALIFIED_PERCENT_LIMIT" --n_base_limit "$FASTP_N_BASE_LIMIT" --length_required "$FASTP_LENGTH_REQUIRED" --disable_trim_poly_g --dont_overwrite)
  if [[ -n $R1_ADAPTER ]]; then fastp_cmd+=(--adapter_sequence "$R1_ADAPTER" --adapter_sequence_r2 "$R2_ADAPTER"); fi
  "${fastp_cmd[@]}" >"$log/${sample}.fastp.stdout.log" 2>"$log/${sample}.fastp.stderr.log"
  python3 "$SCRIPT_DIR/helpers/adapter_evidence.py" summarize --catalog "$ADAPTER_CATALOG" --profile "$ADAPTER_PROFILE" --fastp-json "$report" --sample "$sample" --output "$evidence"
  echo "[INFO] $sample: fastp completed"
}
pids=(); samples=(); status=0
wait_one() { local rc=0; wait "${pids[0]}" || rc=$?; (( rc == 0 )) || { echo "[ERROR] ${samples[0]}: fastp failed (exit $rc)" >&2; status=5; }; pids=("${pids[@]:1}"); samples=("${samples[@]:1}"); }
while IFS= read -r manifest_row; do
  IFS=$'\x1f' read -r sample type raw1 raw2 clean1 clean2 clean_single assembly <<< "${manifest_row//$'\t'/$'\x1f'}"
  [[ $type == pe ]] || continue
  run_sample "$sample" "$raw1" "$raw2" "$clean1" "$clean2" & pids+=("$!"); samples+=("$sample")
  (( ${#pids[@]} < PARALLEL )) || wait_one
done < <(manifest_rows "$MANIFEST")
while (( ${#pids[@]} )); do wait_one; done
exit "$status"
