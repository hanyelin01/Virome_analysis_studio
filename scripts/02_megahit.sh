#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"; load_pipeline_config
MANIFEST='' PARALLEL='' THREADS='' MINLEN='' RESUME=0
while (($#)); do
  case "$1" in
    --manifest|--parallel-samples|--threads-per-sample|--min-contig-len)
      (($#>=2)) || die 2 "Missing value for $1"; case "$1" in --manifest) MANIFEST=$2;; --parallel-samples) PARALLEL=$2;; --threads-per-sample) THREADS=$2;; --min-contig-len) MINLEN=$2;; esac; shift 2;;
    --resume) RESUME=1; shift;; *) die 2 "Unknown option: $1";;
  esac
done
[[ -f $MANIFEST ]] || die 2 '--manifest is required'
positive_int "$PARALLEL" && positive_int "$THREADS" && positive_int "$MINLEN" || die 2 'Numeric arguments must be positive integers'
if ! (( PARALLEL <= MAX_ASSEMBLY_PARALLEL && THREADS <= MAX_THREADS_PER_MEGHIT && PARALLEL * THREADS <= MAX_TOTAL_THREADS )); then
  die 2 "MEGAHIT resource request exceeds configured limits: requested parallel=${PARALLEL}, threads_per_sample=${THREADS}, total=$(( PARALLEL * THREADS )); limits parallel<=${MAX_ASSEMBLY_PARALLEL}, threads_per_sample<=${MAX_THREADS_PER_MEGHIT}, total<=${MAX_TOTAL_THREADS}"
fi
require_command megahit
run_sample() {
  local sample=$1 type=$2 r1=$3 r2=$4 se=$5 out=$6 log final
  # MEGAHIT must create its own output directory. Logs therefore live beside
  # the sample directories rather than inside an output directory created early.
  log="$(dirname "$out")/log/$sample"; final="$out/final.contigs.fa"
  if [[ -s $final ]]; then
    if (( RESUME )); then echo "[INFO] $sample: assembly already complete; skipped"; return 0; fi
    echo "[ERROR] $sample: output exists; use --resume or move it aside" >&2; return 4
  fi
  if [[ -e $out ]] && ! dir_is_empty_or_missing "$out"; then echo "[ERROR] $sample: incomplete/conflicting assembly output: $out" >&2; return 4; fi
  mkdir -p "$log"; echo "[INFO] $sample: MEGAHIT started ($type)"
  if [[ $type == pe ]]; then
    megahit -1 "$r1" -2 "$r2" -o "$out" --min-contig-len "$MINLEN" -t "$THREADS" >"$log/${sample}.megahit.stdout.log" 2>"$log/${sample}.megahit.stderr.log"
  else
    megahit -r "$se" -o "$out" --min-contig-len "$MINLEN" -t "$THREADS" >"$log/${sample}.megahit.stdout.log" 2>"$log/${sample}.megahit.stderr.log"
  fi
  [[ -s $final ]] || { echo "[ERROR] $sample: final.contigs.fa is missing or empty" >&2; return 1; }
  echo "[INFO] $sample: MEGAHIT completed"
}
pids=(); samples=(); status=0
wait_one() { local rc=0; wait "${pids[0]}" || rc=$?; (( rc == 0 )) || { echo "[ERROR] ${samples[0]}: assembly failed (exit $rc)" >&2; status=5; }; pids=("${pids[@]:1}"); samples=("${samples[@]:1}"); }
while IFS= read -r manifest_row; do
  # Tab is whitespace in Bash IFS and would collapse empty TSV cells. Convert
  # it to a non-whitespace delimiter so PE and SE columns retain their position.
  IFS=$'\x1f' read -r sample type raw1 raw2 clean1 clean2 clean_single assembly <<< "${manifest_row//$'\t'/$'\x1f'}"
  run_sample "$sample" "$type" "$clean1" "$clean2" "$clean_single" "$assembly" & pids+=("$!"); samples+=("$sample")
  (( ${#pids[@]} < PARALLEL )) || wait_one
done < <(manifest_rows "$MANIFEST")
while (( ${#pids[@]} )); do wait_one; done
exit "$status"
