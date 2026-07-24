#!/usr/bin/env bash
# Generate NR DIAMOND DAA and, when configured, an MEGAN RMA6 file.
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"; load_pipeline_config

usage() { cat <<'EOF'
Usage: 09_run_diamond_megan.sh --input FASTA --output-dir PATH --threads N
       --taxon-scope virus|none|custom [--taxonlist IDs]
       [--max-target-seqs N] [--resume]
EOF
}

INPUT='' OUTPUT_DIR='' THREADS='' TAXON_SCOPE='virus' TAXONLIST='' MAX_TARGETS="$DIAMOND_NR_MAX_TARGET_SEQS" RESUME=0
while (($#)); do
  case "$1" in
    --input|--output-dir|--threads|--taxon-scope|--taxonlist|--max-target-seqs)
      (($# >= 2)) || die 2 "Missing value for $1"
      case "$1" in
        --input) INPUT=$2;; --output-dir) OUTPUT_DIR=$2;; --threads) THREADS=$2;; --taxon-scope) TAXON_SCOPE=$2;; --taxonlist) TAXONLIST=$2;; --max-target-seqs) MAX_TARGETS=$2;;
      esac
      shift 2;;
    --resume) RESUME=1; shift;;
    -h|--help) usage; exit 0;;
    *) die 2 "Unknown option: $1";;
  esac
done
[[ -s $INPUT && -n $OUTPUT_DIR && -n $THREADS ]] || die 2 'Input FASTA, output directory and threads are required'
[[ $TAXON_SCOPE == virus || $TAXON_SCOPE == none || $TAXON_SCOPE == custom ]] || die 2 'Invalid --taxon-scope'
positive_int "$THREADS" && positive_int "$MAX_TARGETS" || die 2 'Threads and max-target-seqs must be positive integers'
(( THREADS <= MAX_THREADS_PER_VIRAL_TOOL && THREADS <= MAX_TOTAL_THREADS )) || die 2 'DIAMOND thread request exceeds configured limits'
[[ -n $DIAMOND_NR_DB && -f $DIAMOND_NR_DB ]] || die 3 'DIAMOND_NR_DB is missing or not a readable .dmnd file'
require_command diamond

case "$TAXON_SCOPE" in
  virus) TAXONLIST="${TAXONLIST:-$DIAMOND_DEFAULT_TAXONLIST}";;
  none) TAXONLIST='';;
  custom) [[ -n $TAXONLIST ]] || die 2 'Custom taxon scope requires --taxonlist';;
esac
[[ -z $TAXONLIST ]] || valid_taxonlist "$TAXONLIST" || die 2 'Taxon list must contain positive NCBI TaxIDs separated by commas'

OUT="$OUTPUT_DIR/01_diamond_megan"; DAA="$OUT/viral_candidates.nr.daa"; RMA="$OUT/viral_candidates.nr.rma6"
if [[ -s $DAA && -s $RMA ]]; then
  (( RESUME )) && { echo "[INFO] DIAMOND DAA and MEGAN RMA6 already exist; skipped"; exit 0; }
  die 4 "DIAMOND/MEGAN output already exists; use --resume or choose another output directory"
fi
if [[ -e $OUT ]] && ! dir_is_empty_or_missing "$OUT"; then die 4 "DIAMOND/MEGAN output is incomplete or conflicting: $OUT"; fi
[[ -n $MEGAN_DAA2RMA ]] || die 3 'MEGAN_DAA2RMA is not configured; RMA6 generation requires daa2rma'
[[ -x $MEGAN_DAA2RMA || $(command -v "$MEGAN_DAA2RMA" 2>/dev/null) ]] || die 127 "MEGAN daa2rma is not executable: $MEGAN_DAA2RMA"
[[ -n $MEGAN_MAP_DB && -f $MEGAN_MAP_DB ]] || die 3 'MEGAN_MAP_DB is missing or not a readable MEGAN mapping database'

mkdir -p "$OUT"
diamond_args=(blastx --db "$DIAMOND_NR_DB" --query "$INPUT" --out "$DAA" --outfmt 100 --threads "$THREADS" --evalue "$DIAMOND_EVALUE" --max-target-seqs "$MAX_TARGETS")
[[ $DIAMOND_SENSITIVITY == more-sensitive ]] && diamond_args+=(--more-sensitive)
[[ -z $TAXONLIST ]] || diamond_args+=(--taxonlist "$TAXONLIST")
printf '%q ' diamond "${diamond_args[@]}" > "$OUT/diamond_command.sh"; printf '\n' >> "$OUT/diamond_command.sh"
diamond "${diamond_args[@]}" >"$OUT/diamond.stdout.log" 2>"$OUT/diamond.stderr.log"
[[ -s $DAA ]] || die 1 "DIAMOND did not produce a DAA file; inspect: $OUT/diamond.stderr.log"
"$MEGAN_DAA2RMA" -i "$DAA" -o "$RMA" -mdb "$MEGAN_MAP_DB" -sup 1 -t "$THREADS" -v >"$OUT/daa2rma.stdout.log" 2>"$OUT/daa2rma.stderr.log"
[[ -s $RMA ]] || die 1 "daa2rma did not produce RMA6; inspect: $OUT/daa2rma.stderr.log"
printf 'input=%s\ntaxon_scope=%s\ntaxonlist=%s\nmax_target_seqs=%s\n' "$INPUT" "$TAXON_SCOPE" "${TAXONLIST:-all_nr}" "$MAX_TARGETS" > "$OUT/parameters.env"
echo "[INFO] MEGAN RMA6: $RMA"
