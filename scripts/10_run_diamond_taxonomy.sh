#!/usr/bin/env bash
# Generate NR DIAMOND outfmt-6 and a TaxonKit LCA annotation table.
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"; load_pipeline_config

usage() { cat <<'EOF'
Usage: 10_run_diamond_taxonomy.sh --input FASTA --output-dir PATH --threads N
       --taxon-scope virus|none|custom [--taxonlist IDs]
       [--max-target-seqs N] [--resume]
EOF
}

INPUT='' OUTPUT_DIR='' THREADS='' TAXON_SCOPE='virus' TAXONLIST='' MAX_TARGETS="$DIAMOND_NR_MAX_TARGET_SEQS" STAGE_DIR='02_diamond_nr_taxonomy' RESUME=0
while (($#)); do
  case "$1" in
    --input|--output-dir|--threads|--taxon-scope|--taxonlist|--max-target-seqs|--stage-dir)
      (($# >= 2)) || die 2 "Missing value for $1"
      case "$1" in
        --input) INPUT=$2;; --output-dir) OUTPUT_DIR=$2;; --threads) THREADS=$2;; --taxon-scope) TAXON_SCOPE=$2;; --taxonlist) TAXONLIST=$2;; --max-target-seqs) MAX_TARGETS=$2;; --stage-dir) STAGE_DIR=$2;;
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
[[ $STAGE_DIR =~ ^[A-Za-z0-9._-]+$ ]] || die 2 '--stage-dir must be a simple directory name'
(( THREADS <= MAX_THREADS_PER_VIRAL_TOOL && THREADS <= MAX_TOTAL_THREADS )) || die 2 'DIAMOND thread request exceeds configured limits'
[[ -n $DIAMOND_NR_DB && -f $DIAMOND_NR_DB ]] || die 3 'DIAMOND_NR_DB is missing or not a readable .dmnd file'
[[ -n $TAXONKIT_DB && -d $TAXONKIT_DB ]] || die 3 'TAXONKIT_DB is missing or not a directory'
require_command diamond; require_command taxonkit
case "$TAXON_SCOPE" in
  virus) TAXONLIST="${TAXONLIST:-$DIAMOND_DEFAULT_TAXONLIST}";;
  none) TAXONLIST='';;
  custom) [[ -n $TAXONLIST ]] || die 2 'Custom taxon scope requires --taxonlist';;
esac
[[ -z $TAXONLIST ]] || valid_taxonlist "$TAXONLIST" || die 2 'Taxon list must contain positive NCBI TaxIDs separated by commas'

OUT="$OUTPUT_DIR/$STAGE_DIR"; HITS="$OUT/nr_virus_hits.outfmt6.tsv"; TAXIDS="$OUT/query_taxids.tsv"; LCA="$OUT/query_lca.tsv"; LINEAGE="$OUT/query_lca_lineage.tsv"; SUMMARY="$OUT/contig_taxonomy_lca.tsv"
if [[ -s $SUMMARY ]]; then
  (( RESUME )) && { echo "[INFO] DIAMOND/TaxonKit annotation already exists; skipped"; exit 0; }
  die 4 "DIAMOND/TaxonKit output already exists; use --resume or choose another output directory"
fi
if [[ -e $OUT ]] && ! dir_is_empty_or_missing "$OUT"; then
  (( RESUME )) || die 4 "DIAMOND/TaxonKit output is incomplete or conflicting: $OUT"
  echo "[INFO] Restarting incomplete DIAMOND/TaxonKit annotation step"
fi
mkdir -p "$OUT"

# All taxonomic interpretation is based on staxids and TaxonKit below.  Use
# fields supported by the deployed DIAMOND 2.0.x release; rank/lineage output
# fields such as slineages were introduced only in later DIAMOND versions.
fields=(qseqid qlen qstart qend pident length evalue bitscore sstart send sseqid staxids sscinames sskingdoms)
diamond_args=(blastx --db "$DIAMOND_NR_DB" --query "$INPUT" --out "$HITS" --outfmt 6 "${fields[@]}" --threads "$THREADS" --evalue "$DIAMOND_EVALUE" --max-target-seqs "$MAX_TARGETS")
[[ $DIAMOND_SENSITIVITY == more-sensitive ]] && diamond_args+=(--more-sensitive)
[[ -z $TAXONLIST ]] || diamond_args+=(--taxonlist "$TAXONLIST")
printf '%q ' diamond "${diamond_args[@]}" > "$OUT/diamond_command.sh"; printf '\n' >> "$OUT/diamond_command.sh"
diamond "${diamond_args[@]}" >"$OUT/diamond.stdout.log" 2>"$OUT/diamond.stderr.log"
: > "$TAXIDS"
python3 "$SCRIPT_DIR/helpers/diamond_hits_to_taxids.py" --hits "$HITS" --output "$TAXIDS"
if [[ $(wc -l < "$TAXIDS") -gt 1 ]]; then
  TAXONKIT_DB="$TAXONKIT_DB" taxonkit lca --taxids-field 2 --separator ',' --skip-deleted --skip-unfound "$TAXIDS" > "$LCA" 2>"$OUT/taxonkit_lca.stderr.log"
  TAXONKIT_DB="$TAXONKIT_DB" taxonkit lineage --taxid-field 3 --show-name --show-rank "$LCA" > "$LINEAGE" 2>"$OUT/taxonkit_lineage.stderr.log"
else
  : > "$LCA"; : > "$LINEAGE"
fi
python3 "$SCRIPT_DIR/helpers/summarize_diamond_taxonomy.py" --hits "$HITS" --lca "$LCA" --lineage "$LINEAGE" --output "$SUMMARY"
printf 'input=%s\ntaxon_scope=%s\ntaxonlist=%s\nmax_target_seqs=%s\n' "$INPUT" "$TAXON_SCOPE" "${TAXONLIST:-all_nr}" "$MAX_TARGETS" > "$OUT/parameters.env"
echo "[INFO] DIAMOND/TaxonKit annotation: $SUMMARY"
