#!/usr/bin/env bash
# Shared helpers. This file is sourced by the executable scripts.
set -o pipefail

PIPELINE_EXIT_GENERAL=1; PIPELINE_EXIT_USAGE=2; PIPELINE_EXIT_INPUT=3
PIPELINE_EXIT_CONFLICT=4; PIPELINE_EXIT_PARTIAL=5; PIPELINE_EXIT_LOCKED=75; PIPELINE_EXIT_MISSING_TOOL=127
SCRIPT_LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_HOME="$(cd "$SCRIPT_LIB_DIR/../.." && pwd)"
CONFIG_FILE="${CONTIG_PIPELINE_CONFIG:-$PIPELINE_HOME/config/pipeline.env}"

load_pipeline_config() {
  [[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE"
  : "${ALLOWED_DATA_ROOTS:=}"; : "${MAX_TOTAL_THREADS:=96}"; : "${MAX_QC_PARALLEL:=8}"; : "${MAX_ASSEMBLY_PARALLEL:=4}"
  : "${MAX_THREADS_PER_MEGHIT:=48}"; : "${MAX_THREADS_PER_FASTP:=32}"
  : "${ADAPTER_CATALOG:=$PIPELINE_HOME/config/adapter_catalog.tsv}"; : "${FASTP_DEFAULT_ADAPTER_PROFILE:=auto}"
  : "${ADAPTER_SEQUENCE_REFERENCE:=$PIPELINE_HOME/config/adapter_sequence_reference.tsv}"
  : "${ADAPTER_REFERENCE_SCAN_READS:=100000}"
  : "${FASTP_QUALIFIED_QUALITY_PHRED:=6}"; : "${FASTP_UNQUALIFIED_PERCENT_LIMIT:=50}"; : "${FASTP_N_BASE_LIMIT:=15}"; : "${FASTP_LENGTH_REQUIRED:=50}"
  : "${GENOMAD_DB:=}"; : "${CHECKV_DB:=}"; : "${MAX_THREADS_PER_VIRAL_TOOL:=32}"
  : "${VIRAL_MIN_CONTIG_LEN:=1000}"; : "${VOTU_POST_CHECKV_MIN_LEN:=1000}"; : "${VOTU_ANI:=95}"; : "${VOTU_ALIGNED_FRACTION:=85}"
  : "${COVERM_MIN_READ_PERCENT_IDENTITY:=95}"; : "${COVERM_MIN_READ_ALIGNED_PERCENT:=75}"; : "${COVERM_MIN_COVERED_FRACTION:=10}"
  : "${VOTU_IMPORTANCE_RELATIVE_ABUNDANCE:=5}"; : "${BATCH_OVERVIEW_RANK:=family}"
  : "${REPORT_THEME:=quarto-scientific}"; : "${REPORT_TOP_TAXA:=0}"
  : "${REPORT_PRIORITY_FAMILIES:=Coronaviridae,Paramyxoviridae,Orthomyxoviridae,Flaviviridae}"
  : "${ICTV_FAMILY_GENOME_REFERENCE:=$PIPELINE_HOME/config/ictv_family_genome_reference.tsv}"
  : "${PRIORITY_REVIEW_TAXA_REFERENCE:=$PIPELINE_HOME/config/priority_review_taxa.tsv}"
  : "${DIAMOND_NR_DB:=}"; : "${DIAMOND_DEFAULT_TAXONLIST:=10239}"; : "${DIAMOND_EVALUE:=1e-5}"
  : "${DIAMOND_NR_MAX_TARGET_SEQS:=25}"; : "${DIAMOND_SENSITIVITY:=more-sensitive}"
  : "${MEGAN_DAA2RMA:=}"; : "${MEGAN_MAP_DB:=}"; : "${TAXONKIT_DB:=}"
  : "${VIRSORTER_COMMAND:=virsorter}"; : "${VIRSORTER_USE_CONDA_OFF:=0}"
}
die() { local code=$1; shift; printf '[ERROR] %s\n' "$*" >&2; exit "$code"; }
require_command() { command -v "$1" >/dev/null 2>&1 || die "$PIPELINE_EXIT_MISSING_TOOL" "Required command is not in PATH: $1"; }
require_executable() {
  local executable=$1
  if [[ $executable == */* ]]; then
    [[ -x $executable ]] || die "$PIPELINE_EXIT_MISSING_TOOL" "Required executable is unavailable: $executable"
  else
    require_command "$executable"
  fi
}
positive_int() { [[ ${1:-} =~ ^[1-9][0-9]*$ ]]; }
valid_taxonlist() { [[ ${1:-} =~ ^[1-9][0-9]*(,[1-9][0-9]*)*$ ]]; }
path_is_allowed() {
  local candidate=$1 root real_root
  [[ -n "$ALLOWED_DATA_ROOTS" ]] || return 0
  candidate=$(realpath -m -- "$candidate") || return 1
  IFS=':' read -r -a roots <<< "$ALLOWED_DATA_ROOTS"
  for root in "${roots[@]}"; do
    [[ -n "$root" ]] || continue; real_root=$(realpath -m -- "$root") || continue
    [[ "$candidate" == "$real_root" || "$candidate" == "$real_root"/* ]] && return 0
  done
  return 1
}
require_allowed_path() { path_is_allowed "$2" || die "$PIPELINE_EXIT_INPUT" "$1 is outside ALLOWED_DATA_ROOTS: $2"; }
assert_existing_dir() { [[ -d "$2" ]] || die "$PIPELINE_EXIT_INPUT" "$1 does not exist or is not a directory: $2"; }
assert_output_parent() { local parent; parent=$(dirname "$2"); [[ -d "$parent" && -w "$parent" ]] || die "$PIPELINE_EXIT_INPUT" "Parent directory for $1 does not exist or is not writable: $parent"; }
dir_is_empty_or_missing() { [[ ! -e "$1" ]] && return 0; [[ -d "$1" ]] || return 1; [[ -z $(find "$1" -mindepth 1 -print -quit 2>/dev/null) ]]; }
valid_sample_id() { [[ $1 =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; }

R1_SUFFIXES=('_R1.clean.fq.gz' '_R1.clean.fastq.gz' '_1.clean.fq.gz' '_1.clean.fastq.gz' '_R1_001.fq.gz' '_R1_001.fastq.gz' '_R1.fq.gz' '_R1.fastq.gz' '_1.fq.gz' '_1.fastq.gz')
R2_SUFFIXES=('_R2.clean.fq.gz' '_R2.clean.fastq.gz' '_2.clean.fq.gz' '_2.clean.fastq.gz' '_R2_001.fq.gz' '_R2_001.fastq.gz' '_R2.fq.gz' '_R2.fastq.gz' '_2.fq.gz' '_2.fastq.gz')
# Result variables: RESOLVED_R1, RESOLVED_R2, RESOLVED_STEM.
resolve_unique_pe_pair_in_dir() {
  local source_dir=$1 i file name stem expected; local -a r1_matches=() r2_matches=() stems=()
  while IFS= read -r -d '' file; do
    name=${file##*/}
    for ((i=0; i<${#R1_SUFFIXES[@]}; i++)); do
      [[ "$name" == *"${R1_SUFFIXES[$i]}" ]] || continue; stem=${name%"${R1_SUFFIXES[$i]}"}; expected="$source_dir/${stem}${R2_SUFFIXES[$i]}"
      [[ -f "$expected" ]] || continue; r1_matches+=("$file"); r2_matches+=("$expected"); stems+=("$stem")
    done
  done < <(find "$source_dir" -maxdepth 1 -type f \( -iname '*.fq.gz' -o -iname '*.fastq.gz' \) -print0)
  (( ${#r1_matches[@]} == 1 )) || return 1
  RESOLVED_R1=${r1_matches[0]}; RESOLVED_R2=${r2_matches[0]}; RESOLVED_STEM=${stems[0]}
}
# Result variables: RESOLVED_SE, RESOLVED_STEM.
resolve_unique_se_file_in_dir() {
  local source_dir=$1 file name stem; local -a files=() stems=()
  while IFS= read -r -d '' file; do
    name=${file##*/}; case "$name" in *.clean.fq.gz) stem=${name%.clean.fq.gz};; *.clean.fastq.gz) stem=${name%.clean.fastq.gz};; *.fq.gz) stem=${name%.fq.gz};; *.fastq.gz) stem=${name%.fastq.gz};; *) continue;; esac
    files+=("$file"); stems+=("$stem")
  done < <(find "$source_dir" -maxdepth 1 -type f \( -iname '*.fq.gz' -o -iname '*.fastq.gz' \) -print0)
  (( ${#files[@]} == 1 )) || return 1
  RESOLVED_SE=${files[0]}; RESOLVED_STEM=${stems[0]}
}
manifest_header() { printf 'sample_id\tread_type\traw_r1\traw_r2\tclean_r1\tclean_r2\tclean_single\tassembly_dir\n'; }
manifest_add_pe() { printf '%s\tpe\t%s\t%s\t%s\t%s\t\t%s\n' "$2" "$3" "$4" "$5" "$6" "$7" >> "$1"; }
manifest_add_se() { printf '%s\tse\t\t\t\t\t%s\t%s\n' "$2" "$3" "$4" >> "$1"; }
manifest_rows() { tail -n +2 "$1"; }
