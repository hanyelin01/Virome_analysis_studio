#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"
load_pipeline_config

usage() { cat <<'EOF'
Usage: 00_preflight.sh --task qc_only|assembly_only|full --cleandata-dir PATH --manifest PATH
       [--rawdata-dir PATH --raw-layout sample_subdirs|flat]
       [--clean-layout sample_subdirs|flat --read-type pe|se --assembly-dir PATH]
EOF
}
TASK='' RAW_DIR='' RAW_LAYOUT='' CLEAN_DIR='' CLEAN_LAYOUT='' READ_TYPE='' ASSEMBLY_DIR='' MANIFEST=''
while (($#)); do
  case "$1" in
    --task|--rawdata-dir|--raw-layout|--cleandata-dir|--clean-layout|--read-type|--assembly-dir|--manifest)
      (($# >= 2)) || die 2 "Missing value for $1"
      case "$1" in
        --task) TASK=$2;; --rawdata-dir) RAW_DIR=$2;; --raw-layout) RAW_LAYOUT=$2;; --cleandata-dir) CLEAN_DIR=$2;; --clean-layout) CLEAN_LAYOUT=$2;; --read-type) READ_TYPE=$2;; --assembly-dir) ASSEMBLY_DIR=$2;; --manifest) MANIFEST=$2;;
      esac; shift 2;;
    -h|--help) usage; exit 0;; *) die 2 "Unknown option: $1";;
  esac
done
[[ $TASK == qc_only || $TASK == assembly_only || $TASK == full ]] || die 2 'Invalid --task'
[[ -n $CLEAN_DIR && -n $MANIFEST ]] || die 2 '--cleandata-dir and --manifest are required'
require_command find; require_command realpath
require_allowed_path 'cleandata directory' "$CLEAN_DIR"; assert_output_parent 'cleandata directory' "$CLEAN_DIR"
if [[ $TASK == qc_only || $TASK == full ]]; then
  [[ -n $RAW_DIR && ( $RAW_LAYOUT == sample_subdirs || $RAW_LAYOUT == flat ) ]] || die 2 'rawdata arguments are required for QC'
  require_command fastp; require_allowed_path 'rawdata directory' "$RAW_DIR"; assert_existing_dir 'rawdata directory' "$RAW_DIR"
fi
if [[ $TASK == assembly_only || $TASK == full ]]; then
  [[ -n $ASSEMBLY_DIR ]] || die 2 '--assembly-dir is required for assembly'
  require_command megahit; require_allowed_path 'assembly directory' "$ASSEMBLY_DIR"; assert_output_parent 'assembly directory' "$ASSEMBLY_DIR"
  [[ $TASK != full ]] || { CLEAN_LAYOUT=sample_subdirs; READ_TYPE=pe; }
  [[ $CLEAN_LAYOUT == sample_subdirs || $CLEAN_LAYOUT == flat ]] || die 2 'Invalid --clean-layout'
  [[ $READ_TYPE == pe || $READ_TYPE == se ]] || die 2 'Invalid --read-type'
fi

mkdir -p "$(dirname "$MANIFEST")"; TMP_MANIFEST="${MANIFEST}.tmp.$$"; trap 'rm -f "$TMP_MANIFEST"' EXIT
manifest_header > "$TMP_MANIFEST"; declare -A SEEN=()
add_pe() {
  local sample=$1 raw1=$2 raw2=$3 c1 c2 asm=''
  valid_sample_id "$sample" || die 3 "Invalid sample ID: $sample"; [[ -z ${SEEN[$sample]+x} ]] || die 3 "Duplicate sample ID: $sample"; SEEN[$sample]=1
  c1="$CLEAN_DIR/$sample/${sample}_R1.clean.fq.gz"; c2="$CLEAN_DIR/$sample/${sample}_R2.clean.fq.gz"; [[ -z $ASSEMBLY_DIR ]] || asm="$ASSEMBLY_DIR/$sample"
  manifest_add_pe "$TMP_MANIFEST" "$sample" "$raw1" "$raw2" "$c1" "$c2" "$asm"
}
add_se() {
  local sample=$1 clean=$2
  valid_sample_id "$sample" || die 3 "Invalid sample ID: $sample"; [[ -z ${SEEN[$sample]+x} ]] || die 3 "Duplicate sample ID: $sample"; SEEN[$sample]=1
  manifest_add_se "$TMP_MANIFEST" "$sample" "$clean" "$ASSEMBLY_DIR/$sample"
}
scan_flat_pe() {
  local source=$1 mode=$2 file name stem i; declare -A r1=() r2=()
  while IFS= read -r -d '' file; do
    name=${file##*/}
    for ((i=0; i<${#R1_SUFFIXES[@]}; i++)); do
      if [[ $name == *"${R1_SUFFIXES[$i]}" ]]; then stem=${name%"${R1_SUFFIXES[$i]}"}; [[ -z ${r1[$stem]+x} ]] || die 3 "Duplicate R1 for $stem"; r1[$stem]=$file; fi
      if [[ $name == *"${R2_SUFFIXES[$i]}" ]]; then stem=${name%"${R2_SUFFIXES[$i]}"}; [[ -z ${r2[$stem]+x} ]] || die 3 "Duplicate R2 for $stem"; r2[$stem]=$file; fi
    done
  done < <(find "$source" -maxdepth 1 -type f \( -iname '*.fq.gz' -o -iname '*.fastq.gz' \) -print0)
  (( ${#r1[@]} > 0 )) || die 3 "No recognised R1 FASTQ files found: $source"
  for stem in "${!r1[@]}"; do
    [[ -n ${r2[$stem]+x} ]] || die 3 "$stem: R2 is missing"
    if [[ $mode == raw ]]; then
      add_pe "$stem" "${r1[$stem]}" "${r2[$stem]}"
    else
      valid_sample_id "$stem" || die 3 "Invalid sample ID: $stem"
      [[ -z ${SEEN[$stem]+x} ]] || die 3 "Duplicate sample ID: $stem"; SEEN[$stem]=1
      manifest_add_pe "$TMP_MANIFEST" "$stem" '' '' "${r1[$stem]}" "${r2[$stem]}" "$ASSEMBLY_DIR/$stem"
    fi
  done
  for stem in "${!r2[@]}"; do [[ -n ${r1[$stem]+x} ]] || die 3 "$stem: R1 is missing"; done
}

if [[ $TASK == qc_only || $TASK == full ]]; then
  if [[ $RAW_LAYOUT == sample_subdirs ]]; then
    found=0
    while IFS= read -r -d '' dir; do found=1; sample=${dir##*/}; resolve_unique_pe_pair_in_dir "$dir" || die 3 "$sample: expected exactly one recognised R1/R2 pair"; add_pe "$sample" "$RESOLVED_R1" "$RESOLVED_R2"; done < <(find "$RAW_DIR" -mindepth 1 -maxdepth 1 -type d -print0)
    (( found )) || die 3 "No sample subdirectories found: $RAW_DIR"
  else scan_flat_pe "$RAW_DIR" raw; fi
else
  assert_existing_dir 'cleandata directory' "$CLEAN_DIR"
  if [[ $READ_TYPE == pe ]]; then
    if [[ $CLEAN_LAYOUT == sample_subdirs ]]; then
      found=0
      while IFS= read -r -d '' dir; do
        found=1; sample=${dir##*/}; resolve_unique_pe_pair_in_dir "$dir" || die 3 "$sample: expected exactly one PE pair"
        valid_sample_id "$sample" || die 3 "Invalid sample ID: $sample"; [[ -z ${SEEN[$sample]+x} ]] || die 3 "Duplicate sample ID: $sample"; SEEN[$sample]=1
        manifest_add_pe "$TMP_MANIFEST" "$sample" '' '' "$RESOLVED_R1" "$RESOLVED_R2" "$ASSEMBLY_DIR/$sample"
      done < <(find "$CLEAN_DIR" -mindepth 1 -maxdepth 1 -type d -print0)
      (( found )) || die 3 'No clean-data sample subdirectories found'
    else scan_flat_pe "$CLEAN_DIR" clean; fi
  elif [[ $CLEAN_LAYOUT == sample_subdirs ]]; then
    found=0
    while IFS= read -r -d '' dir; do found=1; sample=${dir##*/}; resolve_unique_se_file_in_dir "$dir" || die 3 "$sample: expected exactly one SE FASTQ"; add_se "$sample" "$RESOLVED_SE"; done < <(find "$CLEAN_DIR" -mindepth 1 -maxdepth 1 -type d -print0)
    (( found )) || die 3 'No clean-data sample subdirectories found'
  else
    while IFS= read -r -d '' file; do
      name=${file##*/}; case "$name" in *.clean.fq.gz) sample=${name%.clean.fq.gz};; *.clean.fastq.gz) sample=${name%.clean.fastq.gz};; *.fq.gz) sample=${name%.fq.gz};; *.fastq.gz) sample=${name%.fastq.gz};; *) continue;; esac; add_se "$sample" "$file"
    done < <(find "$CLEAN_DIR" -maxdepth 1 -type f \( -iname '*.fq.gz' -o -iname '*.fastq.gz' \) -print0)
    (( ${#SEEN[@]} > 0 )) || die 3 'No SE clean reads found'
  fi
fi
mv -f -- "$TMP_MANIFEST" "$MANIFEST"; trap - EXIT
printf '[INFO] Preflight passed: %d sample(s). Manifest: %s\n' "$(( $(wc -l < "$MANIFEST") - 1 ))" "$MANIFEST"
