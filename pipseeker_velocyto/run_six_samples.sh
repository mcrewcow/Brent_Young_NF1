#!/usr/bin/env bash
# PIPseeker v5 -> provisional velocyto Loom. No FASTQs needed.
# Original PIPseeker files are NEVER modified or removed.
set -Eeuo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIP_ROOT="${PIP_ROOT:-/mnt/c/Users/rodri/Brent_collab/pipseeker_runs}"
RESULT_ROOT="${RESULT_ROOT:-${PIP_ROOT}/velocyto_results}"
WORK_ROOT="${WORK_ROOT:-${HOME}/projects/Brent_velocyto/batch_work}"
GTF_ORIGINAL="${GTF_ORIGINAL:-${HOME}/references/GRCm39_Cre_source/GRCm39_Cre.gtf}"
THREADS="${THREADS:-8}"
MEM_MB_PER_THREAD="${MEM_MB_PER_THREAD:-1024}"
MASK="${REPEAT_MASK:-}"
SAMPLES="1,2,3,4,5,6"
DRY_RUN=0
ACCEPT=0
KEEP=0
REUSE_LEGACY_SAMPLE1=0
IMPORT_LEGACY_SAMPLE1_LOOM=0

usage() {
  cat <<'EOF'
Usage: bash run_six_samples.sh [options]

  --dry-run                     Preview paths, files and actions (no writes)
  --accept-exploratory-imi      REQUIRED: acknowledge v5 IMI counts are provisional
  --samples 1,2,3               Select samples (default: 1,2,3,4,5,6)
  --reuse-legacy-sample1        Reuse existing sample1 full tagged BAM if present
  --import-legacy-sample1-loom   If corrected sample1 Loom exists, validate/import it
  --keep-intermediates          Do not remove derived BAMs after validated Loom
  --threads N                   samtools sorting threads (default 8)
  --samtools-memory N           MB per thread for samtools sort (default 1024)
  --mask /path/repeats.gtf      Optional appropriate genome repeat mask
  --help                        Show this message

Configurable environment variables:
  PIP_ROOT, RESULT_ROOT, WORK_ROOT, GTF_ORIGINAL, REPEAT_MASK

Defaults:
  Inputs:  /mnt/c/Users/rodri/Brent_collab/pipseeker_runs/sampleN_GRCm39_Cre
  Results: /mnt/c/Users/rodri/Brent_collab/pipseeker_runs/velocyto_results
  Scratch: ~/projects/Brent_velocyto/batch_work (Linux filesystem)

Safeguards:
  - Existing PIPseeker input folders are READ ONLY.
  - Requires exact v5 chemistry and GRCm39_2024A_Cre reference in run_config.csv.
  - Validates barcode encoding before writing any BAM.
  - Skips samples with completed QC and final Loom; never overwrites silently.
  - Removes only pipeline-created derived BAMs after validating and copying Loom.
  - PIPseq v5 molecular-count amplification correction is NOT reproduced.
EOF
}

while (($#)); do
  case "$1" in
    --dry-run) DRY_RUN=1; shift;;
    --accept-exploratory-imi) ACCEPT=1; shift;;
    --samples) SAMPLES="${2:?missing samples}"; shift 2;;
    --reuse-legacy-sample1) REUSE_LEGACY_SAMPLE1=1; shift;;
    --import-legacy-sample1-loom) IMPORT_LEGACY_SAMPLE1_LOOM=1; shift;;
    --keep-intermediates) KEEP=1; shift;;
    --threads) THREADS="${2:?missing threads}"; shift 2;;
    --samtools-memory) MEM_MB_PER_THREAD="${2:?missing memory}"; shift 2;;
    --mask) MASK="${2:?missing repeat mask}"; shift 2;;
    --help|-h) usage; exit 0;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2;;
  esac
done

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo 'Invalid threads' >&2; exit 2; }
[[ "$MEM_MB_PER_THREAD" =~ ^[1-9][0-9]*$ ]] || { echo 'Invalid memory' >&2; exit 2; }
IFS=',' read -r -a IDS <<< "$SAMPLES"
for id in "${IDS[@]}"; do
  [[ "$id" =~ ^[1-6]$ ]] || { echo "Invalid sample ID: $id" >&2; exit 2; }
done

printf '\nPIP input root: %s\nResults root:   %s\nLinux work dir: %s\nOriginal GTF:   %s\nSamples:        %s\n' \
  "$PIP_ROOT" "$RESULT_ROOT" "$WORK_ROOT" "$GTF_ORIGINAL" "$SAMPLES"
if [[ -n "$MASK" ]]; then echo "Repeat mask:    $MASK"; else echo 'Repeat mask:    NOT SUPPLIED (exploratory)'; fi
if (( DRY_RUN )); then
  echo -e '\n[DRY RUN: no files will be created or deleted]'
  for id in "${IDS[@]}"; do
    name="sample${id}_GRCm39_Cre"
    input="${PIP_ROOT}/${name}"
    echo "--- $name ---"
    for f in "${input}/star_out.bam" "${input}/run_config.csv" \
      "${input}/barcodes/barcode_whitelist.txt" \
      "${input}/filtered_matrix/sensitivity_5/barcodes.tsv.gz"; do
      if [[ -f "$f" ]]; then echo "  [OK] $f"; else echo "  [MISSING] $f"; fi
    done
    if [[ -f "${RESULT_ROOT}/${name}/${name}.loom" ]]; then echo '  Existing result detected; check its QC before re-running'; fi
  done
  exit 0
fi
(( ACCEPT )) || { echo 'Refusing to run: add --accept-exploratory-imi after reading README limitations' >&2; exit 2; }
for cmd in python samtools velocyto awk df; do command -v "$cmd" >/dev/null || { echo "Missing $cmd; activate velocyto-py310" >&2; exit 1; }; done
python -c 'import pysam,h5py,numpy,scipy,velocyto' >/dev/null
[[ -r "$GTF_ORIGINAL" ]] || { echo "Missing original GTF: $GTF_ORIGINAL" >&2; exit 1; }
[[ -z "$MASK" || -r "$MASK" ]] || { echo "Missing mask: $MASK" >&2; exit 1; }
mkdir -p "$WORK_ROOT/reference" "$RESULT_ROOT"

# Proven correction for velocyto 0.17.17: internally it normalizes BAM chrM to MT,
# but normalizes GTF chrM to M. Rename only chrM in a *copy* to chrMT.
GTF_FIXED="${WORK_ROOT}/reference/GRCm39_Cre.velocyto_chrMT.gtf"
if [[ ! -s "$GTF_FIXED" || "$GTF_ORIGINAL" -nt "$GTF_FIXED" ]]; then
  awk 'BEGIN{FS=OFS="\t"} /^#/ {print;next} $1=="chrM" {$1="chrMT"} {print}' "$GTF_ORIGINAL" > "${GTF_FIXED}.partial"
  mv -f "${GTF_FIXED}.partial" "$GTF_FIXED"
fi
MT_EXONS="$(awk -F'\t' '$1=="chrMT" && $3=="exon" {n++} END {print n+0}' "$GTF_FIXED")"
[[ "$MT_EXONS" -gt 0 ]] || { echo 'Fixed GTF has zero chrMT exons; check genome naming!' >&2; exit 1; }
echo "Corrected GTF: $GTF_FIXED ($MT_EXONS mitochondrial exons)"

for id in "${IDS[@]}"; do
  SAMPLE="sample${id}_GRCm39_Cre"
  SAMPLE_INPUT="${PIP_ROOT}/${SAMPLE}"
  SAMPLE_WORK="${WORK_ROOT}/${SAMPLE}"
  SAMPLE_RESULT="${RESULT_ROOT}/${SAMPLE}"
  FINAL_LOOM="${SAMPLE_RESULT}/${SAMPLE}.loom"
  QC_FINAL="${SAMPLE_RESULT}/qc_summary.json"
  DONE="${SAMPLE_RESULT}/.complete"
  SOURCE="${SAMPLE_INPUT}/star_out.bam"
  echo -e "\n===== ${SAMPLE} ====="

  if [[ -s "$FINAL_LOOM" && -s "$QC_FINAL" && -f "$DONE" ]]; then
    echo "Skipping completed sample: $FINAL_LOOM"
    continue
  fi
  if [[ -e "$FINAL_LOOM" || -e "$QC_FINAL" || -e "$DONE" ]]; then
    echo "Partial existing output in $SAMPLE_RESULT; inspect and move it before retrying" >&2
    exit 1
  fi
  for f in "$SOURCE" "${SAMPLE_INPUT}/run_config.csv" \
     "${SAMPLE_INPUT}/barcodes/barcode_whitelist.txt" \
     "${SAMPLE_INPUT}/filtered_matrix/sensitivity_5/barcodes.tsv.gz" \
     "${SAMPLE_INPUT}/filtered_matrix/sensitivity_5/matrix.mtx.gz"; do
    [[ -r "$f" ]] || { echo "Missing input $f" >&2; exit 1; }
  done
  # Prevent processing a sample with different chemistry or reference.
  python - "$SAMPLE_INPUT/run_config.csv" <<'PY'
import csv,sys
with open(sys.argv[1], newline='') as f:
    d={r[0]:r[1] for r in csv.reader(f) if len(r)>=2}
if d.get('chemistry','').strip().upper()!='V':
    raise SystemExit(f'Unexpected chemistry: {d.get("chemistry")!r}; expected V')
if 'GRCm39_2024A_Cre' not in d.get('star_index_path',''):
    raise SystemExit('Unexpected reference: ' + d.get('star_index_path',''))
print('PIPseeker chemistry:',d['chemistry'],'| reference:',d['star_index_path'])
PY

  mkdir -p "$SAMPLE_WORK" "$SAMPLE_RESULT/logs"
  # Optional zero-reprocessing import of corrected sample-1 Loom produced earlier.
  LEGACY_LOOM="${HOME}/projects/Brent_velocyto/full/loom_chrMT/sample1_GRCm39_Cre_sens5_chrMT.loom"
  if [[ "$id" == 1 && "$IMPORT_LEGACY_SAMPLE1_LOOM" == 1 && -s "$LEGACY_LOOM" ]]; then
    echo 'Validating/importing the previously corrected sample1 chrMT Loom...'
    python "$HERE/qc_loom.py" --loom "$LEGACY_LOOM" \
      --sample-dir "$SAMPLE_INPUT" --sample "$SAMPLE" \
      --out "${SAMPLE_WORK}/qc_summary.json" > "${SAMPLE_WORK}/qc.log"
    python - "${SAMPLE_WORK}/qc_summary.json" <<'PY_QC'
import json,sys
with open(sys.argv[1]) as f: report=json.load(f)
assert report['mt_total_counts']>0, 'Legacy Loom still has zero mitochondrial counts'
PY_QC
    cp "$LEGACY_LOOM" "${FINAL_LOOM}.partial"
    cmp -s "$LEGACY_LOOM" "${FINAL_LOOM}.partial" || { echo 'Failed to verify copied loom' >&2; exit 1; }
    mv "${FINAL_LOOM}.partial" "$FINAL_LOOM"
    cp "${SAMPLE_WORK}/qc_summary.json" "$QC_FINAL"
    cp "${SAMPLE_WORK}/qc.log" "${SAMPLE_RESULT}/logs/qc.log"
    date -Is > "$DONE"
    python "$HERE/collect_summary.py" --results "$RESULT_ROOT"
    echo "Imported verified sample1 Loom -> $FINAL_LOOM"
    continue
  fi
  TAGGED="${SAMPLE_WORK}/${SAMPLE}.unique_sens5.possorted.bam"
  CELLSORTED="${SAMPLE_WORK}/cellsorted_${SAMPLE}.unique_sens5.possorted.bam"
  SORT_TMP="${SAMPLE_WORK}/sort_tmp"
  REUSED=0

  if [[ ! -s "$TAGGED" ]]; then
    if [[ "$id" == 1 && "$REUSE_LEGACY_SAMPLE1" == 1 && \
       -s "${HOME}/projects/Brent_velocyto/full/sample1_sens5_unique.possorted.bam" ]]; then
      ln -s "${HOME}/projects/Brent_velocyto/full/sample1_sens5_unique.possorted.bam" "$TAGGED"
      REUSED=1
      echo 'Reusing the previously validated sample1 tagged BAM via symbolic link'
    else
      # Sorting by cell requires more scratch space later; estimate a conservative reserve.
      INPUT_KB=$(( ($(stat -c '%s' "$SOURCE") + 1023) / 1024 ))
      REQUIRED_KB=$(( INPUT_KB * 3 ))
      (( REQUIRED_KB < 60*1024*1024 )) && REQUIRED_KB=$((60*1024*1024))
      FREE_KB="$(df -Pk "$WORK_ROOT" | awk 'END {print $4}')"
      if (( FREE_KB < REQUIRED_KB )); then
        echo "Insufficient Linux scratch space (available ${FREE_KB}KB, estimate ${REQUIRED_KB}KB)" >&2
        exit 1
      fi
      echo 'Converting QNAME -> CB/UB/NH; retaining unique primary alignments; coordinate-sorting...'
      # Never stream text to stdout of the converter (it must remain BAM only).
      python -u "$HERE/reconstruct_bam.py" \
        --sample-dir "$SAMPLE_INPUT" \
        --summary "${SAMPLE_WORK}/conversion_stats.json" \
        2> "${SAMPLE_WORK}/conversion.log" \
        | samtools sort -@ "$THREADS" -m "${MEM_MB_PER_THREAD}M" \
            -T "$SORT_TMP" -O BAM \
            -o "${TAGGED}.partial.bam" - \
            2> "${SAMPLE_WORK}/samtools_sort.log"
      samtools quickcheck -v "${TAGGED}.partial.bam"
      mv -f "${TAGGED}.partial.bam" "$TAGGED"
    fi
  else
    echo 'Reusing existing in-work tagged BAM'
  fi
  samtools quickcheck -v "$TAGGED"
  if [[ ! -f "${TAGGED}.bai" ]]; then
    samtools index -@ "$THREADS" "$TAGGED"
  fi

  LOCAL_LOOM_DIR="${SAMPLE_WORK}/loom"
  mkdir -p "$LOCAL_LOOM_DIR"
  LOCAL_LOOM="${LOCAL_LOOM_DIR}/${SAMPLE}.loom"
  if [[ ! -s "$LOCAL_LOOM" ]]; then
    echo 'Running velocyto...'
    OPTIONS=(-b "${SAMPLE_INPUT}/filtered_matrix/sensitivity_5/barcodes.tsv.gz"
             -o "$LOCAL_LOOM_DIR" -e "$SAMPLE" -t uint32
             --samtools-threads "$THREADS" --samtools-memory "$MEM_MB_PER_THREAD")
    if [[ -n "$MASK" ]]; then OPTIONS+=(-m "$MASK"); fi
    velocyto run "${OPTIONS[@]}" -vv "$TAGGED" "$GTF_FIXED" \
      2>&1 | tee "${SAMPLE_WORK}/velocyto.log"
  fi
  [[ -s "$LOCAL_LOOM" ]] || { echo 'velocyto did not create Loom' >&2; exit 1; }

  echo 'Validating recovered cell and gene IDs and computing count discrepancy QC...'
  python "$HERE/qc_loom.py" --loom "$LOCAL_LOOM" \
    --sample-dir "$SAMPLE_INPUT" --sample "$SAMPLE" \
    --out "${SAMPLE_WORK}/qc_summary.json" \
    > "${SAMPLE_WORK}/qc.log"

  # Result directory on Windows; copy and verify before marking complete.
  cp "$LOCAL_LOOM" "${FINAL_LOOM}.partial"
  cmp -s "$LOCAL_LOOM" "${FINAL_LOOM}.partial" || { echo 'Copied loom differs from source' >&2; exit 1; }
  mv "${FINAL_LOOM}.partial" "$FINAL_LOOM"
  cp "${SAMPLE_WORK}/qc_summary.json" "$QC_FINAL"
  for log in conversion.log samtools_sort.log velocyto.log qc.log conversion_stats.json; do
    if [[ -f "${SAMPLE_WORK}/${log}" ]]; then cp "${SAMPLE_WORK}/${log}" "${SAMPLE_RESULT}/logs/${log}"; fi
  done
  date -Is > "$DONE"
  echo "FINISHED ${SAMPLE} -> ${FINAL_LOOM}"

  if (( ! KEEP )); then
    # Remove ONLY generated in-work artifacts. Never the original PIPseeker outputs
    # nor the preexisting sample1 BAM that may have been symlinked above.
    rm -f -- "$TAGGED" "${TAGGED}.bai" "$CELLSORTED" "${CELLSORTED}.bai" \
      "${SAMPLE_WORK}/conversion_stats.json" "$LOCAL_LOOM"
    echo 'Removed in-work BAMs and temporary Loom after validated Windows copy'
  fi
  python "$HERE/collect_summary.py" --results "$RESULT_ROOT"
done

echo -e "\nAll requested samples completed or previously validated. Results: $RESULT_ROOT"
python "$HERE/collect_summary.py" --results "$RESULT_ROOT"
