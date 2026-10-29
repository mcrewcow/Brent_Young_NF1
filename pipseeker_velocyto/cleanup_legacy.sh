#!/usr/bin/env bash
# CONSERVATIVE cleanup of files generated during earlier sample1 testing.
# Dry run by default. Never touches PIPseeker source or output results.
set -Eeuo pipefail
DO_DELETE=0
PURGE_TAGGED=0
PURGE_CORRECTED_LEGACY=0
RESULT_ROOT="${RESULT_ROOT:-/mnt/c/Users/rodri/Brent_collab/pipseeker_runs/velocyto_results}"
LEGACY="${HOME}/projects/Brent_velocyto"
for arg in "$@"; do
  case "$arg" in
    --execute) DO_DELETE=1;;
    --purge-legacy-tagged-bam) PURGE_TAGGED=1;;
    --purge-legacy-corrected-loom) PURGE_CORRECTED_LEGACY=1;;
    --help|-h)
      cat <<'EOF'
Usage: bash cleanup_legacy.sh [--execute] [--purge-legacy-tagged-bam] [--purge-legacy-corrected-loom]
Dry-run by default. Requires validated sample1 result copied to Windows results first.
Removes only exact known pilot, pre-correction and sorting files in ~/projects/Brent_velocyto.
--purge-legacy-tagged-bam removes the large old full tagged BAM; NOT needed if
  you still plan to use --reuse-legacy-sample1.
--purge-legacy-corrected-loom deletes the redundant old chrMT Loom only after it
  has been successfully imported into the Windows results directory.
Never deletes inputs under /mnt/c/.../pipseeker_runs/sampleN_GRCm39_Cre.
EOF
      exit 0;;
    *) echo "Unknown option: $arg" >&2; exit 2;;
  esac
done
SAMPLE1_RESULT="$RESULT_ROOT/sample1_GRCm39_Cre"
if [[ ! -s "$SAMPLE1_RESULT/sample1_GRCm39_Cre.loom" || \
      ! -s "$SAMPLE1_RESULT/qc_summary.json" || ! -f "$SAMPLE1_RESULT/.complete" ]]; then
  echo 'Refusing cleanup: no completed, validated sample1 Windows result found.' >&2
  exit 1
fi
FILES=(
  "$LEGACY/test/sample1_test.unsorted.bam"
  "$LEGACY/test/sample1_test.possorted.bam"
  "$LEGACY/test/sample1_test.possorted.bam.bai"
  "$LEGACY/test/cellsorted_sample1_test.possorted.bam"
  "$LEGACY/test/velocyto_out/sample1_velocity_pilot.loom"
  "$LEGACY/test/velocyto_chrMT_test/sample1_chrMT_pilot.loom"
  "$LEGACY/full/loom/sample1_GRCm39_Cre_sens5.loom"
  "$LEGACY/full/cellsorted_sample1_sens5_unique.possorted.bam"
  "$LEGACY/full/cellsorted_sample1_sens5_unique.possorted.bam.bai"
)
if (( PURGE_TAGGED )); then
  FILES+=("$LEGACY/full/sample1_sens5_unique.possorted.bam"
          "$LEGACY/full/sample1_sens5_unique.possorted.bam.bai")
fi
if (( PURGE_CORRECTED_LEGACY )); then
  FILES+=("$LEGACY/full/loom_chrMT/sample1_GRCm39_Cre_sens5_chrMT.loom")
fi
printf 'Mode: %s\n' "$([[ "$DO_DELETE" == 1 ]] && echo 'EXECUTE' || echo 'PREVIEW ONLY')"
for p in "${FILES[@]}"; do
  if [[ -f "$p" || -L "$p" ]]; then
    if (( DO_DELETE )); then rm -f -- "$p"; echo "DELETED $p";
    else ls -lh -- "$p"; fi
  fi
done
if (( ! DO_DELETE )); then
  echo 'This was a preview. Run with --execute only after reviewing listed files.'
fi
