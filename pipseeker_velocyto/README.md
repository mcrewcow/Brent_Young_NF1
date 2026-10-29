# PIPseeker v5 → velocyto (experimental reconstruction)

**Purpose:** Reuse previously generated PIPseeker 3.3 STAR BAMs for six mouse GRCm39_Cre samples when the original FASTQs are no longer locally available. Creates preliminary `spliced`, `unspliced`, and `ambiguous` Loom layers, retaining PIPseeker sensitivity-5 cell barcodes.

> **Scientific caveat (important):** This pipeline does **not** implement PIPseeker v5's IMI/bin amplification correction or reproduce molecule quantification. The output is **exploratory**, not a validated publication-quality v5 RNA-velocity quantification. Gene/cell expression correlation is a consistency check, **not** a validation of splicing kinetics or velocity direction. Do not blindly rescale counts by a global factor. Reference: [PIPseeker User Guide, v3.2](https://fbs-public.s3.us-east-2.amazonaws.com/public-pipseeker-releases/pipseeker-v3.2.0/PIPseeker-v3.2-User-Guide.pdf) and [velocyto generic BAM workflow](https://velocyto.org/velocyto.py/tutorial/cli.html).

## Tested starting point

Successfully piloted on `sample1_GRCm39_Cre` (chemistry `V`, PIPseeker 3.3, GRCm39_2024A_Cre reference):

- Read-name component 2 encodes 16-base cell barcode as forward 2-bit `ACGT` packed integer; barcode mapping matched all 100,000 mapped reads sampled against whitelist.
- Read-name component 3 is treated as a **candidate** 30-bit packed molecular identifier. Encoding is bijective, but deduplication is **not** equivalent to PIPseeker's v5 chemistry-aware counting.
- BAM lacks auxiliary tags; STAR MAPQ 255 primary alignments are provisionally assigned `NH:i:1`; `CB` and `UB` tags are reconstructed.
- Original GTF built the PIPseeker STAR index: `~/references/GRCm39_Cre_source/GRCm39_Cre.gtf`.
- velocyto 0.17.17 internally normalizes `chrM` differently between BAM and GTF. Pilot restored mitochondrial counts after replacing GTF `chrM` with `chrMT` in a **copy**, without changing the original input files.
- Full sample1 before corrected mitochondrial re-run had 2,012 cells, 33,697 genes, gene-level Spearman 0.9497 and cell-level Spearman 0.9922 against original counts, but 2.3546× overall counts (not resolved by a successful Loom run).

**Do not automatically apply this reconstruction to other PIPseeker chemistries or reference builds.** Each sample must have `run_config.csv` showing chemistry `V` and the expected GRCm39_2024A_Cre STAR reference. The script stops on mismatches. The preflight checks read-name encoding against the corresponding whitelist for each sample.

## Prerequisites (WSL2)

Run using the working Conda environment:

```bash
conda activate velocyto-py310
python -c 'import pysam,h5py,numpy,scipy,velocyto; print("Imports OK")'
samtools --version | head -n 2
velocyto --version
```

If you need a new environment instead, an indicative environment specification is included as `environment.yml`. Since conda dependencies can change, test a clean environment before relying on it. Do not modify your working `lab-py312` environment to install velocyto.

## Files and locations

Six source folders (read-only):

```text
/mnt/c/Users/rodri/Brent_collab/pipseeker_runs/
  sample1_GRCm39_Cre/
  sample2_GRCm39_Cre/
  sample3_GRCm39_Cre/
  sample4_GRCm39_Cre/
  sample5_GRCm39_Cre/
  sample6_GRCm39_Cre/
```

For each sample, the pipeline needs `star_out.bam`, `barcodes/barcode_whitelist.txt`, `run_config.csv`, and `filtered_matrix/sensitivity_5/` with barcode, feature, and matrix files. The FASTQs are **not** required.

Outputs (Windows-accessible):

```text
/mnt/c/Users/rodri/Brent_collab/pipseeker_runs/velocyto_results/
  sample1_GRCm39_Cre/
    sample1_GRCm39_Cre.loom
    qc_summary.json
    .complete
    logs/
  ...
  sample6_GRCm39_Cre/
    sample6_GRCm39_Cre.loom
    qc_summary.json
    .complete
    logs/
  all_samples_qc.csv
```

Scratch/intermediates (on faster Linux filesystem):

```text
~/projects/Brent_velocyto/batch_work/
```

By default, intermediate tagged and cell-sorted BAMs are removed **after** a complete, validated Loom is copied and byte-compared to the Windows results location. On failure, artifacts and logs stay available for debugging. No files inside the six original PIPseeker sample folders are deleted or modified.

## Installation from ZIP into WSL

Download the ZIP archive provided with this project. From WSL (if downloaded to Windows Downloads):

```bash
mkdir -p "$HOME/projects/Brent_velocyto"
unzip '/mnt/c/Users/rodri/Downloads/pipseeker_velocyto_pipeline.zip' \
  -d "$HOME/projects/Brent_velocyto"
cd "$HOME/projects/Brent_velocyto/pipseeker_velocyto"
```

## Running the pipeline

**1. Preview:**

```bash
conda activate velocyto-py310
bash run_six_samples.sh --dry-run
```

**2. Run all six:**

```bash
bash run_six_samples.sh --accept-exploratory-imi
```

**Recommended to avoid rebuilding sample1:** If sample1's **corrected** full Loom has already been produced in `~/projects/Brent_velocyto/full/loom_chrMT/`, import and independently QC it:

```bash
bash run_six_samples.sh --import-legacy-sample1-loom --accept-exploratory-imi
```

If the corrected sample1 Loom has **not** been produced but its previously reconstructed coordinate-sorted tagged BAM is still present:

```bash
bash run_six_samples.sh --reuse-legacy-sample1 --accept-exploratory-imi
```

Only one shortcut is needed. The import option takes priority. The script otherwise rebuilds or processes samples 1–6 normally. Existing completed results are skipped, allowing interrupted batch processing to resume safely. If files exist but the `.complete` marker is absent, the pipeline stops instead of overwriting them; inspect those outputs first.

**3. Optional subsets and more memory control:**

```bash
bash run_six_samples.sh --samples 2,3,4 --threads 8 --samtools-memory 1024 \
  --accept-exploratory-imi
```

Optional repeat mask:

```bash
bash run_six_samples.sh --mask /path/to/compatible_mouse_repeat_mask.gtf \
  --accept-exploratory-imi
```

Only use an annotation/repeat mask built for the same reference assembly and contig conventions. No repeat mask is applied by default. Reprocessing with a mask requires a separately archived output or manual moving of existing results (the pipeline never silently overwrites completed results).

Paths can be overridden with environment variables:

```bash
PIP_ROOT=/mnt/c/.../pipseeker_runs \
RESULT_ROOT=/mnt/c/.../velocyto_results \
WORK_ROOT="$HOME/projects/Brent_velocyto/batch_work" \
GTF_ORIGINAL="$HOME/references/GRCm39_Cre_source/GRCm39_Cre.gtf" \
bash run_six_samples.sh --accept-exploratory-imi
```

## Cleanup

The normal batch run automatically deletes its successful per-sample *temporary* BAMs and Linux Loom copies, keeping final Windows Loom files, QC summaries and logs.

For old pilot intermediates produced before this repository, run:

```bash
bash cleanup_legacy.sh                     # PREVIEW only; deletes nothing
bash cleanup_legacy.sh --execute           # delete only known pilot / old interim files
```

The cleanup script **refuses to operate** until it finds a completed and validated sample1 Loom in the Windows results location. It does not touch the original PIPseeker data. It also preserves the old full sample1 reconstructed BAM and corrected Loom unless those are explicitly included with `--purge-legacy-tagged-bam` / `--purge-legacy-corrected-loom`. Only purge them when you no longer need the reuse options and have verified the final copy.

## QC for later six-sample discrepancy review

The batch automatically writes `all_samples_qc.csv` with all samples' original and velocyto counts, global ratios, median and percentile per-cell ratios, gene- and cell-level Spearman correlations, and spliced/unspliced/ambiguous and mitochondrial counts.

After all six, investigate:

1. Count inflation dependence on depth, cell calling threshold, and gene identity, including the global and per-cell ratios. No universal multiplicative correction.
2. Candidate IMI/bin semantics and collision/duplicate handling versus the PIPseeker v5 `molecule_info.h5` records and official molecular counting rules.
3. Nuclear splicing-classification behavior, intron annotation, and sensitivity to repeat masking.
4. Whether velocity directions are stable under reasonable counting/quality-control choices. Until then, mark the Loom outputs *exploratory*.

## Uploading to GitHub

GitHub should contain only these source files, **not the BAM/FASTQ/LOOM/MTX datasets or log outputs**. The `.gitignore` blocks typical generated data files.

```bash
git init
git add README.md run_six_samples.sh reconstruct_bam.py qc_loom.py \
  collect_summary.py cleanup_legacy.sh environment.yml .gitignore
git commit -m 'Add PIPseeker v5 exploratory velocyto reconstruction pipeline'
```

You can then push to a **private** repository through your preferred GitHub method. No original sequencing data is included in this code archive.
