#!/usr/bin/env python3
"""Collect per-sample QC JSON files into a compact CSV."""
import argparse
import csv
import json
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument('--results', type=Path, required=True)
args = ap.parse_args()
rows = []
for p in sorted(args.results.glob('sample*_GRCm39_Cre/qc_summary.json')):
    with p.open() as f:
        d = json.load(f)
    rows.append({
        'sample': d['sample'], 'cells': d['cells'], 'genes': d['genes'],
        'pipseeker_total': d['original_total_counts'],
        'velocyto_total': d['velocyto_total_counts'],
        'ratio_global': d['global_count_ratio'],
        'ratio_median_cell': d['per_cell_count_ratio_median'],
        'ratio_p5_cell': d['per_cell_count_ratio_p5'],
        'ratio_p95_cell': d['per_cell_count_ratio_p95'],
        'gene_spearman': d['gene_spearman_orig_ge10'],
        'cell_spearman': d['cell_spearman'],
        'spliced': d['layer_totals']['spliced'],
        'unspliced': d['layer_totals']['unspliced'],
        'ambiguous': d['layer_totals']['ambiguous'],
        'mt_total': d['mt_total_counts'],
    })
out = args.results / 'all_samples_qc.csv'
if rows:
    with out.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f'Updated {out} ({len(rows)} samples)')
else:
    print('No per-sample reports yet')
