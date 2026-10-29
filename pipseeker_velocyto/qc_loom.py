#!/usr/bin/env python3
"""Validate full velocyto loom against original sensitivity-5 PIPseeker matrix.

Reports agreement; does NOT certify RNA velocity validity or v5 IMI counting.
"""
import argparse
import gzip
import json
import re
from pathlib import Path

import h5py
import numpy as np
from scipy.io import mmread
from scipy.stats import spearmanr


def find(folder, names):
    for n in names:
        p = folder / n
        if p.is_file():
            return p
    raise FileNotFoundError(f'Cannot find {names} under {folder}')


def lines(path):
    op = gzip.open if path.suffix == '.gz' else open
    with op(path, 'rt') as f:
        return [line.strip() for line in f if line.strip()]


def strs(xs):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in xs]


def rho(a, b):
    if len(a) < 3 or np.unique(a).size < 2 or np.unique(b).size < 2:
        return None
    value = float(spearmanr(a, b).statistic)
    return value if np.isfinite(value) else None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--loom', type=Path, required=True)
    ap.add_argument('--sample-dir', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--sample', required=True)
    args = ap.parse_args()
    folder = args.sample_dir / 'filtered_matrix/sensitivity_5'
    gene_records = [v.split('\t') for v in lines(find(folder, ['features.tsv.gz', 'features.tsv', 'genes.tsv.gz', 'genes.tsv']))]
    original_genes = [r[0] for r in gene_records]
    original_names = [r[1] if len(r) > 1 else r[0] for r in gene_records]
    bcs = lines(find(folder, ['barcodes.tsv.gz', 'barcodes.tsv']))
    X = mmread(str(find(folder, ['matrix.mtx.gz', 'matrix.mtx']))).tocsr()
    if X.shape == (len(bcs), len(original_genes)):
        X = X.T.tocsr()
    if X.shape != (len(original_genes), len(bcs)):
        raise ValueError(f'Unexpected original matrix shape {X.shape}')
    original_gene_totals = np.asarray(X.sum(axis=1), dtype=np.float64).ravel()
    original_cell_totals = np.asarray(X.sum(axis=0), dtype=np.float64).ravel()
    original_total = float(original_cell_totals.sum())
    with h5py.File(args.loom, 'r') as f:
        genes = strs(f['row_attrs/Accession'][:])
        cells_raw = strs(f['col_attrs/CellID'][:])
        extracted = [re.search(r'[ACGT]{16}', c) for c in cells_raw]
        if any(v is None for v in extracted):
            raise ValueError('Cannot parse 16bp barcode from all Loom CellID values')
        cells = [m.group() for m in extracted]
        if len(set(cells)) != len(cells):
            raise ValueError('Duplicate cell barcodes in Loom')
        if set(cells) != set(bcs):
            raise ValueError(f'Cell barcode mismatch: loom={len(set(cells))}, original={len(set(bcs))}')
        if len(genes) != len(original_genes) or set(genes) != set(original_genes):
            raise ValueError('Loom genes do not match original PIPseeker feature IDs')
        n_genes, n_cells = len(genes), len(cells)
        gene_totals = np.zeros(n_genes, dtype=np.float64)
        cell_totals = np.zeros(n_cells, dtype=np.float64)
        layer_totals, layer_cell_totals = {}, {}
        for layer in ('spliced', 'unspliced', 'ambiguous'):
            ds = f['layers'][layer]
            if ds.shape != (n_genes, n_cells):
                raise ValueError(f'Unexpected layer {layer}: {ds.shape}')
            layer_cells = np.zeros(n_cells, dtype=np.float64)
            for start in range(0, n_genes, 512):
                block = np.asarray(ds[start:start+512, :], dtype=np.float64)
                gene_totals[start:start+block.shape[0]] += block.sum(axis=1)
                layer_cells += block.sum(axis=0)
            layer_totals[layer] = float(layer_cells.sum())
            layer_cell_totals[layer] = layer_cells
            cell_totals += layer_cells
    gm = {g: i for i, g in enumerate(original_genes)}
    bm = {bc: i for i, bc in enumerate(bcs)}
    original_gene_aligned = original_gene_totals[[gm[g] for g in genes]]
    original_cell_aligned = original_cell_totals[[bm[c] for c in cells]]
    gene_mask = original_gene_aligned >= 10
    cell_mask = original_cell_aligned > 0
    ratios = np.divide(cell_totals[cell_mask], original_cell_aligned[cell_mask])
    mito_ids = {original_genes[i] for i, n in enumerate(original_names) if n.lower().startswith('mt-')}
    mito_idx = [i for i, g in enumerate(genes) if g in mito_ids]
    velocity_total = float(cell_totals.sum())
    report = {
        'sample': args.sample,
        'status': 'complete_exploratory',
        'genes': n_genes,
        'cells': n_cells,
        'matching_genes': len(genes),
        'matching_cells': len(cells),
        'original_total_counts': original_total,
        'velocyto_total_counts': velocity_total,
        'global_count_ratio': velocity_total/original_total if original_total else None,
        'layer_totals': layer_totals,
        'layer_fractions': {k: v/velocity_total if velocity_total else None for k,v in layer_totals.items()},
        'cells_spliced_and_unspliced': int(np.count_nonzero((layer_cell_totals['spliced']>0) & (layer_cell_totals['unspliced']>0))),
        'gene_spearman_orig_ge10': rho(original_gene_aligned[gene_mask], gene_totals[gene_mask]),
        'gene_spearman_n': int(gene_mask.sum()),
        'cell_spearman': rho(original_cell_aligned[cell_mask], cell_totals[cell_mask]),
        'cell_spearman_n': int(cell_mask.sum()),
        'per_cell_count_ratio_median': float(np.median(ratios)),
        'per_cell_count_ratio_p5': float(np.percentile(ratios, 5)),
        'per_cell_count_ratio_p95': float(np.percentile(ratios, 95)),
        'mt_gene_count': len(mito_idx),
        'mt_total_counts': float(gene_totals[mito_idx].sum()),
        'disclaimer': 'PIPseq v5 IMI/bin amplification correction is NOT reproduced by velocyto; exploratory splicing layers only; repeat mask absent unless explicitly supplied.',
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('w') as f:
        json.dump(report, f, indent=2)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
