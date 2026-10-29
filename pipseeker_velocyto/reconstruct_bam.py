#!/usr/bin/env python3
"""Stream a PIPseeker v5 STAR BAM as a tagged BAM on stdout.

Research-only reconstruction, validated initially on sample1_GRCm39_Cre.
CB = verified 2-bit forward encoding of 16bp barcode in QNAME field 2.
UB = reversible 15-base encoding of candidate packed MI in field 3.
Only STAR MAPQ=255 primary alignments retained; NH=1 assigned provisionally.

This does NOT reproduce PIPseeker v5 molecule-level amplification correction.
"""
import argparse
import csv
import gzip
import sys
from collections import Counter
from pathlib import Path

import pysam

BASES = 'ACGT'
ENC = {b: i for i, b in enumerate(BASES)}


def pack_barcode(seq: str) -> int:
    v = 0
    for b in seq:
        v = (v << 2) | ENC[b]
    return v


def decode_packed_mi(v: int) -> str:
    if v < 0 or v >= (1 << 30):
        raise ValueError('candidate MI outside 30-bit range')
    return ''.join(BASES[(v >> (2 * i)) & 3] for i in range(14, -1, -1))


def load_selected(path):
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt') as f:
        seqs = [x.strip() for x in f if x.strip()]
    if any(len(x) != 16 or any(ch not in ENC for ch in x) for x in seqs):
        raise RuntimeError('Filtered barcodes are not 16bp A/C/G/T strings')
    if len(seqs) != len(set(seqs)):
        raise RuntimeError('Duplicate filtered barcodes')
    return {pack_barcode(x): x for x in seqs}


def config_check(folder):
    with (folder / 'run_config.csv').open(newline='') as f:
        rows = {row[0]: row[1] for row in csv.reader(f) if len(row) >= 2}
    if rows.get('chemistry', '').strip().upper() != 'V':
        raise RuntimeError(f"Unsupported chemistry: {rows.get('chemistry')!r}; validated only for V")
    if 'GRCm39_2024A_Cre' not in rows.get('star_index_path', ''):
        raise RuntimeError('Unexpected STAR reference; verify against exact GTF before using this pipeline')
    return rows


def validate_names(bam, wl_path, n=5000):
    # Preflight on a small read sample: barcode-encoding hypothesis must generalize.
    with open(wl_path) as f:
        valid_codes = {pack_barcode(x.strip()) for x in f if x.strip()}
    total = matched = 0
    with pysam.AlignmentFile(str(bam), 'rb') as fh:
        for aln in fh.fetch(until_eof=True):
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            parts = aln.query_name.split(':')
            total += 1
            if len(parts) == 3:
                try:
                    matched += int(int(parts[1]) in valid_codes)
                except ValueError:
                    pass
            if total >= n:
                break
    rate = matched / total if total else 0
    print(f'Preflight mapped primary reads: {total:,}; valid whitelist barcode codes: {matched:,} ({rate:.2%})', file=sys.stderr)
    if total < 100 or rate < 0.99:
        raise RuntimeError('Barcode encoding failed preflight validation; no BAM emitted')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--sample-dir', type=Path, required=True)
    ap.add_argument('--summary', type=Path, required=True)
    ap.add_argument('--preflight-reads', type=int, default=5000)
    args = ap.parse_args()
    sample_dir = args.sample_dir
    config_check(sample_dir)
    bam = sample_dir / 'star_out.bam'
    selected_bc = sample_dir / 'filtered_matrix/sensitivity_5/barcodes.tsv.gz'
    whitelist = sample_dir / 'barcodes/barcode_whitelist.txt'
    for p in (bam, selected_bc, whitelist):
        if not p.is_file():
            raise FileNotFoundError(p)
    selected = load_selected(selected_bc)
    print(f'Sensitivity 5 selected cells: {len(selected):,}', file=sys.stderr)
    validate_names(bam, whitelist, args.preflight_reads)

    stats = Counter()
    with pysam.AlignmentFile(str(bam), 'rb') as src:
        with pysam.AlignmentFile('-', 'wb', header=src.header) as dst:
            for read in src.fetch(until_eof=True):
                stats['examined_alignments'] += 1
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                stats['primary_mapped'] += 1
                if read.mapping_quality != 255:
                    continue
                stats['unique_primary_mapq255'] += 1
                parts = read.query_name.split(':')
                if len(parts) != 3:
                    stats['invalid_qname'] += 1
                    continue
                try:
                    bc_code, mi_code = int(parts[1]), int(parts[2])
                except ValueError:
                    stats['invalid_qname'] += 1
                    continue
                bc = selected.get(bc_code)
                if bc is None:
                    continue
                stats['reads_in_selected_cells'] += 1
                try:
                    ub = decode_packed_mi(mi_code)
                except ValueError:
                    stats['invalid_mi'] += 1
                    continue
                read.set_tag('CB', bc, value_type='Z')
                read.set_tag('UB', ub, value_type='Z')
                read.set_tag('NH', 1, value_type='i')
                dst.write(read)
                stats['written_alignments'] += 1
                if stats['examined_alignments'] % 5_000_000 == 0:
                    print(f"Examined {stats['examined_alignments']:,} | written {stats['written_alignments']:,}", file=sys.stderr, flush=True)
    if stats['written_alignments'] == 0:
        raise RuntimeError('No alignments selected')
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    import json
    with args.summary.open('w') as f:
        json.dump({'sample': sample_dir.name, 'sensitivity': 5,
                   'called_cells': len(selected), **stats}, f, indent=2)
    print('Final reconstruction statistics:', dict(stats), file=sys.stderr)


if __name__ == '__main__':
    main()
