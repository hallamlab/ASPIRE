#!/usr/bin/env python3
"""Select an analysis cohort without changing upstream biological sample roles."""
import argparse
import json
from pathlib import Path
import pandas as pd


def select_cohort(metadata, counts, sample_col, group_col, exclude_groups):
    if sample_col not in metadata or group_col not in metadata:
        raise ValueError(f'Metadata must contain {sample_col!r} and {group_col!r}')
    metadata = metadata.copy()
    metadata[sample_col] = metadata[sample_col].astype(str)
    if metadata[sample_col].duplicated().any():
        raise ValueError('Duplicate metadata sample identifiers')
    counts = counts.copy()
    counts.columns = counts.columns.astype(str)
    if set(counts.columns) != set(metadata[sample_col]):
        raise ValueError('Cohort metadata and count-table samples must match exactly')
    keep = ~metadata[group_col].isin(exclude_groups)
    audit = metadata[[sample_col, group_col]].copy()
    audit['included'] = keep
    audit['reason'] = keep.map({True: 'analysis cohort', False: 'excluded study group'})
    selected_meta = metadata.loc[keep].copy()
    if selected_meta.empty:
        raise ValueError('Analysis cohort selection removed every sample')
    selected = counts.loc[:, selected_meta[sample_col]]
    selected = selected.loc[selected.sum(axis=1) > 0]
    return selected_meta, selected, audit


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('metadata', 'counts', 'settings', 'outdir'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--long', type=Path)
    p.add_argument('--clr', type=Path)
    args = p.parse_args()
    cfg = json.loads(args.settings.read_text())
    sample_col = cfg['sample_col']
    meta = pd.read_csv(args.metadata, sep='\t', dtype={sample_col: str})
    counts = pd.read_csv(args.counts, sep='\t', index_col=0)
    selected_meta, selected, audit = select_cohort(meta, counts, sample_col, cfg['group_col'], cfg['exclude_groups'])
    args.outdir.mkdir(parents=True, exist_ok=True)
    selected_meta.to_csv(args.outdir/'metadata.tsv', sep='\t', index=False)
    selected.to_csv(args.outdir/'counts.tsv', sep='\t')
    audit.to_csv(args.outdir/'sample_selection.tsv', sep='\t', index=False)
    if args.long:
        long = pd.read_csv(args.long, sep='\t', dtype={sample_col: str, 'ASV_ID': str})
        positive_samples = set(selected.columns[selected.sum(axis=0) > 0])
        if not positive_samples.issubset(set(long[sample_col])):
            raise ValueError('Selected samples are missing from the long ASV table')
        long = long.loc[long[sample_col].isin(selected_meta[sample_col]) & long.ASV_ID.isin(selected.index.astype(str))]
        long.to_csv(args.outdir/'asv_metadata.tsv', sep='\t', index=False)
    if args.clr:
        clr = pd.read_csv(args.clr, sep='\t', index_col=0)
        clr.index = clr.index.astype(str)
        # Preserve the full feature basis of the supplied CLR transformation.
        clr.loc[selected_meta[sample_col]].to_csv(args.outdir/'clr.tsv', sep='\t')
    (args.outdir/'summary.json').write_text(json.dumps({
        'input_samples': len(meta), 'analysis_samples': len(selected_meta),
        'input_asvs': len(counts), 'analysis_asvs': len(selected),
        'excluded_groups': cfg['exclude_groups'],
        'diversity_includes_excluded_groups': True,
    }, indent=2)+'\n')


if __name__ == '__main__':
    main()
