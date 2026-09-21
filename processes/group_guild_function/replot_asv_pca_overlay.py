#!/usr/bin/env python3
"""Create the labeled variant from existing outputs without rerunning analyses."""
import argparse
import json
import hashlib
from pathlib import Path
import pandas as pd
from group_guild_function import (
    read_table, merge_basin_pca_background, select_pca_label_asvs,
    plot_asv_module_pca_overlay, attach_label_taxonomy, attach_label_module_peaks,
)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--group-guild-dir', type=Path, required=True)
    ap.add_argument('--abundance-table', type=Path, required=True,
                    help='Full selected abundance matrix, with samples in rows and ASVs in columns.')
    ap.add_argument('--pca-scores', type=Path, required=True)
    ap.add_argument('--pca-explained', type=Path, required=True)
    ap.add_argument('--hybrid-assignments', type=Path, required=True)
    ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--taxonomy', type=Path, help='SILVA taxonomy table; defaults to the sibling taxonomy output.')
    ap.add_argument('--label-table', type=Path, help='Optional path for the selection audit TSV.')
    ap.add_argument('--min-mean-abundance', type=float, default=0.01)
    ap.add_argument('--formats', default='pdf,png,svg')
    args = ap.parse_args()
    tables = args.group_guild_dir / 'tables'
    counts = pd.read_csv(args.abundance_table, sep='\t', index_col=0).T
    positions = read_table(tables / 'ecological_module_asv_pca_positions.tsv')
    labels = select_pca_label_asvs(counts, positions, args.min_mean_abundance)
    taxonomy_path = args.taxonomy or args.group_guild_dir.parent / 'taxonomy/tables/ASV_SILVA_tax.full-length.vsearch.tsv'
    labels = attach_label_taxonomy(labels, read_table(taxonomy_path))
    module_stats_path = tables / 'compartment_matched_ecological_module_association.tsv'
    labels = attach_label_module_peaks(labels, read_table(module_stats_path))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    label_table = args.label_table or args.output_dir / 'ecological_module_asv_pca_overlay_labels.tsv'
    label_table.parent.mkdir(parents=True, exist_ok=True)
    labels.to_csv(label_table, sep='\t', index=False)
    background = merge_basin_pca_background(read_table(args.pca_scores), read_table(args.hybrid_assignments))
    for suffix, annotations in [('asv_level', None), ('asv_level_abundant_asvs_labeled', labels)]:
        plot_asv_module_pca_overlay(
            background,
            read_table(tables / 'ecological_module_pca_hybrid_centroids.tsv'),
            positions,
            read_table(args.pca_explained),
            args.output_dir / ('ecological_module_asv_pca_overlay_' + suffix),
            args.formats.split(','), label_positions=annotations,
        )
    source_paths = [module_stats_path, taxonomy_path, args.abundance_table, args.pca_scores, args.pca_explained,
                    args.hybrid_assignments, tables/'ecological_module_asv_pca_positions.tsv',
                    tables/'ecological_module_pca_hybrid_centroids.tsv', Path(__file__),
                    Path(__file__).with_name('group_guild_function.py')]
    provenance = {'arguments': vars(args), 'selected_asvs': labels.ASV_ID.tolist(),
                  'sources': [{'path': str(p.resolve()), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
                              for p in source_paths]}
    (args.output_dir/'ecological_module_asv_pca_overlay_render.json').write_text(
        json.dumps(provenance, indent=2, default=str))
    print(f'[done] Labeled {labels.plotted.sum()} of {len(labels)} selected ASVs; output: {args.output_dir}')


if __name__ == '__main__':
    main()
