"""Preserve both filtering layers and count/sequence agreement."""
from pathlib import Path
import subprocess
import importlib.util
import sys
import tempfile
import unittest

import pandas as pd
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]


class CombinedFilterTest(unittest.TestCase):
    def test_unassigned_blank_and_missing_taxonomy_are_removed(self):
        spec = importlib.util.spec_from_file_location('taxonomy_filter', ROOT/'processes/filter_counts/filter_nontarget.py')
        module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        labels = ['Unassigned', ' unASSIGNED ', '', None, 'd__; p__; g__',
                  'd__Unassigned; p__unknown', 'd__Bacteria; g__unclassified',
                  'd__Archaea', 'd__Eukaryota; p__Vertebrata; c__Mammalia']
        ids = [f'ASV{i}' for i in range(len(labels))]
        counts = pd.DataFrame({'S1': range(1, len(ids)+2)}, index=ids+['absent'])
        taxonomy = pd.DataFrame({'Taxon': labels, 'Consensus': 1.0}, index=ids)
        actual = module.filter_by_taxonomy(counts, taxonomy, 'Taxon', 'Consensus', .05,
                                          module.parse_exclude_taxa(['Class:Mammalia']))
        pd.testing.assert_frame_equal(actual, counts.loc[['ASV6', 'ASV7']])

    def test_both_layers_and_final_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            ids = ['keep', 'nontarget', 'mito', 'initial_rare', 'micro_rare', 'bad_tax', 'excluded', 'rare_group', 'zero']
            counts = pd.DataFrame({'A': [1000, 500, 500, 1, 100, 500, 500, 0, 0],
                                   'B': [1000, 500, 500, 1, 100, 500, 500, 0, 0],
                                   'R': [1000, 500, 500, 1, 100, 500, 500, 1000, 0]}, index=ids)
            counts.to_csv(p / 'raw.tsv', sep='\t')
            pd.DataFrame({'Sample': ['A', 'B', 'R', 'not_in_counts'], 'Group': ['main', 'main', 'rare', 'rare']}).to_csv(p / 'meta.tsv', sep='\t', index=False)
            pd.DataFrame({'BioFactorial': [1, 0, 1, 1, 1, 1, 1, 1, 1],
                          'MITOMASTER': [1, 1, 0, 1, 1, 1, 1, 1, 1]}, index=ids).to_csv(p / 'non.tsv', sep='\t')
            tax = pd.DataFrame({'Taxon': ['d__Bacteria'] * len(ids), 'Consensus': [1.] * len(ids)}, index=ids)
            tax.loc['bad_tax', 'Consensus'] = 0.1
            tax.loc['excluded', 'Taxon'] = 'd__Eukaryota'
            tax.to_csv(p / 'tax.tsv', sep='\t')
            (p / 'asvs.fa').write_text(''.join(f'>{i}\nACGT\n' for i in ids))
            def run(script, *args):
                subprocess.run([sys.executable, str(ROOT / script), *map(str, args)], check=True, capture_output=True, text=True)
            # Depth=0 because upstream biological inclusion has already run.
            run('processes/filter_table/filter_ASV_table.py', p/'raw.tsv', p/'initial.tsv', 0, 0.05, p/'asvs.fa', p/'initial.fa')
            run('processes/filter_counts/filter_nontarget.py', '--count-table', p/'initial.tsv',
                '--nontarget-table', p/'non.tsv', '--taxonomy-table', p/'tax.tsv',
                '--metadata', p/'meta.tsv', '--group-col', 'Group', '--min-group-size', 2,
                '--sample-id-col', 'Sample', '--mito-cols', 'MITOMASTER',
                '--abundance-threshold', 10, '--min-consensus', 0.8,
                '--exclude-taxon', 'Domain:Eukaryota', '--save-intermediates',
                '--mito-output-dir', p, '--output', p/'final.tsv')
            run('processes/filter_counts/finalize_filter.py', '--input-counts', p/'raw.tsv',
                '--table-counts', p/'initial.tsv', '--counts', p/'final.tsv',
                '--fasta', p/'initial.fa', '--output-fasta', p/'final.fa')
            final = pd.read_csv(p/'final.tsv', sep='\t', index_col=0)
            pd.testing.assert_frame_equal(final, counts.loc[['keep'], ['A', 'B']])
            self.assertEqual([r.id for r in SeqIO.parse(p/'final.fa', 'fasta')], ['keep'])
            self.assertNotIn('initial_rare', pd.read_csv(p/'initial.tsv', sep='\t', index_col=0).index)
            self.assertIn('micro_rare', pd.read_csv(p/'final.micro.tsv', sep='\t', index_col=0).index)
            self.assertTrue((p/'final.mito.tsv').exists())
            self.assertTrue((p/'filter_removed_asvs.tsv').exists())
            self.assertEqual(pd.read_csv(p/'filter_audit.tsv', sep='\t').iloc[-1]['samples'], 2)
            pd.testing.assert_frame_equal(pd.read_csv(p/'raw.tsv', sep='\t', index_col=0), counts)


if __name__ == '__main__':
    unittest.main()
