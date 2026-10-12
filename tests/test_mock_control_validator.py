"""The mock validator rejects cohort leakage and numerical count subtraction."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CONTROL = load('control', 'processes/control_decontam/control_decontam.py')
VALIDATOR = load('validator', 'examples/mock_test/validate_results.py')


class ControlValidatorTest(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.directory = self.root / 'modules/contamination_filtering/tables/three_tier_results'
        raw_dir = self.root / 'intermediates/ASVs'
        raw_dir.mkdir(parents=True)
        self.counts = pd.DataFrame({'airway': [5000, 10, 10, 10], 'oral': [5000, 0, 0, 0],
                                   'skin': [0, 0, 10, 10], 'blank': [0, 10, 0, 10],
                                   'low_airway': [4999, 0, 0, 0]},
                                  index=['clean', 'tech', 'bio', 'both'])
        CONTROL.save_counts(self.counts, raw_dir / 'ASV_counts.tsv')
        self.meta = pd.DataFrame({'sample_id': self.counts.columns,
                                  'Type_Group': ['Airways', 'Oral', 'Skin', 'Control', 'Airways'],
                                  'Participant_ID': ['P1', 'P1', 'P1', 'blank', 'P2']})
        self.meta.to_csv(self.root / 'metadata.tsv', sep='\t', index=False)
        config = dict(sample_col='sample_id', class_col='Type_Group', min_biological_reads=5000,
                      labels=dict(biological=['Airways', 'Oral'], technical=['Control'],
                                  bio_control=['Skin'], positive=['Positive']),
                      technical=dict(enabled=True, threshold=0.1), bio_control=dict(enabled=True, threshold=0.1))
        CONTROL.prepare(raw_dir/'ASV_counts.tsv', self.root/'metadata.tsv', config, self.directory)
        for arm, flags in [('TECH', [False, True, False, True]), ('BIO', [False, False, True, True])]:
            pd.DataFrame({'ASV_ID': self.counts.index, 'score': [0.01 if x else 0.9 for x in flags],
                          'contaminant': flags, 'status': 'tested'}).to_csv(self.directory/f'{arm}_scores.tsv', sep='\t', index=False)
        pd.DataFrame({'Feature ID': self.counts.index, 'Taxon': 'd__Bacteria'}).to_csv(self.root/'tax.tsv', sep='\t', index=False)
        CONTROL.finalize(self.root/'tax.tsv', self.directory)
        final_dir = self.root/'modules/non_target_filtering/tables'
        final_dir.mkdir(parents=True)
        CONTROL.save_counts(self.counts.loc[['clean'], ['airway', 'oral']], final_dir/'ASV_target.tsv')
        micro = self.counts.loc[['clean'], ['airway','oral']]
        CONTROL.save_counts(micro, final_dir/'ASV_target.micro.tsv')
        filtering = load('feature_filter', 'processes/filter_counts/filter_nontarget.py')
        qc = filtering.biological_feature_qc(micro, 0.1, 0.05)
        qc['retained_final'] = True
        qc.to_csv(final_dir/'ASV_target.feature_qc.tsv',sep='\t')
        config_dir = self.root/'summary/tables'
        config_dir.mkdir(parents=True)
        (config_dir/'run_config.yml').write_text('standard:\n  filter_counts:\n    min_relative_abundance_pct: 0.1\n    min_prevalence_fraction: 0.05\n')
        import gzip
        with gzip.open(raw_dir/'ASVs_target.fasta.gz', 'wt') as handle:
            handle.write('>clean\nACGT\n')

    def validate(self):
        audit = VALIDATOR.Audit()
        VALIDATOR.validate_control_filtering(self.meta, self.root, audit, cami=True)
        return audit.failures

    def test_valid_fixture_and_low_depth_skin_retained(self):
        self.assertEqual(self.validate(), [])

    def test_detects_technical_controls_in_bio_arm(self):
        path = self.directory/'BIO_counts.tsv'
        table = pd.read_csv(path, sep='\t', index_col=0)
        table['blank'] = self.counts.blank
        CONTROL.save_counts(table, path)
        self.assertTrue(any('BIO_cohort' in failure for failure in self.validate()))

    def test_detects_control_in_prevalence_denominator(self):
        path = self.root/'modules/non_target_filtering/tables/ASV_target.feature_qc.tsv'
        table = pd.read_csv(path,sep='\t',index_col=0)
        table['n_biological_samples'] += 1
        table.to_csv(path,sep='\t')
        self.assertTrue(any('biological_prevalence_abundance' in x for x in self.validate()))

    def test_detects_subtraction(self):
        path = self.directory/'ASV_cleaned.tsv'
        table = pd.read_csv(path, sep='\t', index_col=0)
        table.loc['clean', 'airway'] -= 10
        CONTROL.save_counts(table, path)
        self.assertTrue(any('contamination_union' in failure for failure in self.validate()))


if __name__ == '__main__':
    unittest.main()
