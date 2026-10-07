import importlib.util
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('nontarget', ROOT/'processes/filter_counts/filter_nontarget.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_five_percent_rounding_and_independent_abundance_gate():
    counts = pd.DataFrame(0, index=['background','two_samples','three_samples','too_rare','zero'],
                          columns=[f'bio{i}' for i in range(50)])
    counts.loc['background'] = 10000
    counts.loc['two_samples', ['bio0','bio1']] = [100, 1]
    counts.loc['three_samples', ['bio0','bio1','bio2']] = [100, 1, 1]
    counts.loc['too_rare', ['bio0','bio1','bio2']] = 1
    original = counts.copy()
    audit = module.biological_feature_qc(counts, 0.1, 0.05)
    assert not audit.loc['two_samples', 'passes_feature_filters']
    assert audit.loc['three_samples', 'passes_feature_filters']
    assert audit.loc['three_samples', 'n_nonzero_samples'] == 3
    assert not audit.loc['too_rare', 'passes_feature_filters']
    assert not audit.loc['zero', 'passes_feature_filters']
    assert audit.n_biological_samples.eq(50).all()
    pd.testing.assert_frame_equal(counts, original)


def test_zero_depth_biological_sample_stays_in_denominator():
    audit = module.biological_feature_qc(pd.DataFrame({'bio1':[5], 'bio2':[0]}), 0, 0.75)
    assert audit.iloc[0].prevalence_fraction == 0.5
    assert not audit.iloc[0].passes_feature_filters


@pytest.mark.parametrize('ra,prev', [(0.1,5), (-1,0.05), (101,0), (0,float('nan'))])
def test_invalid_units(ra, prev):
    with pytest.raises(ValueError):
        module.biological_feature_qc(pd.DataFrame({'bio1':[1]}),ra,prev)
