import importlib.util
from pathlib import Path
import pandas as pd
import pytest
pytest.importorskip("matplotlib")
pytest.importorskip("seaborn")
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('isa_plot',ROOT/'processes/indicspecies_plots/plot_indicspecies.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def test_multipatt_combination_index_is_not_a_bitmask():
    df=pd.DataFrame({'index':[3,5,11], 's.BAL':[0,1,1], 's.Bronchial Brush':[0,1,1], 's.Oral Rinse':[1,0,1], 's.Skin Brush':[0,0,0]})
    mapping=m.infer_index_map_from_sign_table(df,'index','p.value','stat')
    assert mapping[3]=='Oral Rinse'
    assert mapping[5]=='BAL+Bronchial Brush'
    assert mapping[11]=='BAL+Bronchial Brush+Oral Rinse'
    palette=m.augment_combo_palette(m.EXPLICIT_GROUP_TYPE_PALETTE,list(mapping.values()))
    assert palette[mapping[3]]=='#6A3D9A'
    assert palette[mapping[5]]=='#5CC8C8'
    assert palette[mapping[11]]=='#CBB6E9'


@pytest.mark.parametrize('empty_groups', [(1,), (2,), (1, 2)])
def test_cli_accepts_header_only_isa_results(tmp_path, empty_groups):
    """Full-union removal can legitimately leave either grouping with no rows."""
    import subprocess
    import sys
    for group in (1, 2):
        table = pd.DataFrame({
            'ASV': ['ASV_1'], 's.Cancer': [1], 's.Control': [0],
            'index': [1], 'stat': [.8], 'p.value': [.01], 'q.value': [.02],
        })
        if group in empty_groups:
            table = table.iloc[:0]
        table.to_csv(tmp_path / f'group{group}.tsv', sep='\t', index=False)
    output = tmp_path / 'plots'
    subprocess.run([
        sys.executable, str(ROOT/'processes/indicspecies_plots/plot_indicspecies.py'),
        '--group1-results', str(tmp_path/'group1.tsv'),
        '--group2-results', str(tmp_path/'group2.tsv'),
        '--group1-name', 'Site', '--group2-name', 'Case', '--outdir', str(output),
    ], check=True, capture_output=True, text=True)
    for group, name in [(1, 'site'), (2, 'case')]:
        result = pd.read_csv(output/f'{name}_ISA_enriched.tsv', sep='\t')
        assert result.empty == (group in empty_groups)
        if group not in empty_groups:
            assert (output/f'{name}_ISA_plot.svg').is_file()
