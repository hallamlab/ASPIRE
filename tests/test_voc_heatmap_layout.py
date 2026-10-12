"""Verify exported ASV-VOC plot axes, clustering and annotation alignment."""
import importlib.util
from pathlib import Path

import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / 'processes/voc_correlation/plot_voc_corr.py'
spec = importlib.util.spec_from_file_location('voc_heatmap', SCRIPT)
voc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(voc)


@pytest.mark.parametrize('single_voc', [False, True])
def test_voc_rows_asv_columns_clustered_with_matching_annotations(tmp_path, monkeypatch, single_voc):
    matrix = pd.DataFrame({'acetone': [.2, .8, -.4], 'ethanol': [.9, -.1, .2]}, index=['ASV1', 'ASV2', 'ASV3'])
    if single_voc:
        matrix = matrix[['acetone']]
    colors = pd.Series(['#00ff00', '#ff0000', '#0000ff'], index=['ASV3','ASV1','ASV2'], name='ISA group')
    calls = []
    real = voc.sns.clustermap
    def capture(data, **kwargs):
        grid = real(data, **kwargs)
        calls.append((data.copy(), kwargs, grid))
        return grid
    monkeypatch.setattr(voc.sns, 'clustermap', capture)
    monkeypatch.setattr(voc, 'FIGURE_FORMATS', ('.svg',))
    stem = tmp_path / 'heatmap'
    voc.save_asv_voc_clustermap(matrix, stem, row_colors=colors,
        row_color_legend=[('Airways','#00ff00')], correlation_direction='both')
    data, settings, grid = calls[0]
    pd.testing.assert_frame_equal(data, matrix.T.rename_axis(index='VOC',columns='ASV'))
    assert settings['row_cluster'] is (not single_voc)
    assert settings['col_cluster'] is True
    assert settings['row_colors'] is None
    pd.testing.assert_series_equal(settings['col_colors'].iloc[:,0], colors)
    assert grid.ax_heatmap.get_xlabel() == 'ASV'
    assert grid.ax_heatmap.get_ylabel() == 'VOC'
    assert settings['figsize'][0] > settings['figsize'][1]
    assert grid.dendrogram_col is not None
    if not single_voc:
        assert grid.dendrogram_row is not None
    text = stem.with_suffix('.svg').read_text()
    assert '>ASV<' in text and '>VOC<' in text
    assert 'ISA group' in text


def test_output_rename_removes_only_known_retired_files(tmp_path):
    old = tmp_path / 'sample_voc_brush_clustermap.svg'
    unrelated = tmp_path / 'custom_brush_measurements.tsv'
    current = tmp_path / 'sample_voc_clustermap.svg'
    for path in (old, unrelated, current):
        path.write_text('keep or retire')
    voc.remove_legacy_outputs(tmp_path)
    assert not old.exists()
    assert unrelated.is_file() and current.is_file()
