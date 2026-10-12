import importlib.util
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('cohort', ROOT/'processes/analysis_cohort/select_cohort.py')
cohort = importlib.util.module_from_spec(spec); spec.loader.exec_module(cohort)


def test_skin_is_removed_only_from_selected_copy():
    meta = pd.DataFrame({'Sample':['a','b','c'], 'Type_Group':['BAL','Skin Brush','Oral Rinse']})
    counts = pd.DataFrame({'a':[3,0], 'b':[0,7], 'c':[2,0]}, index=['shared','skin_only'])
    selected_meta, selected, audit = cohort.select_cohort(meta, counts, 'Sample','Type_Group',['Skin Brush'])
    assert list(selected.columns) == ['a','c']
    assert list(selected.index) == ['shared']
    assert selected.values.sum() == 5
    assert len(meta) == 3 and counts.values.sum() == 12  # diversity input remains intact
    assert audit.included.tolist() == [True,False,True]
    assert selected_meta.Sample.tolist() == ['a','c']


def test_mismatched_ids_fail():
    with pytest.raises(ValueError, match='match exactly'):
        cohort.select_cohort(pd.DataFrame({'Sample':['x'],'Type_Group':['BAL']}), pd.DataFrame({'y':[1]}), 'Sample','Type_Group',[])
