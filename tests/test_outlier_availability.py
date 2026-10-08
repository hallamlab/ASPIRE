import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('outlier',ROOT/'processes/outlier_checker/outlier_checker.py')
AVAILABLE = all(importlib.util.find_spec(name) is not None for name in ['hdbscan','skbio','matplotlib','seaborn','sklearn'])
m = None
if AVAILABLE:
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class NoClusters:
    def __init__(self,**kwargs): pass
    def fit(self,x): self.labels_=np.full(len(x),-1); return self

@unittest.skipUnless(AVAILABLE, "Outlier conda environment is required")
class AvailabilityTests(unittest.TestCase):
    def test_no_clusters_do_not_supply_votes(self):
        with patch.object(m.hdbscan,'HDBSCAN',NoClusters), patch.object(m,'approximate_predict') as predict:
            result=m.run_for_group('BAL',pd.DataFrame([[0,1],[1,0],[1,1]],index=['a','b','c']),pd.DataFrame(index=['a','b','c']),False,(True,True,True),{'random_state':1},{},{},3)
        predict.assert_not_called()
        self.assertTrue(result.HDBSCAN.isna().all())
        self.assertTrue((result.available_detectors==2).all())
        self.assertTrue(result.is_outlier.isna().all())
        self.assertTrue((result.outlier_votes<=2).all())

    def test_singleton_is_unavailable(self):
        result=m.run_for_group('tiny',pd.DataFrame([[1]],index=['a']),pd.DataFrame(index=['a']),False,(True,True,True),{},{},{},3)
        self.assertTrue(result.is_outlier.isna().all())
        self.assertEqual(result.available_detectors.iloc[0],0)

if __name__=='__main__':unittest.main()
