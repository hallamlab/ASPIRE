import ast
from pathlib import Path
import tempfile
import unittest

import pandas as pd


SCRIPT = Path(__file__).parents[1] / "processes" / "voc_correlation" / "plot_voc_corr.py"


def load_normalizer():
    tree = ast.parse(SCRIPT.read_text())
    function = next(
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "normalize_sample_id"
    )
    module = ast.Module(body=[function], type_ignores=[])
    namespace = {}
    exec(compile(module, str(SCRIPT), "exec"), namespace)
    return namespace["normalize_sample_id"]


class VocSampleIdTests(unittest.TestCase):
    def test_exact_preserves_complete_sample_id(self):
        normalize = load_normalizer()
        self.assertEqual(normalize(" CASE_001_BRUSH_CONTRA ", "exact"), "CASE_001_BRUSH_CONTRA")

    def test_legacy_and_prefix_modes_remain_supported(self):
        normalize = load_normalizer()
        self.assertEqual(normalize("P001_BRUSH_LEFT", "legacy_patient_pair"), "P001_BRUSH")
        self.assertEqual(normalize("CASE_001_BRUSH", "prefix1"), "CASE")

    def test_voc_annotations_are_optional(self):
        tree = ast.parse(SCRIPT.read_text())
        functions = [
            node for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name in {"normalize_sample_id", "load_voc_metadata"}
        ]
        namespace = {"pd": pd}
        exec(compile(ast.Module(body=functions, type_ignores=[]), str(SCRIPT), "exec"), namespace)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "voc.tsv"
            pd.DataFrame({"sample_id": ["S1"], "compound": [1.0]}).to_csv(path, sep="\t", index=False)
            metadata = namespace["load_voc_metadata"](str(path), "sample_id", "exact")
        self.assertEqual(metadata.loc["S1", "subclass2"], "missing")
        self.assertEqual(metadata.loc["S1", "Type"], "missing")


if __name__ == "__main__":
    unittest.main()


def test_airway_union_is_nondistinct_even_when_skin_is_biological():
    tree = ast.parse(SCRIPT.read_text())
    names = {'parse_csv_list', 'canonicalize_sample_type', 'isa_group_is_focus_specific'}
    funcs = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    ns = {'pd': pd}
    exec(compile(ast.Module(body=funcs, type_ignores=[]), str(SCRIPT), 'exec'), ns)
    airway = {'Oral Rinse','BAL','Bronchial Brush'}
    f = ns['isa_group_is_focus_specific']
    assert f('BAL+Bronchial Brush',airway,airway,True)
    assert not f('Oral Rinse+BAL+Bronchial Brush',airway,airway,True)
    assert not f('Skin Brush+Oral Rinse+BAL+Bronchial Brush',airway,airway,True)
