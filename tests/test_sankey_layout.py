"""HTML and SVG must preserve conservation, node order and separate loss lanes."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('sankey', ROOT/'processes/sankey/sankey_builder.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

class SankeyLayoutTest(unittest.TestCase):
    def test_layout_conservation_and_shared_svg(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'plot.html'
            module.build_sankey(['Input','QC','Cohort','Finished'],[1000,400,350,300],
                {'Study':900,'TECH':100,'Empty':0},{'BAL':100,'Oral':200,'TECH':0},
                {'Study':'blue','TECH':'red'},'Test',path,True,'snap',
                loss_groups={1:[('TECH controls',30),('BIO controls',20)]})
            flow=json.loads(path.with_suffix('.flow.json').read_text())
            self.assertIn(path.with_suffix('.svg').read_text().strip(),path.read_text())
            nodes=flow['nodes'];edges=flow['links']
            self.assertFalse(any(n['value']==0 for n in nodes))
            for i,n in enumerate(nodes):
                if n['lane']=='removed': self.assertGreaterEqual(n['y'],560)
                else: self.assertLess(n['y']+n['h'],520)
                incoming=sum(e['value'] for e in edges if e['target']==i)
                outgoing=sum(e['value'] for e in edges if e['source']==i)
                if incoming: self.assertEqual(incoming,n['value'])
                if outgoing: self.assertEqual(outgoing,n['value'])
            final=next(i for i,n in enumerate(nodes) if n['label']=='Finished')
            split=[e for e in edges if e['source']==final]
            self.assertEqual([e['sy'] for e in split], sorted(e['sy'] for e in split))
            self.assertIn("flow.arrangement==='snap'",path.read_text())

if __name__=='__main__':unittest.main()
