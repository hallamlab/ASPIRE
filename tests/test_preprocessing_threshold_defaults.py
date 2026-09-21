"""Check production threshold expressions using Nextflow's Groovy semantics."""
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which('nextflow')

@unittest.skipUnless(NEXTFLOW, 'Nextflow required')
class ThresholdDefaultsTests(unittest.TestCase):
    def test_explicit_zero_and_missing_values(self):
        source = (ROOT / 'workflow/modules/preprocessing.nf').read_text()
        cases = [('tableCfg', 'min_sample_sum', 5000, 100),
                 ('tableCfg', 'min_asv_sum', 0.01, 0.05),
                 ('unoiseCfg', 'min_size', 8, 3)]
        checks = []
        for variable, key, default, positive in cases:
            expr = re.search(r'\$\{(' + variable + r'\.' + key + r'[^}]*)\}', source).group(1)
            for value, expected in [('[:]', default), (f'[{key}:null]', default),
                                    (f'[{key}:{positive}]', positive)]:
                checks.append(f'{variable} = {value}; assert ({expr}) == {expected}')
            if variable == 'tableCfg':
                checks.append(f'{variable} = [{key}:0]; assert ({expr}) == 0')
        script = 'nextflow.enable.dsl=2\nworkflow {\ndef tableCfg\ndef unoiseCfg\n' + '\n'.join(checks) + '\nprintln "THRESHOLDS_OK"\n}\n'
        with tempfile.TemporaryDirectory() as d:
            (Path(d) / 'main.nf').write_text(script)
            result = subprocess.run([NEXTFLOW, 'run', 'main.nf', '-ansi-log', 'false'],
                                    cwd=d, text=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, timeout=90)
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn('THRESHOLDS_OK', result.stdout)

if __name__ == '__main__':
    unittest.main()
