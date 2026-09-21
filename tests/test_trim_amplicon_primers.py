import importlib.util
import os
from pathlib import Path
import shutil
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/trim_amplicon_primers.py'
spec = importlib.util.spec_from_file_location('primer_trim', SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def write_pairs(root, pairs, sample='sample'):
    for mate in (1, 2):
        with (root / f'{sample}_R{mate}.fastq').open('w') as f:
            for i, pair in enumerate(pairs):
                seq = pair[mate-1]
                f.write(f'@read{i}/{mate}\n{seq}\n+\n'+ 'I'*len(seq)+'\n')


class PrimerTests(unittest.TestCase):
    def test_discovery(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            for name in ['a_R1_001.fastq.gz','a_R2_001.fastq.gz','b_1.fq','b_2.fq']:
                (p/name).touch()
            self.assertEqual(set(m.discover(p)), {'a_001','b'})
            (p/'orphan_R1.fastq').touch()
            with self.assertRaisesRegex(ValueError,'Missing paired'):m.discover(p)

    def test_mismatched_ids(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);write_pairs(p,[('ACGT','ACGT')])
            r2=p/'sample_R2.fastq';r2.write_text(r2.read_text().replace('read0','other'))
            with self.assertRaisesRegex(ValueError,'out of sync'):list(m.paired_records(p/'sample_R1.fastq',r2))

    @unittest.skipUnless(shutil.which(os.environ.get('CUTADAPT','cutadapt')), 'Cutadapt required')
    def test_unsupported(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);raw=p/'raw';raw.mkdir();out=p/'out'
            write_pairs(raw,[('A'*80,'C'*80)]*10)
            code=m.main(['--fastq-dir',str(raw),'--output-dir',str(out),'--cutadapt',os.environ.get('CUTADAPT','cutadapt')])
            self.assertEqual(code,2)
            self.assertEqual(list((out/'trimmed').iterdir()),[])
            self.assertIn('unsupported',(out/'summary.tsv').read_text())

    @unittest.skipUnless(shutil.which(os.environ.get('CUTADAPT','cutadapt')), 'Cutadapt required')
    def test_cutadapt_end_to_end(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);raw=p/'raw';raw.mkdir();out=p/'out'
            forward='GTGCCAGCAGCCGCGGTAA'
            reverse='GGACTACCGGGGTTTCTAAT'
            reverse2='CCGTCAATT CMTTTRAGTTT'.replace(' ','').replace('M','A').replace('R','A')
            insert='TACGATCGATGCTAGCTAGCATCGATCGAT'
            rc=lambda s:s.translate(m.COMPLEMENT)[::-1]
            pairs=[('GAATT'+forward+insert,'T'+reverse+insert),
                   (forward+insert+rc(reverse),reverse+insert+rc(forward)),
                   (forward+insert+rc(reverse2),reverse2+insert+rc(forward)),
                   (reverse+insert,forward+insert),
                   (forward,reverse),
                   (forward+insert,'A'*50),
                   (forward[:5]+'T'+forward[6:]+insert,reverse+insert)]
            write_pairs(raw,pairs)
            code=m.main(['--fastq-dir',str(raw),'--output-dir',str(out),'--cutadapt',os.environ.get('CUTADAPT','cutadapt'),'--threads','1','--min-family-reads','1'])
            self.assertEqual(code,0)
            for mate in (1,2):
                records=list(m.fastq(out/'trimmed'/f'sample_R{mate}.fastq.gz'))
                self.assertEqual(len(records),5)
                self.assertEqual([r[1] for r in records],[insert]*5)
                self.assertEqual(len(list(m.fastq(out/'unmatched'/f'sample_R{mate}.fastq.gz'))),1)
                # Empty processed reads are valid in this bin; count records directly.
                import gzip
                with gzip.open(out/'too_short'/f'sample_R{mate}.fastq.gz','rt') as f:self.assertEqual(len(f.readlines()),4)
            with self.assertRaisesRegex(ValueError,'not empty'):m.main(['--fastq-dir',str(raw),'--output-dir',str(out)])

if __name__=='__main__':unittest.main()
