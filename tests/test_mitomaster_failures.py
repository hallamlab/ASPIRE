import importlib.util
import json
from pathlib import Path
from threading import Lock
import unittest
from unittest.mock import patch, Mock
import tempfile
try:
    import requests
except ImportError:
    requests = None

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('mitomaster',ROOT/'processes/mitomaster/mitomaster.py')
m = None
if requests is not None:
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

@unittest.skipIf(requests is None, "MITOMASTER requests environment required")
class MitomasterFailures(unittest.TestCase):
    def test_all_failed_continue_writes_audited_empty_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'a.fasta').write_text('>ASV1\nACGT\n');out=root/'result.tsv'
            with patch('sys.argv',['mitomaster','--data-dir',tmp,'--output-file',str(out),'--failure-policy','continue']), patch.object(m,'post_one',side_effect=requests.ConnectionError('offline')):
                m.main()
            self.assertEqual(out.read_text(),'Sequence_ID\thaplo\n')
            self.assertEqual(json.loads(out.with_suffix('.status.json').read_text())['status'],'unavailable')
            self.assertIn('offline',out.with_suffix('.failures.tsv').read_text())
            self.assertFalse(Path(str(out)+'.done').exists())

    def test_first_chunk_failure_does_not_strip_successful_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);out=root/'result.tsv';failures=[]
            with patch.object(m,'post_one',side_effect=[requests.Timeout('timeout'),'Sequence_ID\thaplo\nASV2\tH1\n']):
                counts=m.process_first_then_pool(None,[root/'a',root/'b'],'url','sequences','hsd',out,root/'done',1,Lock(),'first',failures=failures)
            self.assertEqual(counts,(1,1))
            self.assertTrue(out.read_text().startswith('Sequence_ID\thaplo\nASV2'))
            self.assertEqual((root/'done').read_text(),'b\n')

    def test_fail_policy_stops_after_auditing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'a.fasta').write_text('>ASV1\nACGT\n');out=root/'result.tsv'
            with patch('sys.argv',['mitomaster','--data-dir',tmp,'--output-file',str(out)]), patch.object(m,'post_one',side_effect=requests.Timeout('timeout')):
                with self.assertRaises(SystemExit) as raised: m.main()
            self.assertEqual(raised.exception.code,1)
            self.assertEqual(json.loads(out.with_suffix('.status.json').read_text())['failed_chunks'],1)

    def test_html_error_page_is_not_accepted_as_sequence_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta=Path(tmp)/'a.fasta';fasta.write_text('>ASV1\nACGT\n')
            session=Mock();session.post.return_value.text='<html>Service unavailable</html>'
            with self.assertRaises(ValueError):
                m.post_one(session,'https://example.invalid',fasta,'sequences','hsd')

    def test_retry_budget_and_timeout_are_applied(self):
        session=m.build_session(4,1,90,'test')
        retry=session.get_adapter('https://').max_retries
        self.assertEqual(retry.total,4)
        self.assertFalse(retry.respect_retry_after_header)
        captured={}
        def request(method,url,**kwargs): captured.update(kwargs)
        m._timeout_wrapper(request,90)('POST','https://example.invalid')
        self.assertEqual(captured['timeout'],90)
        session.close()

if __name__=='__main__':unittest.main()
