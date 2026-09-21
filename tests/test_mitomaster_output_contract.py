from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "processes" / "mitomaster" / "mitomaster.py"
SPEC = importlib.util.spec_from_file_location("mitomaster_output_contract", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class MitomasterOutputContractTests(unittest.TestCase):
    def run_case(self, replies, files=2):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        for i in range(files):
            (root / f"chunk{i}.fasta").write_text(f">ASV{i}\nACGT\n")
        output = root / "mitomaster_output.tsv"
        argv = [str(SCRIPT), "--data-dir", str(root), "--output-file", str(output), "--overwrite"]
        with patch.object(sys, "argv", argv), patch.object(MODULE, "post_one", side_effect=replies) as post:
            rc = MODULE.main()
        import json
        status = json.loads(Path(str(output) + ".status.json").read_text())
        return rc, output, status, post.call_count

    def test_no_chunks_stops(self):
        rc, output, status, calls = self.run_case([], files=0)
        self.assertEqual(rc, 2)
        self.assertFalse(output.exists())
        self.assertEqual(status["status"], "failed")

    def test_connection_failure_stops_after_first_probe(self):
        rc, output, status, calls = self.run_case([MODULE.requests.ConnectionError("Connection refused")])
        self.assertEqual(rc, 2)
        self.assertFalse(output.exists())
        self.assertEqual(calls, 1)
        self.assertIn("Connection refused", status["chunks"][0]["error"])

    def test_partial_failure_is_not_a_complete_result(self):
        rc, output, status, calls = self.run_case([
            "Sequence_ID\thaplo\nASV0\tH1\n", MODULE.requests.Timeout("Read timed out")])
        self.assertEqual(rc, 2)
        self.assertFalse(output.exists())
        self.assertTrue(Path(str(output) + ".partial").exists())
        self.assertEqual(calls, 2)
        self.assertEqual(status["chunks"][1]["error_type"], "Timeout")

    def test_html_http_success_is_rejected(self):
        rc, output, status, calls = self.run_case(["<!DOCTYPE html><html>Service unavailable</html>"])
        self.assertEqual(rc, 2)
        self.assertFalse(output.exists())

    def test_valid_header_only_is_successful_no_hits(self):
        rc, output, status, calls = self.run_case(["Sequence_ID\thaplo\n"] * 2)
        self.assertEqual(rc, 0)
        self.assertEqual(status["status"], "success")
        self.assertTrue(output.exists())
        self.assertEqual(sum(r["result_rows"] for r in status["chunks"]), 0)

    def test_success_has_one_header_and_both_results(self):
        rc, output, status, calls = self.run_case([
            "Sequence_ID\thaplo\nASV0\tH1\n", "Sequence_ID\thaplo\nASV1\tH2\n"])
        self.assertEqual(rc, 0)
        text = output.read_text()
        self.assertEqual(text.count("Sequence_ID"), 1)
        self.assertIn("ASV0", text)
        self.assertIn("ASV1", text)

    def test_retry_policy_handles_connection_read_and_http_errors(self):
        session = MODULE.build_session(4, 1.0, 90, "test")
        self.addCleanup(session.close)
        retry = session.get_adapter("https://").max_retries
        self.assertEqual((retry.total, retry.connect, retry.read), (4, 4, 4))
        self.assertIn(429, retry.status_forcelist)
        self.assertIn(503, retry.status_forcelist)

    def test_empty_and_malformed_results_are_rejected(self):
        for text in ["", "Connection refused", "Error\tDatabase unavailable", "Sequence_ID\thaplo\nmalformed"]:
            with self.subTest(text=text), self.assertRaises(ValueError):
                MODULE.validate_response(text)

    def test_default_user_agent_is_browser_compatible(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            argv = [
                str(SCRIPT),
                "--data-dir", tmp,
                "--output-file", str(Path(tmp) / "out.tsv"),
            ]
            with patch.object(sys, "argv", argv):
                args = MODULE.parse_args()
        self.assertIn("Mozilla/5.0", args.user_agent)


if __name__ == "__main__":
    unittest.main()
