"""Behavioral contracts for biological QC and independent control cohorts."""
import importlib.util
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("control_decontam", ROOT / "processes/control_decontam/control_decontam.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class ControlDecontamTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.out = self.root / "audit"
        self.config = dict(sample_col="Sample", class_col="class", min_biological_reads=5000,
                           labels={c: [c] for c in MODULE.CLASSES},
                           technical=dict(enabled=True, threshold=0.1),
                           bio_control=dict(enabled=True, threshold=0.1))
        self.counts = pd.DataFrame({
            "bio_pass": [5000, 10, 20, 30], "bio_boundary": [5000, 0, 0, 0],
            "bio_fail": [4999, 0, 0, 0], "blank": [0, 10, 0, 5],
            "bio_control": [0, 0, 10, 5], "zero_blank": [0, 0, 0, 0],
            "positive": [5000, 5000, 5000, 5000]}, index=["clean", "tech", "bio_control_asv", "both"])
        MODULE.save_counts(self.counts, self.root / "counts.tsv")
        self.meta = pd.DataFrame(dict(
            Sample=[*self.counts.columns, "missing_zero_scope"],
            **{"class": ["biological"] * 3 + ["technical", "bio_control", "technical", "positive", "bio_control"]},
            Participant_ID=["P1"] * 8, Procedure_ID=["PR1"] * 8,
            custom_pair_id=["PAIR1"] * 8))
        self.meta.to_csv(self.root / "metadata.tsv", sep="\t", index=False)
        pd.DataFrame({"Feature ID": self.counts.index, "Taxon": ["d__Bacteria"] * 4}).to_csv(
            self.root / "taxonomy.tsv", sep="\t", index=False)

    def prepare(self):
        MODULE.prepare(self.root / "counts.tsv", self.root / "metadata.tsv", self.config, self.out)

    def fake_scores(self):
        # Isolate the union operation from the independently tested R scorer.
        for arm, flags in [("TECH", [False, True, False, True]),
                           ("BIO", [False, False, True, True])]:
            pd.DataFrame(dict(ASV_ID=self.counts.index, score=[0.9 if not f else 0.01 for f in flags],
                              contaminant=flags, status="tested")).to_csv(
                                  self.out / f"{arm}_scores.tsv", sep="\t", index=False)

    def test_cutoff_only_biological_and_disjoint_control_cohorts(self):
        self.prepare()
        tech = MODULE.read_counts(self.out / "TECH_counts.tsv")
        bio_control = MODULE.read_counts(self.out / "BIO_counts.tsv")
        self.assertEqual(tech.columns.tolist(), ["bio_pass", "bio_boundary", "blank"])
        self.assertEqual(bio_control.columns.tolist(), ["bio_pass", "bio_boundary", "bio_control"])
        self.assertEqual(MODULE.read_counts(self.out / "positive_qc_counts.tsv").columns.tolist(), ["positive"])
        qc = pd.read_csv(self.out / "sample_qc.tsv", sep="\t").set_index("Sample")
        self.assertEqual(qc.loc["missing_zero_scope", "post_qc_reads"], 0)
        self.assertEqual(qc.loc["bio_control", "custom_pair_id"], "PAIR1")
        self.assertEqual(qc.loc["bio_control", "Procedure_ID"], "PR1")
        prevalence = pd.read_csv(self.out / "prevalence.tsv", sep="\t").set_index("ASV_ID")
        self.assertEqual(prevalence.loc["tech", "prevalence_biological"], 0.5)
        self.assertEqual(prevalence.loc["tech", "n_technical"], 1)

    def test_union_categories_and_no_count_subtraction(self):
        original = (self.root / "counts.tsv").read_bytes()
        self.prepare()
        self.fake_scores()
        MODULE.finalize(self.root / "taxonomy.tsv", self.out)
        calls = pd.read_csv(self.out / "contamination_calls.tsv", sep="\t")
        self.assertEqual(calls.final_category.tolist(), ["CLEAN", "TECH", "BIO", "TECH+BIO"])
        cleaned = MODULE.read_counts(self.out / "ASV_cleaned.tsv")
        pd.testing.assert_frame_equal(cleaned, self.counts.loc[["clean"], ["bio_pass", "bio_boundary"]], check_names=False)
        self.assertEqual((self.root / "counts.tsv").read_bytes(), original)

    def test_independent_arms_and_configurable_depth(self):
        for enabled in ("technical", "bio_control"):
            with self.subTest(enabled=enabled):
                self.out = self.root / enabled
                self.config["technical"]["enabled"] = enabled == "technical"
                self.config["bio_control"]["enabled"] = enabled == "bio_control"
                self.config["min_biological_reads"] = 4999
                self.prepare()
                arm = "TECH" if enabled == "technical" else "BIO"
                other = "BIO" if enabled == "technical" else "TECH"
                self.assertIn("bio_fail", MODULE.read_counts(self.out / f"{arm}_counts.tsv"))
                self.assertFalse((self.out / f"{other}_counts.tsv").exists())
                self.fake_scores()
                MODULE.finalize(self.root / "taxonomy.tsv", self.out)
                calls = pd.read_csv(self.out / "contamination_calls.tsv", sep="\t")
                self.assertTrue((calls[f"{other}_status"] == "disabled").all())
                self.assertFalse(calls[f"{other}_contaminant"].any())

    def test_missing_controls_fail_instead_of_silently_clean(self):
        self.counts["blank"] = 0
        MODULE.save_counts(self.counts, self.root / "counts.tsv")
        with self.assertRaisesRegex(ValueError, "TECH enabled but no nonzero"):
            self.prepare()

    def test_unknown_metadata_class_rejected(self):
        self.meta.loc[0, "class"] = "typo"
        self.meta.to_csv(self.root / "metadata.tsv", sep="\t", index=False)
        with self.assertRaisesRegex(ValueError, "Unmapped sample classes"):
            self.prepare()

    def test_missing_or_ambiguous_metadata_rejected(self):
        for metadata in [self.meta.iloc[1:], pd.concat([self.meta, self.meta.iloc[[0]]])]:
            metadata.to_csv(self.root / "metadata.tsv", sep="\t", index=False)
            with self.assertRaises(ValueError):
                self.prepare()

    def test_metadata_classes_override_any_old_boolean_flags(self):
        self.meta["is_negative_control"] = True
        self.meta["is_positive_control"] = False
        self.meta.to_csv(self.root / "metadata.tsv", sep="\t", index=False)
        self.prepare()
        self.assertNotIn("positive", MODULE.read_counts(self.out / "TECH_counts.tsv"))

    def test_noninteger_counts_rejected(self):
        (self.root / "counts.tsv").write_text("ASV_ID\tS1\nASV1\t0.5\n")
        with self.assertRaisesRegex(ValueError, "integer"):
            MODULE.read_counts(self.root / "counts.tsv")

    def test_duplicate_headers_rejected(self):
        (self.root / "counts.tsv").write_text("ASV_ID\tS1\tS1\nASV1\t1\t2\n")
        with self.assertRaisesRegex(ValueError, "unique sample IDs"):
            MODULE.read_counts(self.root / "counts.tsv")

    @unittest.skipUnless(os.environ.get("ASPIRE_TEST_RSCRIPT"), "Set ASPIRE_TEST_RSCRIPT to test real decontam")
    def test_real_r_prevalence_and_single_control_occurrence(self):
        rscript = os.environ["ASPIRE_TEST_RSCRIPT"]
        # Large enough fixture to distinguish persistent control contamination
        # from an ASV present in all biological samples and a single blank.
        ids = [f"b{i}" for i in range(12)] + [f"n{i}" for i in range(6)]
        counts = pd.DataFrame([
            [100] * 12 + [10, 0, 0, 0, 0, 0],
            [0] * 12 + [10] * 6,
            [100] * 12 + [1] * 6,
            [0] * 18], index=["single_blank", "contaminant", "background", "absent"], columns=ids)
        MODULE.save_counts(counts, self.root / "r_counts.tsv")
        pd.DataFrame(dict(Sample=ids, is_negative=[False] * 12 + [True] * 6)).to_csv(
            self.root / "r_metadata.tsv", sep="\t", index=False)
        subprocess.run([rscript, str(ROOT / "processes/control_decontam/run_prevalence.R"),
                        "--counts", str(self.root / "r_counts.tsv"), "--metadata", str(self.root / "r_metadata.tsv"),
                        "--threshold", "0.1", "--output", str(self.root / "r_scores.tsv")], check=True)
        result = pd.read_csv(self.root / "r_scores.tsv", sep="\t").set_index("ASV_ID")
        self.assertFalse(result.loc["single_blank", "contaminant"])
        self.assertTrue(result.loc["contaminant", "contaminant"])
        self.assertEqual(result.loc["absent", "status"], "absent_in_arm")
        MODULE.save_counts(counts.loc[["background"]], self.root / "r_counts.tsv")
        subprocess.run([rscript, str(ROOT / "processes/control_decontam/run_prevalence.R"),
                        "--counts", str(self.root / "r_counts.tsv"), "--metadata", str(self.root / "r_metadata.tsv"),
                        "--threshold", "0.25", "--output", str(self.root / "r_scores.tsv")], check=True)
        result = pd.read_csv(self.root / "r_scores.tsv", sep="\t")
        self.assertEqual(result.ASV_ID.tolist(), ["background"])
        self.prepare()
        subprocess.run([rscript, str(ROOT / "processes/control_decontam/plot_qc.R"),
                        str(self.out / "sample_qc.tsv"), str(self.out)], check=True)
        self.assertTrue((self.out / "read_depth_by_class.svg").exists())


if __name__ == "__main__":
    unittest.main()
