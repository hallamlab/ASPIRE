from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "processes" / "asv_mag_curtains" / "asv_mag_curtains.py"
SPEC = importlib.util.spec_from_file_location("asv_mag_curtains", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_nextflow_wrapper_escapes_shell_environment_variables() -> None:
    wrapper = (ROOT / "workflow/modules/networks_reporting.nf").read_text(
        encoding="utf-8"
    )
    process = wrapper[
        wrapper.index("process ASV_MAG_CURTAINS") : wrapper.index(
            "process GROUP_GUILD_FUNCTION"
        )
    ]
    assert 'export MPLCONFIGDIR="\\$PWD/.matplotlib"' in process
    assert 'mkdir -p "\\$MPLCONFIGDIR"' in process


def test_curtain_environment_declares_spearman_dependency() -> None:
    environment = (ROOT / "processes/asv_mag_curtains/env.yml").read_text(
        encoding="utf-8"
    )
    assert "  - scipy\n" in environment


class AsvMagCurtainTest(unittest.TestCase):
    def test_renewal_onsets_are_loaded_and_drawn_as_dashed_lines(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "renewals.tsv"
            pd.DataFrame({"start_date": ["2010-01-02", "2010-01-02", "2011-03-04"]}).to_csv(
                path, sep="\t", index=False
            )
            onsets = MODULE.load_renewal_onsets(path, "start_date")
        self.assertEqual(len(onsets), 2)
        figure, axis = MODULE.plt.subplots()
        MODULE.draw_renewal_onsets(
            axis, onsets, pd.Timestamp("2010-01-01"), pd.Timestamp("2010-12-31")
        )
        self.assertEqual(len(axis.lines), 1)
        self.assertEqual(axis.lines[0].get_linestyle(), "--")
        MODULE.plt.close(figure)

    def test_species_aggregation_uses_mean_across_linked_genomes(self) -> None:
        profiles = pd.DataFrame({
            "ASV_ID": ["ASV1", "ASV1", "ASV1"],
            "species": ["SUP05 sp1", "SUP05 sp1", "Other sp2"],
            "mag_match_id": ["MAG1", "MAG2", "MAG3"],
            "modality": ["metagenome"] * 3,
            "sample": ["S1"] * 3,
            "date": pd.to_datetime(["2020-01-01"] * 3),
            "depth_m": [100.0] * 3,
            "normalization": ["input_fragment_fpm"] * 3,
            "normalized_recruitment": [10.0, 30.0, 500.0],
            "raw_recruitment_count": [100.0, 300.0, 5000.0],
        })
        eligible = pd.DataFrame({
            "ASV_ID": ["ASV1", "ASV1"],
            "species": ["SUP05 sp1", "SUP05 sp1"],
            "mag_match_id": ["MAG1", "MAG2"],
        })
        result = MODULE.aggregate_species_points(profiles, eligible)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["linked_genome_count"], 2)
        self.assertAlmostEqual(result.iloc[0]["mean_normalized_recruitment_fpm"], 20.0)

    def test_point_area_is_linear_and_strongly_separates_abundance(self) -> None:
        sizes = MODULE.size_scale(
            pd.Series([100.0, 10_000.0, 100_000.0]), 100_000.0, 1.0, 320.0
        )
        self.assertTrue(np.all(np.diff(sizes) > 0))
        self.assertAlmostEqual(sizes[0], 1.0)
        self.assertAlmostEqual(sizes[1], 32.0)
        self.assertAlmostEqual(sizes[-1], 320.0)

    def test_pairwise_contrast_uses_matched_samples_and_documented_ratio(self) -> None:
        points = pd.DataFrame({
            "ASV_ID": ["ASV1"] * 4,
            "species": ["A", "B", "A", "B"],
            "modality": ["metagenome"] * 4,
            "sample": ["S1", "S1", "S2", "S2"],
            "date": pd.to_datetime(["2020-01-01", "2020-01-01", "2020-02-01", "2020-02-01"]),
            "depth_m": [100.0, 100.0, 150.0, 150.0],
            "normalization": ["input_fragment_fpm"] * 4,
            "mean_normalized_recruitment_fpm": [9.0, 4.0, 1.0, 7.0],
        })
        result = MODULE.calculate_species_contrasts(points, pseudocount_fpm=1.0, fold_threshold=2.0)
        self.assertEqual(len(result), 2)
        first = result.loc[result["sample"].eq("S1")].iloc[0]
        self.assertAlmostEqual(first["log2_fpm_ratio"], 1.0)
        self.assertEqual(first["contrast_class"], "numerator_at_least_threshold")
        second = result.loc[result["sample"].eq("S2")].iloc[0]
        self.assertEqual(second["contrast_class"], "denominator_at_least_threshold")

        summary = MODULE.summarize_species_contrasts(result)
        self.assertEqual(summary.iloc[0]["matched_samples_n"], 2)
        self.assertAlmostEqual(abs(summary.iloc[0]["spearman_rho"]), 1.0)

    def test_sample_anchors_use_explicit_ids_and_surface_cohort(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "crosswalk.tsv"
            pd.DataFrame({
                "sampleID": ["S1", "S2", "S3", "S3"],
                "Date": ["2020-01-01", "2020-02-01", "2020-03-01", "2020-03-01"],
                "Depth": [10, 20, 30, 30],
                "hybrid": ["A", "B", None, None],
            }).to_csv(path, sep="\t", index=False)
            result = MODULE.load_sample_anchors(
                path, "sampleID", "Date", "Depth", "hybrid", ["A"], "hybrid"
            )
        self.assertEqual(result["sample_ID"].tolist(), ["S1"])
        self.assertEqual(result["state"].tolist(), ["A"])
        self.assertEqual(result["state_surface"].tolist(), ["hybrid"])

    def test_workflow_wires_curtains_after_both_required_stages(self) -> None:
        pipeline = (ROOT / "asv_pipeline.nf").read_text()
        processes = (ROOT / "workflow" / "modules" / "networks_reporting.nf").read_text()
        self.assertIn("process ASV_MAG_CURTAINS", processes)
        self.assertIn("asv_mag_network_stage.done", pipeline)
        self.assertIn("microbial_state_interpretation_stage.done", pipeline)
        # Script-global resources are consumed directly by the process and
        # must not be referenced from the method-scoped configuration parser.
        self.assertNotIn(
            "asvMagCurtainsCondaEnvPath: asvMagCurtainsCondaEnvPath", pipeline
        )
        self.assertNotIn("asvMagCurtainsScriptPath: asvMagCurtainsScriptPath", pipeline)
        self.assertNotIn("asvMagCurtainsScriptHash: asvMagCurtainsScriptHash", pipeline)
        self.assertIn("asvMagCurtainsRenewalEvents", pipeline)
        self.assertIn("--renewal-events", processes)


if __name__ == "__main__":
    unittest.main()
