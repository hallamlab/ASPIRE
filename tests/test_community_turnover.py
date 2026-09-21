from __future__ import annotations

import importlib.util
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).parents[1]
PREPARE = ROOT / "processes" / "microbial_compartments" / "prepare_microbial_compartments.py"
TURNOVER = ROOT / "processes" / "community_turnover" / "community_turnover.py"
LCBD = ROOT / "processes" / "community_turnover" / "lcbd_analysis.R"
GENERATOR = ROOT / "processes" / "community_turnover" / "tests" / "generate_synthetic_turnover_data.py"


def load_generator():
    spec = importlib.util.spec_from_file_location("synthetic_turnover", GENERATOR)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def adespatial_available() -> bool:
    if not shutil.which("Rscript"):
        return False
    result = subprocess.run(
        ["Rscript", "-e", 'quit(status=ifelse(all(vapply(c("adespatial","optparse"), requireNamespace, logical(1), quietly=TRUE)),0,1))'],
        check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    return result.returncode == 0


class CommunityTurnoverTest(unittest.TestCase):
    def run_scenario(self, root: Path, scenario: str) -> tuple[Path, Path, Path]:
        fixtures = root / "fixtures"
        if not fixtures.exists():
            load_generator().generate(fixtures)
        fixture = fixtures / scenario
        prepared = root / f"prepared_{scenario}"
        output = root / f"output_{scenario}"
        subprocess.run([
            sys.executable, str(PREPARE), "--counts", str(fixture / "asv_counts.tsv"),
            "--spieceasi-filter-audit", str(fixture / "spieceasi_filtering_audit.csv"),
            "--metadata", str(fixture / "metadata.tsv"), "--outdir", str(prepared),
            "--sample-col", "sampleID", "--transpose", "--zero-replacement", "multiplicative",
        ], check=True)
        subprocess.run([
            sys.executable, str(TURNOVER), "--prepared-dir", str(prepared),
            "--metadata", str(fixture / "metadata.tsv"), "--outdir", str(output),
            "--sample-col", "sampleID", "--profile-col", "Cruise", "--date-col", "date",
            "--depth-col", "Depth", "--environmental-compartment-cols",
            "o2_compartment|gmm_component|o2_subcompartment_final",
            "--primary-environmental-compartment-col", "o2_subcompartment_final",
            "--distance-metrics", "aitchison,braycurtis",
            "--fixed-depth-min-profile-fraction", "0.5",
        ], check=True)
        return fixture, prepared, output

    def test_known_depth_transition_and_irregular_spacing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            _, _, output = self.run_scenario(Path(tmp), "sharp_depth")
            vertical = pd.read_csv(output / "tables" / "vertical_turnover.tsv", sep="\t")
            primary = vertical.loc[vertical.distance_metric.eq("aitchison")]
            peaks = primary.loc[primary.groupby("profile_ID").community_distance.idxmax()]
            self.assertTrue(peaks.shallow_depth.eq(50).all())
            expected_deep = {"P1": 100, "P2": 100, "P3": 150, "P4": 100}
            self.assertEqual(dict(zip(peaks.profile_ID, peaks.deep_depth)), expected_deep)
            self.assertTrue((primary.turnover_per_meter * primary.depth_difference).sub(primary.community_distance).abs().lt(1e-10).all())
            depth_audit = pd.read_csv(output / "audit" / "fixed_depth_prevalence_audit.tsv", sep="\t")
            depth100 = depth_audit.loc[depth_audit.depth.eq(100)].iloc[0]
            self.assertAlmostEqual(depth100.profile_fraction, 0.75)
            self.assertTrue(depth100.retained_as_common_depth)
            duplicate_audit = pd.read_csv(output / "audit" / "duplicate_profile_depth_centers.tsv", sep="\t")
            self.assertEqual(len(duplicate_audit), 1)
            self.assertEqual(duplicate_audit.iloc[0].number_of_member_samples, 2)
            boundary_summary = pd.read_csv(
                output / "tables" / "vertical_turnover_boundary_summary.tsv", sep="\t"
            )
            boundary_cruises = pd.read_csv(
                output / "tables" / "vertical_turnover_boundary_cruise_contrasts.tsv", sep="\t"
            )
            self.assertEqual(set(boundary_summary.compartment_strategy), {"O2", "GMM", "Hybrid"})
            self.assertTrue((boundary_summary.paired_profile_n >= 0).all())
            self.assertTrue({
                "mean_paired_turnover_difference", "mean_paired_difference_ci95_low",
                "mean_paired_difference_ci95_high", "paired_sign_flip_pvalue",
                "paired_sign_flip_qvalue",
            }.issubset(boundary_summary.columns))
            inferential = boundary_summary.loc[boundary_summary.paired_profile_n.gt(0)]
            self.assertTrue(inferential.paired_sign_flip_pvalue.between(0, 1).all())
            self.assertTrue(inferential.paired_sign_flip_qvalue.between(0, 1).all())
            self.assertTrue(
                np.allclose(
                    boundary_cruises.paired_turnover_difference,
                    boundary_cruises.median_crossing_turnover
                    - boundary_cruises.median_within_turnover,
                )
            )

    def test_missing_dates_are_excluded_and_audited(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixtures = root / "fixtures"
            load_generator().generate(fixtures)
            fixture = fixtures / "sharp_depth"
            metadata = pd.read_csv(fixture / "metadata.tsv", sep="\t")
            metadata.loc[metadata.sampleID.eq("P2_50m"), "date"] = ""
            metadata.to_csv(fixture / "metadata.tsv", sep="\t", index=False)
            _, _, output = self.run_scenario(root, "sharp_depth")
            audit = pd.read_csv(output / "audit" / "community_turnover_sample_audit.tsv", sep="\t")
            excluded = audit.loc[audit.sampleID.eq("P2_50m")].iloc[0]
            self.assertFalse(excluded.retained_for_turnover)
            self.assertEqual(excluded.turnover_exclusion_reason, "missing_or_invalid_date")

    def test_moving_boundary_and_within_compartment_change_are_distinguishable(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, _, moving_output = self.run_scenario(root, "moving_boundary")
            moving_fixed = pd.read_csv(moving_output / "tables" / "temporal_fixed_depth_turnover.tsv", sep="\t")
            moving_comp = pd.read_csv(moving_output / "tables" / "temporal_environmental_compartment_turnover.tsv", sep="\t")
            moving_fixed = moving_fixed.loc[moving_fixed.distance_metric.eq("aitchison")]
            moving_comp = moving_comp.loc[
                moving_comp.distance_metric.eq("aitchison") & moving_comp.environmental_variable.eq("o2_compartment")
            ]
            self.assertGreater(moving_fixed.community_distance.max(), moving_comp.community_distance.max())

            _, _, change_output = self.run_scenario(root, "within_compartment_change")
            changed = pd.read_csv(change_output / "tables" / "temporal_environmental_compartment_turnover.tsv", sep="\t")
            changed = changed.loc[
                changed.distance_metric.eq("aitchison") & changed.environmental_variable.eq("o2_compartment")
            ]
            self.assertGreater(changed.community_distance.max(), moving_comp.community_distance.max())
            self.assertTrue(changed.time_difference_days.gt(0).all())
            temporal_summary = pd.read_csv(
                change_output / "tables" / "temporal_compartment_turnover_summary.tsv", sep="\t"
            )
            self.assertFalse(temporal_summary.empty)

    @unittest.skipUnless(adespatial_available(), "adespatial is unavailable")
    def test_lcbd_global_and_within_profile_schemas(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fixture, prepared, _ = self.run_scenario(Path(tmp), "sharp_depth")
            out = Path(tmp) / "lcbd"
            subprocess.run([
                "Rscript", str(LCBD), "--prepared-dir", str(prepared),
                "--metadata", str(fixture / "metadata.tsv"), "--outdir", str(out),
                "--sample-col", "sampleID", "--profile-col", "Cruise", "--date-col", "date",
                "--depth-col", "Depth", "--environmental-compartment-cols", "o2_compartment",
                "--permutations", "49", "--seed", "42", "--within-profile-enabled", "TRUE",
            ], check=True)
            lcbd = pd.read_csv(out / "lcbd.tsv", sep="\t")
            self.assertEqual(set(lcbd.reference_population), {"global", "within_profile"})
            self.assertAlmostEqual(lcbd.loc[lcbd.reference_population.eq("global"), "LCBD"].sum(), 1.0)
            self.assertTrue(lcbd.loc[lcbd.reference_population.eq("within_profile"), "LCBD_pvalue"].isna().all())
            parameters = pd.read_csv(out / "lcbd_parameters.tsv", sep="\t").iloc[0]
            self.assertGreaterEqual(parameters.adespatial_input_translation, 0)
            self.assertLess(parameters.direct_CLR_LCBD_max_abs_difference, 1e-10)

    def test_workflow_contract(self) -> None:
        pipeline = (ROOT / "asv_pipeline.nf").read_text()
        workflow = (ROOT / "workflow" / "modules" / "indicators_diagnostics.nf").read_text()
        config = (ROOT / "examples" / "si.local.yml").read_text()
        for process in ("COMMUNITY_TURNOVER_PREPARE", "COMMUNITY_TURNOVER_ANALYSIS", "COMMUNITY_TURNOVER_LCBD", "COMMUNITY_TURNOVER_PLOTS"):
            self.assertIn(f"process {process}", workflow)
            self.assertIn(process, pipeline)
        self.assertIn("spieceasi_stage.filter_audit.map", pipeline)
        self.assertIn("fixed_depth_min_profile_fraction: 0.50", config)
        self.assertIn("primary_environmental_compartment_col: o2_subcompartment_final", config)


if __name__ == "__main__":
    unittest.main()
