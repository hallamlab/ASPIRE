from __future__ import annotations

import importlib.util
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).parents[1]
MODULE = ROOT / "processes" / "microbial_compartments"
PREPARE = MODULE / "prepare_microbial_compartments.py"
INFERENCE = MODULE / "infer_microbial_compartments.R"
POSTHOC = MODULE / "posthoc_microbial_compartments.py"
GENERATOR = MODULE / "tests" / "generate_synthetic_microbial_compartments.py"


def load_generator():
    spec = importlib.util.spec_from_file_location("synthetic_microbial_compartments", GENERATOR)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def r_dependencies_available() -> bool:
    if not shutil.which("Rscript"):
        return False
    probe = subprocess.run([
        "Rscript", "-e",
        'quit(status=ifelse(all(vapply(c("cluster","optparse"), requireNamespace, logical(1), quietly=TRUE)),0,1))',
    ], check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return probe.returncode == 0


class MicrobialCompartmentsTest(unittest.TestCase):
    def prepare(self, root: Path, scenario: str = "clear") -> tuple[Path, Path]:
        fixtures = root / "fixtures"
        load_generator().generate(fixtures)
        fixture = fixtures / scenario
        prepared = root / f"prepared_{scenario}"
        subprocess.run([
            sys.executable, str(PREPARE),
            "--counts", str(fixture / "asv_counts.tsv"),
            "--spieceasi-filter-audit", str(fixture / "spieceasi_filtering_audit.csv"),
            "--metadata", str(fixture / "metadata.tsv"),
            "--outdir", str(prepared), "--sample-col", "sampleID", "--transpose",
            "--stability-block-col", "Cruise",
            "--stability-stratum-col", "Season",
            "--zero-replacement", "multiplicative", "--multiplicative-delta", "0",
        ], check=True)
        return fixture, prepared

    def test_multiplicative_replacement_clr_and_aitchison_are_auditable(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            _, prepared = self.prepare(Path(tmp))
            replaced = pd.read_csv(prepared / "asv_zero_replaced_composition.tsv", sep="\t", index_col=0)
            clr = pd.read_csv(prepared / "asv_clr_matrix.tsv", sep="\t", index_col=0)
            distances = pd.read_csv(prepared / "aitchison_distance_matrix.tsv", sep="\t", index_col=0)
            audit = pd.read_csv(prepared / "zero_replacement_audit.tsv", sep="\t")
            config = json.loads((prepared / "microbial_clustering_preparation_config.json").read_text())
            self.assertTrue(np.allclose(replaced.sum(axis=1), 1.0))
            self.assertTrue((replaced > 0).all().all())
            self.assertTrue(np.allclose(clr.mean(axis=1), 0.0, atol=1e-10))
            expected = np.linalg.norm(clr.iloc[0].to_numpy() - clr.iloc[1].to_numpy())
            self.assertAlmostEqual(expected, distances.iloc[0, 1], places=10)
            self.assertTrue(np.allclose(audit.resolved_replacement_value, 1 / replaced.shape[1] ** 2))
            self.assertEqual(config["inference_variables"], ["ASV relative abundance only"])

    @unittest.skipUnless(r_dependencies_available(), "R clustering dependencies are unavailable")
    def test_clear_states_supported_continuous_states_unresolved_and_seed_reproducible(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for scenario, expected_status in (("clear", "supported"), ("continuous", "weak_unresolved")):
                _, prepared = self.prepare(root / scenario, scenario)
                outputs = []
                for run_name in ("run1", "run2"):
                    output = root / scenario / run_name
                    subprocess.run([
                        "Rscript", str(INFERENCE), "--prepared-dir", str(prepared), "--outdir", str(output),
                        "--sample-col", "sampleID", "--k-min", "2", "--k-max", "6",
                        "--min-cluster-size", "5", "--min-cluster-fraction", "0.02",
                        "--min-mean-silhouette", "0.25", "--min-stability-ari", "0.70",
                        "--min-cluster-jaccard", "0.75",
                        "--stability-replicates", "30", "--stability-sample-fraction", "0.8",
                        "--stability-primary-quantile", "0.25", "--stability-near-tie-tolerance", "0.02",
                        "--blocked-robustness-enabled", "TRUE",
                        "--prediction-strength-enabled", "TRUE",
                        "--hierarchical-enabled", "TRUE", "--seed", "42", "--ncpus", "1",
                    ], check=True)
                    outputs.append(output)
                decision = pd.read_csv(outputs[0] / "microbial_cluster_selection_decision.tsv", sep="\t").iloc[0]
                self.assertEqual(decision.support_status, expected_status)
                self.assertEqual(decision.silhouette_selection_role, "descriptive_only")
                self.assertTrue((outputs[0] / "stability_cluster_recovery_summary.tsv").exists())
                self.assertTrue((outputs[0] / "season_balanced_block_replicates.tsv").exists())
                self.assertTrue((outputs[0] / "season_balanced_block_strata_audit.tsv").exists())
                self.assertEqual(
                    (outputs[0] / "microbial_compartments.tsv").read_bytes(),
                    (outputs[1] / "microbial_compartments.tsv").read_bytes(),
                )

    def test_workflow_keeps_inference_and_posthoc_inputs_separate(self) -> None:
        pipeline = (ROOT / "asv_pipeline.nf").read_text()
        workflow = (ROOT / "workflow" / "modules" / "indicators_diagnostics.nf").read_text()
        inference_block = workflow.split("process MICROBIAL_COMPARTMENT_INFERENCE", 1)[1].split(
            "process MICROBIAL_COMPARTMENT_POSTHOC", 1
        )[0]
        self.assertIn("spieceasi_stage.filter_audit.map", pipeline)
        self.assertIn("MICROBIAL_COMPARTMENT_INFERENCE", pipeline)
        self.assertNotIn("metadata_table", inference_block)
        self.assertNotIn("taxonomy_table", inference_block)
        self.assertNotIn("environmental", inference_block.lower())
        self.assertIn("MICROBIAL_COMPARTMENT_POSTHOC", workflow)
        self.assertIn("hierarchical-enabled", inference_block)
        self.assertIn("stability-near-tie-tolerance", inference_block)
        self.assertIn("blocked-robustness-enabled", inference_block)
        self.assertIn("stability-stratum-col", workflow)
        self.assertIn('path("microbial_compartment_posthoc_metadata.tsv"), emit: metadata', workflow)
        self.assertIn("indicspeciesRequiresMicrobialCompartments", pipeline)
        self.assertIn("microbial_compartment_posthoc_stage.metadata", pipeline)


if __name__ == "__main__":
    unittest.main()
