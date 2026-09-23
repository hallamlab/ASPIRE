from __future__ import annotations

import unittest
from pathlib import Path

import yaml


PROJECT = Path(__file__).parents[1]


class EnvironmentContractTest(unittest.TestCase):
    def test_controller_pins_supported_nextflow_release(self) -> None:
        environment = yaml.safe_load(
            (PROJECT / "processes/controller/env.yml").read_text()
        )
        self.assertIn("nextflow=25.10.4", environment["dependencies"])

    def test_outlier_environment_includes_plotting_dependencies(self) -> None:
        environment = yaml.safe_load(
            (PROJECT / "processes/shared_envs/outlier_checker.yml").read_text()
        )
        dependencies = {
            item.split("=", 1)[0]
            for item in environment["dependencies"]
            if isinstance(item, str)
        }
        self.assertTrue(
            {"matplotlib", "seaborn"} <= dependencies,
            "The outlier checker writes plots and must include its plotting libraries.",
        )

    def test_diversity_environment_supports_python_and_r_steps(self) -> None:
        environment = yaml.safe_load(
            (PROJECT / "processes/diversity_analysis/env.yml").read_text()
        )
        dependencies = {
            item.split("=", 1)[0]
            for item in environment["dependencies"]
            if isinstance(item, str)
        }
        self.assertTrue(
            {"r-base", "r-optparse", "r-vegan", "r-permute"} <= dependencies,
            "Diversity analysis runs the Python plots and an R vegan PERMANOVA step.",
        )

    def test_batch_correction_uses_staged_nextflow_inputs(self) -> None:
        workflow = (PROJECT / "asv_pipeline.nf").read_text()
        process_body = workflow.split("process ASV_BATCH_CORRECTION {", 1)[1].split(
            "process ASV_META_FROM_CORRECTED {", 1
        )[0]
        self.assertIn('--data-dir "."', process_body)
        self.assertIn('--asv "${asv_counts}"', process_body)
        self.assertIn('--metadata "${metadata_table}"', process_body)
        self.assertIn('--asv-meta "${asv_meta}"', process_body)
        self.assertNotIn('"ASVs/${asv_counts}"', process_body)


if __name__ == "__main__":
    unittest.main()
