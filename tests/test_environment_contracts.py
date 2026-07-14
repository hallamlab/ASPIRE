from __future__ import annotations

import unittest
from pathlib import Path

import yaml


PROJECT = Path(__file__).parents[1]


class EnvironmentContractTest(unittest.TestCase):
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

    def test_diversity_environment_includes_patient_aware_r_dependencies(self) -> None:
        environment = yaml.safe_load(
            (PROJECT / "processes/diversity_analysis/env.yml").read_text()
        )
        dependencies = {
            item.split("=", 1)[0]
            for item in environment["dependencies"]
            if isinstance(item, str)
        }
        required = {
            "r-base",
            "r-optparse",
            "r-dplyr",
            "r-tidyr",
            "r-vegan",
            "r-permute",
        }
        self.assertTrue(
            required <= dependencies,
            "Patient-aware diversity requires all directly imported R packages.",
        )


if __name__ == "__main__":
    unittest.main()
