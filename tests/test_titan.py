from __future__ import annotations

import importlib.util
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).parents[1]
PREPARE = ROOT / "processes" / "titan" / "prepare_titan.py"
RUN_TITAN = ROOT / "processes" / "titan" / "run_titan.R"
GENERATOR = ROOT / "processes" / "titan" / "tests" / "generate_synthetic_titan_data.py"


def load_generator():
    spec = importlib.util.spec_from_file_location("generate_synthetic_titan_data", GENERATOR)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class TitanModuleTest(unittest.TestCase):
    def prepare(self, root: Path, extra_args: list[str] | None = None) -> Path:
        fixture = root / "fixture"
        prepared = root / "prepared"
        load_generator().generate(fixture)
        command = [
            sys.executable, str(PREPARE),
            "--counts", str(fixture / "asv_counts.tsv"),
            "--filter-audit", str(fixture / "spieceasi_filtering_audit.csv"),
            "--measurement-matrix", str(fixture / "measurement_matrix.samples_by_measurement.tsv"),
            "--taxonomy", str(fixture / "taxonomy.tsv"),
            "--outdir", str(prepared),
            "--sample-col", "sampleID",
            "--asv-id-col", "ASV_ID",
            "--variables", "synthetic_gradient|invariant_gradient|sparse_gradient",
            "--transpose", "--minimum-occurrence", "3",
            "--minimum-prevalence", "0", "--minimum-mean-relative-abundance", "0",
        ]
        command.extend(extra_args or [])
        subprocess.run(command, check=True)
        return prepared

    def test_preparation_uses_exact_spieceasi_cohort_and_variablewise_screens(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            prepared = self.prepare(Path(tmp))
            manifest = pd.read_csv(prepared / "titan_input_manifest.tsv", sep="\t")
            status = dict(zip(manifest.environmental_variable, manifest.status))
            self.assertEqual(status["synthetic_gradient"], "ready")
            self.assertEqual(status["invariant_gradient"], "skipped")
            self.assertEqual(status["sparse_gradient"], "skipped")
            cohort = pd.read_csv(prepared / "titan_spieceasi_cohort.tsv", sep="\t")
            self.assertEqual(
                cohort.ASV_ID.tolist(),
                ["ASV_declining", "ASV_increasing", "ASV_background", "ASV_rare"],
            )
            sample_audit = pd.read_csv(prepared / "titan_sample_matching_audit.tsv", sep="\t")
            self.assertEqual(len(sample_audit), 24)
            self.assertTrue(sample_audit.matched_by_explicit_identifier.all())
            relative = pd.read_csv(
                prepared / "variables" / "synthetic_gradient" / "taxa_relative_abundance.tsv",
                sep="\t",
            )
            self.assertNotIn("ASV_rare", relative.columns)
            self.assertNotIn("ASV_filtered", relative.columns)
            retained_mass = relative.drop(columns="sampleID").sum(axis=1)
            self.assertTrue(retained_mass.le(1.0 + 1e-12).all())
            self.assertAlmostEqual(retained_mass.iloc[0], 605.0 / 625.0)
            self.assertAlmostEqual(retained_mass.iloc[-1], 1.0)
            audit = pd.read_csv(prepared / "titan_filtered_asv_audit.tsv", sep="\t")
            rare = audit.loc[
                audit.environmental_variable.eq("synthetic_gradient") & audit.identifier.eq("ASV_rare")
            ].iloc[0]
            self.assertEqual(rare.exclusion_reason, "occurrence_below_titan_minimum")

    def test_optional_future_study_abundance_filter_is_audited(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            prepared = self.prepare(
                Path(tmp), ["--minimum-mean-relative-abundance", "0.20"]
            )
            relative = pd.read_csv(
                prepared / "variables" / "synthetic_gradient" / "taxa_relative_abundance.tsv",
                sep="\t",
            )
            self.assertNotIn("ASV_background", relative.columns)
            audit = pd.read_csv(prepared / "titan_filtered_asv_audit.tsv", sep="\t")
            background = audit.loc[
                audit.environmental_variable.eq("synthetic_gradient")
                & audit.identifier.eq("ASV_background")
            ].iloc[0]
            self.assertIn("mean_relative_abundance_below_configured_minimum", background.exclusion_reason)

    @unittest.skipUnless(shutil.which("Rscript"), "Rscript is unavailable")
    def test_fixed_seed_reproducibility_and_known_transition(self) -> None:
        probe = subprocess.run(
            ["Rscript", "-e", 'quit(status=ifelse(requireNamespace("TITAN2", quietly=TRUE),0,1))'],
            check=False,
        )
        if probe.returncode != 0:
            self.skipTest("TITAN2 is unavailable")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            prepared = self.prepare(root)
            outputs = []
            for name in ("run1", "run2"):
                out = root / name
                subprocess.run([
                    "Rscript", str(RUN_TITAN), "--prepared-dir", str(prepared),
                    "--outdir", str(out), "--permutations", "50",
                    "--bootstrap-count", "50", "--seed", "42", "--ncpus", "1",
                ], check=True)
                outputs.append(out)
            first = outputs[0] / "variables" / "synthetic_gradient" / "titan_taxon_results.tsv"
            second = outputs[1] / "variables" / "synthetic_gradient" / "titan_taxon_results.tsv"
            self.assertEqual(first.read_bytes(), second.read_bytes())
            result = pd.read_csv(first, sep="\t")
            directions = dict(zip(result.ASV_ID, result.response_direction))
            self.assertEqual(directions["ASV_declining"], "z-")
            self.assertEqual(directions["ASV_increasing"], "z+")

    def test_workflow_contract(self) -> None:
        pipeline = (ROOT / "asv_pipeline.nf").read_text()
        module = (ROOT / "workflow" / "modules" / "indicators_diagnostics.nf").read_text()
        config = (ROOT / "examples" / "si.local.yml").read_text()
        self.assertIn("retained_final", PREPARE.read_text())
        self.assertIn("spieceasi_stage.filter_audit.map", pipeline)
        self.assertIn("process TITAN_PREPARE", module)
        self.assertIn("process TITAN_ANALYSIS", module)
        self.assertIn("process TITAN_COLLECT", module)
        self.assertIn("process TITAN_PLOTS", module)
        self.assertIn("environmental_variables: []", config)


if __name__ == "__main__":
    unittest.main()
