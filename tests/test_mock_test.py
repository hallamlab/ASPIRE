from __future__ import annotations

import hashlib
import importlib.util
import gzip
import tempfile
import unittest
from pathlib import Path

import yaml


SCRIPT = Path(__file__).parents[1] / "examples/mock_test/configure_mock_run.py"
SPEC = importlib.util.spec_from_file_location("configure_mock_run", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)

VALIDATOR_SCRIPT = Path(__file__).parents[1] / "examples/mock_test/validate_results.py"
VALIDATOR_SPEC = importlib.util.spec_from_file_location("validate_results", VALIDATOR_SCRIPT)
VALIDATOR = importlib.util.module_from_spec(VALIDATOR_SPEC)
assert VALIDATOR_SPEC.loader is not None
VALIDATOR_SPEC.loader.exec_module(VALIDATOR)


class MockTestConfigTest(unittest.TestCase):
    def test_build_config_replaces_all_fixture_paths(self) -> None:
        project = Path(__file__).parents[1]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            dataset = root / "fixture"
            output = root / "output"
            runtime = output / ".aspire"
            config = MODULE.build_config(
                project / "examples/mock.local.yml", dataset, output, runtime, project
            )
            primers = MODULE.config_section(config, "primer_trimming")
            self.assertTrue(primers['enabled'])
            self.assertIn('Class:Mammalia', MODULE.config_section(config, 'filter_counts')['exclude_taxa'])
            self.assertEqual(primers['primers'][0]['name'], '515F_806R')
            self.assertTrue(all(v == 0 for v in MODULE.config_section(config, 'fastp').values()))
            paths = MODULE.config_section(config, "paths")
            self.assertEqual(paths["input_dir"], str(dataset / "fastq"))
            self.assertEqual(paths["manifest"], str(dataset / "fastq_manifest.tsv"))
            self.assertEqual(MODULE.config_section(config, "metadata_plots")["metadata"], str(dataset / "sample_metadata.tsv"))
            self.assertEqual(MODULE.config_section(config, "voc_correlation")["voc_table"], str(dataset / "chemistry.tsv"))
            voc_config = MODULE.config_section(config, "voc_correlation")
            self.assertTrue(voc_config["patient_inference"])
            self.assertEqual(voc_config["patient_permutations"], 999)
            self.assertEqual(voc_config["clr_pseudocount"], 0.5)
            self.assertEqual(paths["output_dir"], str(output))
            self.assertEqual(paths["runtime_dir"], str(runtime))
            self.assertTrue(paths["keep_runtime_dir"])
            self.assertTrue(Path(MODULE.config_section(config, "table_filter")["script"]).is_absolute())
            self.assertTrue(all(Path(path).is_absolute() for path in config["environments"].values()))

    def test_checksum_validation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            dataset = Path(tmp)
            payload = dataset / "payload.tsv"
            payload.write_text("fixture\n")
            digest = hashlib.sha256(payload.read_bytes()).hexdigest()
            (dataset / "checksums.tsv").write_text(
                f"relative_path\tsha256\npayload.tsv\t{digest}\n"
            )
            MODULE.validate_checksums(dataset)
            payload.write_text("changed\n")
            with self.assertRaises(SystemExit):
                MODULE.validate_checksums(dataset)

    def test_generated_yaml_is_serializable(self) -> None:
        project = Path(__file__).parents[1]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = MODULE.build_config(
                project / "examples/mock.local.yml",
                root / "fixture",
                root / "output",
                root / "runtime",
                project,
            )
            reparsed = yaml.safe_load(yaml.safe_dump(config, sort_keys=False))
            self.assertEqual(MODULE.config_section(reparsed, "indicspecies")["perms"], 999)
            self.assertEqual(MODULE.config_section(reparsed, "power_analysis")["sample_sizes_cancer"], "4,6,8,10,15,20,30,40,50")
            self.assertEqual(MODULE.config_section(reparsed, "power_analysis")["n_simulations"], 1000)
            self.assertEqual(MODULE.config_section(reparsed, "power_analysis")["n_perm"], 999)

    def test_truth_mapping_accepts_trimmed_inferred_sequences(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            dataset = root / "dataset"
            results = root / "results"
            dataset.mkdir()
            fasta = results / "intermediates/ASVs/ASVs.fasta.gz"
            fasta.parent.mkdir(parents=True)
            (dataset / "ground_truth_feature_registry.tsv").write_text(
                "ASV_ID\tv4_sequence\ntruth_1\tAAAACCCCGGGGTTTT\n"
            )
            with gzip.open(fasta, "wt") as handle:
                handle.write(">ASV1;size=20\nCCCCGGGG\n")
            self.assertEqual(
                VALIDATOR.map_inferred_to_truth(dataset, results), {"ASV1": "truth_1"}
            )

    def test_portable_manifest_recovers_stale_fastq_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            dataset = Path(tmp) / "mock_dataset"
            fastq = dataset / "fastq"
            fastq.mkdir(parents=True)
            (fastq / "S1_R1.fastq.gz").write_bytes(b"fixture")
            (fastq / "S1_R2.fastq.gz").write_bytes(b"fixture")
            (dataset / "fastq_manifest.tsv").write_text(
                "sample_id\tfastq_r1\tfastq_r2\n"
                "S1\t/old/machine/S1_R1.fastq.gz\t/old/machine/S1_R2.fastq.gz\n"
            )
            destination = Path(tmp) / "mock_run.manifest.tsv"
            path, ids = MODULE.write_portable_manifest(dataset, destination)
            content = path.read_text()
            self.assertEqual(ids, {"S1"})
            self.assertIn(str((fastq / "S1_R1.fastq.gz").resolve()), content)
            self.assertIn(str((fastq / "S1_R2.fastq.gz").resolve()), content)


if __name__ == "__main__":
    unittest.main()


class CamiShowcaseTests(unittest.TestCase):
    def test_cami_profile_uses_body_sites_and_established_colors(self):
        import pandas as pd
        project=Path(__file__).parents[1]
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);dataset=root/'dataset';dataset.mkdir()
            pd.DataFrame([
                dict(sample_id=site+'1',source_sample=site+'_source',body_site=site,
                     study_role='biological_control' if site=='Skin' else 'comparison',
                     Type_Group=site,batch='plate_1',DNA_conc=20,
                     is_negative_control=False,is_positive_control=False)
                for site in ('Airways','Oral','Skin')
            ]).to_csv(dataset/'sample_metadata.tsv',sep='\t',index=False)
            pd.DataFrame({'sample_id':['Airways1','Oral1','Skin1'],'acetone':[1.,2.,3.]}).to_csv(dataset/'chemistry.tsv',sep='\t',index=False)
            MODULE.validate_tabular_inputs(dataset,{'Airways1','Oral1','Skin1'})
            cfg=MODULE.build_config(project/'examples/mock.local.yml',dataset,root/'output',root/'runtime',project,4)
            section=lambda name: MODULE.config_section(cfg,name)
            self.assertEqual(section('table_filter')['min_sample_reads'],5000)
            self.assertEqual(section('filter_counts')['min_prevalence_fraction'],0.05)
            self.assertEqual(section('filter_counts')['min_relative_abundance_pct'],0.1)
            self.assertNotIn('subtraction_groups', section('metadata_plots'))
            self.assertEqual(section('metadata_plots')['keep_types'],['Airways','Oral'])
            self.assertFalse(section('voc_correlation')['patient_inference'])
            self.assertFalse(section('diversity')['patient_aware']['enabled'])
            for name in ('power_analysis','lung_status_analysis','taxonomy_patient_aware'):
                self.assertFalse(section(name)['enabled'])
            self.assertTrue(section('control_decontam')['enabled'])
            self.assertEqual(section('control_decontam')['biological_labels'], ['Airways','Oral'])
            self.assertEqual(section('control_decontam')['bio_control_labels'], ['Skin'])
            self.assertEqual(section('control_decontam')['technical_labels'], ['Control'])
            self.assertTrue(section('control_decontam')['bio_control_enabled'])
            self.assertTrue(section('control_decontam')['technical_enabled'])
            self.assertEqual(section('batch_correction')['biological_covariates'],'Type_Group')
            self.assertIn('Oral=#6A3D9A',section('diversity')['group1_palette'])
            palette=pd.read_csv(section('metadata_plots')['palette_file'],sep='\t').set_index('value').color.to_dict()
            self.assertEqual(palette['Skin'],'#CC79A7')
            self.assertEqual(palette['Airways'],'#009E73')
            serialized=yaml.safe_dump(cfg)
            for stale in ('Bronchial Brush','Cancer=#','Control=#8C8C8C,Cancer','biological_covariates: Type_Group,Case'):
                self.assertNotIn(stale,serialized)

    def test_clinical_cami_enables_patient_and_lung_modules(self):
        project=Path(__file__).parents[1]
        cfg=yaml.safe_load((project/'examples/mock.local.yml').read_text())
        import pandas as pd
        clinical=MODULE.has_synthetic_clinical_design(pd.DataFrame({'synthetic_clinical':[True,False]}))
        self.assertIs(type(clinical),bool)
        cfg=MODULE.configure_cami_profile(cfg,project,clinical=clinical)
        yaml.safe_dump(cfg)
        section=lambda name:MODULE.config_section(cfg,name)
        for name in ('power_analysis','lung_status_analysis','taxonomy_patient_aware'):
            self.assertTrue(section(name)['enabled'])
        self.assertTrue(section('diversity')['patient_aware']['enabled'])
        self.assertTrue(section('voc_correlation')['patient_inference'])
        self.assertEqual(section('voc_correlation')['patient_col'],'Participant_ID')
        self.assertEqual(section('lung_status_analysis')['sample_types'],'Airways')
        self.assertEqual(section('lung_status_analysis')['case_col'],'Case')
        self.assertEqual(section('diversity')['patient_aware']['contralateral_value'],'Contralateral')
        self.assertEqual(section('metadata_plots')['keep_types'],['Airways','Oral'])
        self.assertNotIn('subtraction_groups', section('metadata_plots'))
        self.assertEqual(section('indicspecies')['group_orders']['Case'],['Control','Cancer'])
        self.assertIn('Oral=#6A3D9A',section('diversity')['group1_palette'])
