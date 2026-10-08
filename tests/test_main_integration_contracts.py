"""Contracts protecting both local and upstream functionality during integration."""
import os
import re
import subprocess
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def test_optional_graph_routes_new_modules_through_filtered_tables():
    source = (ROOT / 'asv_pipeline.nf').read_text()
    body = source.split('workflow optional {', 1)[1].split('workflow RUN_METADATA_ANALYSES {', 1)[0]
    assert 'CONTROL_DECONTAM(' not in body
    core = source.split('workflow core {', 1)[1].split('workflow standard {', 1)[0]
    assert core.index('TAXONOMY(') < core.index('CONTROL_DECONTAM(')
    assert 'FILTER_TABLE(' not in source and 'FILTER_COUNTS(' not in source
    standard = source.split('workflow standard {', 1)[1].split('workflow optional {', 1)[0]
    assert standard.index('MITO_DECONTAM(') < standard.index('FILTER_ASVS(') < standard.index('PLOT_METADATA(')
    assert 'def metadataMicroInput = filter_counts_stage.filtered_counts' in standard
    assert 'filtered_fasta = filter_counts_stage.filtered_fasta' in standard
    assert 'counts_for_feature_filter = decontam_stage.cleaned' in core
    assert "binding.setVariable('metadataPlotsSubtractionGroups', [])" in source
    assert 'controlDecontamEnabled ? 0 :' in source
    for name in ['Batch', 'MeasurementAssociation', 'GroupingDiagnostics']:
        assert f'asvFinalFor{name} = baseAsvFinal.map' in body
    assert 'asvMetaForGroupAugmentation = baseAsvMeta.map' in body
    assert 'asvMetaForMeasurementAssociation = baseAsvMeta.map' in body
    assert 'batch_stage.asv_selected_counts_int' in body
    assert 'taxonomy_stage.' not in body
    assert 'metadata_analysis_stage.' not in body
    for process in ['GROUPING_DIAGNOSTICS', 'GROUP_LABEL_AUGMENTATION',
                    'MEASUREMENT_ASSOCIATION', 'ASV_MAG_NETWORK',
                    'GROUP_POWER_ANALYSIS', 'TAXONOMY_GROUP_ASSOCIATION',
                    'PAIRED_GROUP_CONTRAST', 'VOC_CORRELATION', 'NETWORK_TOPOLOGY']:
        assert re.search(r'\b' + process + r'\(', body), process


def test_launcher_lists_every_stage_with_tier_and_preserves_aliases():
    env = dict(os.environ, IN_CONTROLLER_ENV='1')
    result = subprocess.run(['bash', str(ROOT / 'run_asv_pipeline.sh'), '--list-stages'],
                            env=env, check=True, capture_output=True, text=True)
    stages = result.stdout.split('# aliases')[0].splitlines()
    assert stages and all(re.fullmatch(r'(core|standard|optional):[A-Z_]+', row) for row in stages)
    assert 'POWER_ANALYSIS_PIPELINE -> GROUP_POWER_ANALYSIS' in result.stdout
    assert 'core:CONTROL_DECONTAM' in stages


def test_mock_preserves_existing_correction_and_network_choices():
    for name in ['examples/mock.local.yml', 'asv_pipeline_nextflow.yml']:
        config = yaml.safe_load((ROOT / name).read_text())
        assert config['optional']['batch_correction']['correction_policy'] == 'always'
        assert config['optional']['spieceasi']['pulsar_criterion'] == 'stars'
    validator = (ROOT / 'examples/mock_test/validate_results.py').read_text()
    assert 'within_truth >= 3' in validator


def test_no_tracked_environment_points_to_removed_shared_environment():
    for path in (ROOT / 'processes').rglob('env.yml'):
        assert path.exists(), f'Broken environment symlink: {path}'


def test_analysis_cohort_preserves_diversity_and_uses_staged_outlier_inputs():
    source = (ROOT / 'asv_pipeline.nf').read_text()
    optional = source.split('workflow optional {', 1)[1].split('workflow RUN_METADATA_ANALYSES {', 1)[0]
    assert 'asvFinalForDiversity = analysis_cohort_stage' not in optional
    assert 'metaMicroForDiversity = analysis_cohort_stage' not in optional
    for name in ['Indicspecies', 'VocCorrelation', 'Network', 'PowerAnalysis', 'MasterSummary']:
        assert f'asvFinalFor{name} = analysis_cohort_stage.counts' in optional
    assert optional.index('ANALYSIS_COHORT(') < optional.index('PLOT_UPSET(')
    outlier = source.split('process OUTLIER_CHECKER {', 1)[1].split('process COLLECTORS_CURVE {', 1)[0]
    assert r'\$PWD/${asv_clr}' in outlier
    assert r'\$PWD/${metadata_table}' in outlier


def test_staged_task_paths_are_resolved_by_the_task_shell():
    source = (ROOT / "asv_pipeline.nf").read_text()
    assert ".toAbsolutePath()" not in source
