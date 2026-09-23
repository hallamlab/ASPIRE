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
    assert 'baseAsvMeta = three_tier_stage.filtered_long' in body
    assert 'baseAsvFinal = three_tier_stage.filtered_wide' in body
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
    assert 'optional:THREE_TIER_DECONTAM' in stages


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
