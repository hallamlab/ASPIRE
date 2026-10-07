"""Verify the portable installation fixture against reads, truth and full-run settings."""
import gzip
import importlib.util
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


def module(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'examples/mock_test' / f'{name}.py')
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


@pytest.fixture(scope='module')
def dataset(tmp_path_factory):
    path = tmp_path_factory.mktemp('quickstart') / 'dataset'
    module('build_quickstart').build(path)
    return path


def test_reads_reconstruct_truth_and_have_correct_primer_trims(dataset):
    builder = module('build_quickstart')
    features = pd.read_csv(dataset / 'ground_truth_feature_registry.tsv', sep='\t')
    sequence_ids = dict(zip(features.v4_sequence, features.ASV_ID))
    expected = pd.read_csv(dataset / 'asv_counts.tsv', sep='\t', index_col=0)
    manifest = pd.read_csv(dataset / 'fastq_manifest.tsv', sep='\t')
    total = 0
    for row in manifest.itertuples(index=False):
        observed = {key: 0 for key in expected.index}
        with gzip.open(dataset / row.fastq_r1, 'rt') as r1, gzip.open(dataset / row.fastq_r2, 'rt') as r2:
            while name := r1.readline():
                assert r2.readline() == name
                forward = r1.readline().strip()[19:]
                reverse = builder.reverse_complement(r2.readline().strip()[20:])
                matches = [seq for seq in sequence_ids if seq.startswith(forward) and seq.endswith(reverse)]
                assert len(matches) == 1
                observed[sequence_ids[matches[0]]] += 1
                for handle in (r1, r2):
                    assert handle.readline().strip() == '+'
                    assert set(handle.readline().strip()) == {'I'}
            assert r2.readline() == ''
        assert observed == expected[row.sample_id].to_dict()
        total += sum(observed.values())
    assert total == 218000


def test_reproducible_checksums_and_full_configuration(dataset, tmp_path):
    builder = module('build_quickstart')
    rebuilt = tmp_path / 'rebuilt'
    builder.build(rebuilt)
    assert (dataset / 'checksums.tsv').read_bytes() == (rebuilt / 'checksums.tsv').read_bytes()
    configure = module('configure_mock_run')
    configure.validate_checksums(dataset)
    metadata = pd.read_csv(dataset / 'sample_metadata.tsv', sep='\t')
    assert metadata.Type_Group.value_counts().to_dict() == {'Airways': 18, 'Oral': 12, 'Skin': 12, 'Control': 4}
    configure.validate_tabular_inputs(dataset, set(metadata.sample_id))
    config = configure.build_config(ROOT/'examples/mock.local.yml', dataset, tmp_path/'results', tmp_path/'runtime', ROOT, 4)
    section = lambda key: configure.config_section(config,key)
    assert section('control_decontam')['technical_enabled']
    assert section('control_decontam')['bio_control_enabled']
    assert section('control_decontam')['min_biological_reads'] == 5000
    assert section('filter_counts')['min_prevalence_fraction'] == .05
    assert section('filter_counts')['min_relative_abundance_pct'] == .1
    assert section('indicspecies')['perms'] == 999
    assert section('spieceasi')['rep_num'] == 50
    assert section('power_analysis')['n_simulations'] == 10
    assert section('power_analysis')['n_perm'] == 49
    assert section('voc_correlation')['patient_inference']
    assert all(section(key)['enabled'] for key in ('diversity','power_analysis','lung_status_analysis','taxonomy_patient_aware','network_topology'))
    with pytest.raises(SystemExit, match='already exists'):
        builder.build(dataset)


def test_configurator_preserves_unchanged_resume_inputs(dataset, tmp_path):
    config = tmp_path / 'run.yml'
    command = [sys.executable, str(ROOT / 'examples/mock_test/configure_mock_run.py'),
               '--dataset', str(dataset), '--output', str(tmp_path / 'results'),
               '--config-out', str(config), '--threads', '4']
    subprocess.run(command, check=True, capture_output=True, text=True)
    generated = [config, config.with_suffix('.metadata.tsv'), config.with_suffix('.manifest.tsv')]
    before = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in generated}
    subprocess.run(command, check=True, capture_output=True, text=True)
    assert before == {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in generated}
