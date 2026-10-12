"""Parallel execution must preserve scientific results and indexed resumption."""
import importlib
import json
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / 'processes/power_analysis_pipeline'))
import power_parallel as execution


def cases():
    rng = np.random.RandomState(81)
    patients = np.repeat(np.arange(12).astype(str), 2)
    types = np.tile(['A', 'B'], 12)
    groups = np.repeat(['Cancer'] * 6 + ['Control'] * 6, 2)
    counts = rng.poisson(12, (24, 5)).astype(float)
    counts[groups == 'Cancer', 0] *= 4
    names = list('abcde')
    for null in (False, True):
        yield 'power_permanova_stratified', 'run_power_simulation', (
            counts, patients, groups, names, [], 1., 4, 4), dict(n_perm=9, use_true_null=null)
        yield 'power_shannon_stratified', 'run_power_simulation', (
            counts, patients, groups, 4, 4), dict(use_true_null=null)
        yield 'power_taxonomic_abundance', 'run_power_simulation', (
            counts, patients, groups, names, {'taxa_indices': [0], 'fold_change': 2}, 4, 4), dict(use_true_null=null)
    yield 'power_sample_type_permanova', 'run_power_simulation_omnibus', (
        counts, patients, types, 6), dict(n_perm=9)
    yield 'power_sample_type_permanova', 'run_power_simulation_pairwise', (
        counts, patients, types, 'A', 'B', 6), dict(n_perm=9)
    yield 'power_sample_type_shannon', 'run_power_simulation', (
        counts, patients, types, ['A', 'B'], 6), {}


def assert_same(a, b):
    if isinstance(a, (list, tuple)):
        assert len(a) == len(b)
        for x, y in zip(a, b):
            assert_same(x, y)
    else:
        np.testing.assert_allclose(a, b, rtol=0, atol=0, equal_nan=True)


@pytest.mark.parametrize('module,name,inputs,options', list(cases()))
def test_simulations_parallel_resume(module, name, inputs, options, tmp_path):
    function = getattr(importlib.import_module(module), name)
    options.update(n_simulations=8, seed=12)
    execution.configure_parallel(1)
    serial = function(*inputs, **options)
    execution.configure_parallel(2, tmp_path)
    parallel = function(*inputs, **options)
    assert_same(serial, parallel)
    checkpoint, = tmp_path.glob('*.json')
    results = json.loads(checkpoint.read_text())
    checkpoint.write_text(json.dumps(results[:3]))
    execution.configure_parallel(3, tmp_path)
    assert_same(serial, function(*inputs, **options))
    assert len(json.loads(checkpoint.read_text())) == 8
    execution.configure_parallel(1)


def test_checkpoint_identity_and_completed_reuse(tmp_path):
    execution.configure_parallel(2, tmp_path)
    def replicate(i):
        return i * i
    assert list(execution.run_replicates(replicate, 6, {'seed': 1})) == [0, 1, 4, 9, 16, 25]
    # A completed checkpoint must not start workers or execute a replicate.
    def fail(i):
        raise AssertionError('cached replicate executed')
    fail.__qualname__ = replicate.__qualname__
    assert list(execution.run_replicates(fail, 6, {'seed': 1})) == [0, 1, 4, 9, 16, 25]
    list(execution.run_replicates(replicate, 6, {'seed': 2}))
    assert len(list(tmp_path.glob('*.json'))) == 2
    execution.configure_parallel(1)


def test_worker_failure_preserves_prefix(tmp_path):
    execution.configure_parallel(2, tmp_path)
    def replicate(i):
        if i == 3:
            raise RuntimeError('intentional failure')
        return i
    with pytest.raises(RuntimeError, match='intentional failure'):
        list(execution.run_replicates(replicate, 8, {}))
    checkpoint, = tmp_path.glob('*.json')
    assert json.loads(checkpoint.read_text()) == [0, 1, 2]
    execution.configure_parallel(1)


def test_r_isa_parallel_resume():
    import os
    import subprocess
    runner = os.environ.get('ASPIRE_TEST_RSCRIPT')
    if not runner:
        pytest.skip('Set ASPIRE_TEST_RSCRIPT to the power-analysis Rscript')
    root = Path(__file__).parents[1]
    result = subprocess.run([runner, str(root/'tests/power_parallel_checks.R'), str(root)],
                            capture_output=True, text=True, timeout=180)
    assert result.returncode == 0, result.stdout + result.stderr


def test_effect_bootstrap_parallel_resume(tmp_path):
    import estimate_effects
    a = np.arange(1, 12, dtype=float)
    b = np.arange(4, 15, dtype=float)
    execution.configure_parallel(1)
    expected = estimate_effects.bootstrap_cohens_d(a, b, n_bootstrap=30, seed=41)
    execution.configure_parallel(2, tmp_path)
    assert_same(expected, estimate_effects.bootstrap_cohens_d(a, b, n_bootstrap=30, seed=41))
    checkpoint, = tmp_path.glob('*.json')
    checkpoint.write_text(json.dumps(json.loads(checkpoint.read_text())[:5]))
    execution.configure_parallel(3, tmp_path)
    assert_same(expected, estimate_effects.bootstrap_cohens_d(a, b, n_bootstrap=30, seed=41))
    execution.configure_parallel(1)


def test_nextflow_worker_budget(tmp_path):
    import os
    import subprocess
    runner = os.environ.get('ASPIRE_TEST_NEXTFLOW')
    if not runner:
        pytest.skip('Set ASPIRE_TEST_NEXTFLOW to exercise actual worker configuration')
    root = Path(__file__).parents[1]
    source = (root/'asv_pipeline.nf').read_text()
    start = source.index('    def powerAnalysisWorkersRaw =')
    end = source.index('    if (powerAnalysisEnabled) log.info', start)
    fragment = source[start:end].replace('exit 1, ', 'throw new IllegalArgumentException(')
    # Convert the two exit statements into catchable errors for this isolated guard test.
    fragment = fragment.replace('"power_analysis.workers must be a positive integer or null"',
                                '"power_analysis.workers must be a positive integer or null")')
    fragment = fragment.replace('"Power-analysis CPU budget must be positive"',
                                '"Power-analysis CPU budget must be positive")')
    script = tmp_path/'workers.nf'
    script.write_text('nextflow.enable.dsl=2\n' + '''
def budget(value, threads, host) {
    def powerAnalysisConfig = [workers:value]
    def config = [resources:[threads:threads]]
    int pipelineThreads = host
''' + fragment + '''
    return powerAnalysisWorkers
}
def checkBudget() {
    assert budget(null, 32, 64) == 32
    assert budget(8, 32, 64) == 8
    assert budget(64, 32, 64) == 32
    assert budget(32, 32, 8) == 8
    assert budget(1, 32, 64) == 1
    for (bad in [0, -1, true, 'all', 1.5]) {
        boolean failed = false
        try { budget(bad, 32, 64) } catch (IllegalArgumentException e) { failed = true }
        assert failed
    }
}
workflow { checkBudget() }
''')
    result = subprocess.run([runner, 'run', str(script)], cwd=tmp_path,
                            env=dict(os.environ, NXF_OFFLINE='true', NXF_SYNTAX_PARSER='v1'),
                            capture_output=True, text=True, timeout=180)
    assert result.returncode == 0, result.stdout + result.stderr
