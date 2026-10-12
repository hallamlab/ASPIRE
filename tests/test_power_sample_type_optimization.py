"""Regression tests against the unchanged serial simulation implementation."""
import sys
from pathlib import Path
from itertools import combinations

import numpy as np
import pytest
from scipy.stats import wilcoxon

sys.path.insert(0, str(Path(__file__).parents[1] / 'processes' / 'power_analysis_pipeline'))
import power_taxonomic_sample_type as power


def fixture_data():
    rng = np.random.RandomState(123)
    patients = np.repeat(np.arange(8), 4)
    types = np.tile(['A', 'B', 'B', 'C'], 8)
    counts = rng.randint(0, 10, (32, 5)).astype(float)
    counts[rng.random_sample(counts.shape) < .35] = 0
    # Incomplete pairing and repeated samples of the same patient/type.
    keep = np.arange(32) % 7 != 0
    return counts[keep], patients[keep], types[keep], list('abcde')


def test_exact_probability_matches_scipy():
    rng = np.random.RandomState(17)
    for n in [3, 5, 10, 13, 14, 20]:
        vectors = [np.zeros(n), np.ones(n), np.arange(n, dtype=float),
                   rng.normal(size=n)]
        vectors += [rng.randint(-3, 4, n).astype(float) for _ in range(8)]
        for d in vectors:
            expected = wilcoxon(d, zero_method='wilcox', correction=False,
                                alternative='two-sided', method='auto').pvalue
            np.testing.assert_allclose(power.fast_wilcoxon(d, np.zeros(n)), expected,
                                       rtol=0, atol=1e-14, equal_nan=True)


def test_precomputed_profiles_preserve_individual_pvalues():
    counts, patients, types, names = fixture_data()
    unique, profiles = power._prepare_profiles(counts, patients, types)
    for seed in range(5):
        boot = power.bootstrap_patients_sample_types(counts, patients, types, 10, seed)
        abundance, stypes = power.patient_level_abundance_by_type(*boot)
        expected = []
        for taxon in range(len(names)):
            pairs = []
            for a, b in combinations(stypes, 2):
                common = sorted(set(abundance[a]) & set(abundance[b]))
                if len(common) >= 3:
                    x = [abundance[a][p][taxon] for p in common]
                    y = [abundance[b][p][taxon] for p in common]
                    pairs.append(wilcoxon(x, y, method='auto').pvalue)
            if pairs:
                expected.append(min(pairs))
        np.testing.assert_allclose(power._simulation_pvalues(unique, profiles, 10, seed),
                                   expected, rtol=0, atol=1e-14, equal_nan=True)


def test_parallel_and_checkpoint_equivalence(tmp_path, monkeypatch):
    data = fixture_data()
    options = dict(n_patients=10, n_simulations=8, seed=9, alpha=.2)
    expected = power.run_power_simulation_reference(*data, **options)
    assert power.run_power_simulation(*data, **options, workers=1) == expected
    assert power.run_power_simulation(*data, **options, workers=2,
                                      checkpoint_dir=tmp_path) == expected
    checkpoint = next(tmp_path.glob('*.json'))
    import json
    completed = json.loads(checkpoint.read_text())
    checkpoint.write_text(json.dumps(completed[:3]))
    assert power.run_power_simulation(*data, **options, workers=2,
                                      checkpoint_dir=tmp_path) == expected
    def should_not_run(index):
        raise AssertionError('Completed simulations must be reused')
    monkeypatch.setattr(power, '_run_replicate', should_not_run)
    assert power.run_power_simulation(*data, **options, checkpoint_dir=tmp_path) == expected


def test_checkpoint_changes_with_input_and_seed(tmp_path):
    data = fixture_data()
    for seed in [1, 2]:
        power.run_power_simulation(*data, n_patients=3, n_simulations=1,
                                   seed=seed, checkpoint_dir=tmp_path)
    assert len(list(tmp_path.glob('*.json'))) == 2
    data[0][0, 0] += 1
    power.run_power_simulation(*data, n_patients=3, n_simulations=1,
                               seed=1, checkpoint_dir=tmp_path)
    assert len(list(tmp_path.glob('*.json'))) == 3


def test_invalid_simulation_counts():
    with pytest.raises(ValueError):
        power.run_power_simulation(*fixture_data(), n_patients=3, n_simulations=0)


def test_command_line_outputs_and_resume(tmp_path):
    import subprocess
    import pandas as pd
    counts, patients, types, taxa = fixture_data()
    rows = [dict(sample=f's{i}', Participant_ID=patients[i], type_group=types[i],
                 Phylum=taxon, Family=taxon, count=counts[i, j])
            for i in range(len(counts)) for j, taxon in enumerate(taxa)]
    source = tmp_path / 'long.tsv'
    pd.DataFrame(rows).to_csv(source, sep='\t', index=False)
    command = [sys.executable, power.__file__, '--data-long', str(source),
               '--outdir', str(tmp_path), '--n-simulations', '2', '--sample-sizes', '3',
               '--workers', '2']
    subprocess.run(command, check=True, capture_output=True, text=True)
    output = tmp_path / 'taxonomic_sample_type_power.tsv'
    before = output.read_bytes()
    resumed = subprocess.run(command, check=True, capture_output=True, text=True)
    assert 'Resuming 2/2' in resumed.stdout
    assert output.read_bytes() == before
    results = pd.read_csv(output, sep='\t')
    assert results.tax_level.tolist() == ['Phylum', 'Family']
    assert (results.n_simulations == 2).all()


def test_permanova_cached_gram_preserves_seeded_permutations():
    import power_permanova as permanova
    rng = np.random.RandomState(72)
    matrix = permanova.bray_curtis_from_counts(rng.randint(1, 20, (12, 8)))
    patients = np.repeat(np.arange(6), 2)
    groups = np.repeat(['A', 'A', 'A', 'B', 'B', 'B'], 2)
    def original_r2(labels):
        n = len(labels)
        h = np.eye(n) - np.ones((n, n)) / n
        gram = -.5 * h @ (matrix ** 2) @ h
        between = 0
        for group in np.unique(labels):
            mask = labels == group
            between += np.sum(gram[np.ix_(mask, mask)]) / mask.sum()
        return between / np.trace(gram)
    observed = original_r2(groups)
    rng = np.random.RandomState(42)
    perms = [original_r2(np.repeat(rng.permutation(groups[::2]), 2)) for _ in range(199)]
    expected = (observed, (1 + np.sum(np.array(perms) >= observed)) / 200)
    assert permanova.permanova_permutation_test(matrix, groups, patients) == expected
