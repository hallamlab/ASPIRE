"""Exercise the real Groovy configuration boundary without scheduling analyses."""
import os
from pathlib import Path
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_parameter_aliases_and_units(tmp_path):
    runner = os.environ.get('ASPIRE_TEST_NEXTFLOW')
    if not runner:
        pytest.skip('Set ASPIRE_TEST_NEXTFLOW to exercise Nextflow configuration normalization')
    source = (ROOT / 'asv_pipeline.nf').read_text()
    boundary = source.split('def flattenTieredConfig = ', 1)[1].split('\n\ndef config\n', 1)[0]
    script = tmp_path / 'aliases.nf'
    script.write_text('nextflow.enable.dsl=2\ndef normalize = ' + boundary + '''
assert normalize([core:[table_filter:[min_sample_sum:5000, min_asv_sum:0]], standard:[filter_counts:[abundance_threshold:0.1]], optional:[three_tier_decontam:[technical_threshold:0.1, bio_control_threshold:0.2], spieceasi:[min_prevalence:0.05, min_rel_abund:0.001]]]) == normalize([core:[table_filter:[min_sample_reads:5000, min_relative_abundance_pct:0], control_decontam:[technical_score_threshold:0.1, bio_control_score_threshold:0.2]], standard:[filter_counts:[min_relative_abundance_pct:0.1]], optional:[spieceasi:[min_prevalence_fraction:0.05, min_relative_abundance_fraction:0.001]]])
assert normalize([voc_correlation:[isa_brush_groups:['Airways'], isa_exclude_all_types_from_brush:false]]) == normalize([voc_correlation:[isa_focus_groups:['Airways'], isa_exclude_all_types_from_focus:false]])
assert normalize([spieceasi:[min_prevalence:5]]).spieceasi.min_prevalence_fraction == 0.05
assert normalize([spieceasi:[min_prevalence_fraction:0]]).spieceasi.min_prevalence_fraction == 0
for (bad in [
    [filter_counts:[abundance_threshold:0.1, min_relative_abundance_pct:0.1]],
    [control_decontam:[technical_threshold:0.1, technical_score_threshold:0.1]],
    [voc_correlation:[isa_brush_groups:['Airways'], isa_focus_groups:['Airways']]],
    [spieceasi:[min_prevalence_fraction:5]],
    [spieceasi:[min_relative_abundance_fraction:-0.1]],
    [filter_counts:[min_relative_abundance_pct:101]],
    [table_filter:[min_relative_abundance_pct:Double.NaN]],
    [control_decontam:[:], three_tier_decontam:[:]]
]) {
    def rejected = false
    try { normalize(bad) } catch (IllegalArgumentException expected) { rejected = true }
    assert rejected: "Invalid or ambiguous configuration was accepted: ${bad}"
}
println 'Parameter aliases and explicit units verified'
workflow { }
''')
    result = subprocess.run([runner, 'run', str(script), '-preview'], cwd=tmp_path,
                            env=dict(os.environ, NXF_OFFLINE='true', NXF_SYNTAX_PARSER='v1'),
                            capture_output=True, text=True, timeout=180)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'Parameter aliases and explicit units verified' in result.stdout
