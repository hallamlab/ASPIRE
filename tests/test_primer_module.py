"""Exercise integrated primer trimming, family gates and real Nextflow wiring."""
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / 'processes/primer_trimming/primer_trimming.py'
UTILITY = ROOT / 'scripts/trim_amplicon_primers.py'
spec = importlib.util.spec_from_file_location('primer_module', RUNNER)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
F = 'GTGCCAGCAGCCGCGGTAA'
R = 'GGACTACCGGGGTTTCTAAT'
R2 = 'CCGTCAATTYMTTTRAGTTT'.replace('Y', 'C').replace('M', 'A').replace('R', 'A')
INSERT = 'TACGATCGATGCTAGCTAGCATCGATCGAT' * 3


def fixture(tmp_path, monkeypatch, mixed=False):
    cutadapt = shutil.which(os.environ.get('CUTADAPT', 'cutadapt'))
    if not cutadapt:
        pytest.skip('Cutadapt required')
    monkeypatch.setenv('PATH', str(Path(cutadapt).parent) + os.pathsep + os.environ['PATH'])
    pairs = [(F + INSERT, R + INSERT)] * 10
    if mixed:
        pairs += [(F + INSERT, R2 + INSERT)] * 10
    raw = tmp_path / 'raw'; raw.mkdir()
    for mate in (1, 2):
        (raw / f'original_R{mate}.fastq').write_text(''.join(
            f'@read{i}/{mate}\n{pair[mate-1]}\n+\n' + 'I' * len(pair[mate-1]) + '\n'
            for i, pair in enumerate(pairs)))
    settings = tmp_path / 'settings.json'
    settings.write_text(json.dumps({'primers': [
        {'name': 'V4', 'forward': F, 'reverse': R},
        {'name': 'V4V5', 'forward': F, 'reverse': R2}], 'threads': 1}))
    return raw, settings


def test_manifest_id_counts_and_cohort(tmp_path, monkeypatch):
    raw, settings = fixture(tmp_path, monkeypatch)
    monkeypatch.chdir(tmp_path)
    m.run_sample('manifest-01_A', raw/'original_R1.fastq', raw/'original_R2.fastq', settings, UTILITY)
    summary = tmp_path/'manifest-01_A.primer_summary.json'
    row = json.loads(summary.read_text())
    assert row['sample_id'] == 'manifest-01_A'
    assert row['family'] == 'V4'
    assert row['input_pairs'] == row['trimmed_pairs'] == 10
    assert row['unmatched_pairs'] == row['too_short_pairs'] == 0
    assert (tmp_path/'reads/trimmed/manifest-01_A_R1.fastq.gz').exists()
    assert (tmp_path/'audit/primer_detection.tsv').exists()
    m.check_cohort([summary], tmp_path/'cohort.json')
    second = tmp_path/'second.json'
    second.write_text(json.dumps(dict(row, sample_id='other', family='V4V5')))
    with pytest.raises(ValueError, match='Different amplicon families'):
        m.check_cohort([summary, second], tmp_path/'bad.json')
    with pytest.raises(ValueError, match='unique'):
        m.check_cohort([summary, summary], tmp_path/'bad.json')


def test_mixed_sample_stops(tmp_path, monkeypatch):
    raw, settings = fixture(tmp_path, monkeypatch, mixed=True)
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match='mixed'):
        m.run_sample('mixed', raw/'original_R1.fastq', raw/'original_R2.fastq', settings, UTILITY)
    assert not (tmp_path/'mixed.primer_summary.json').exists()


def test_nextflow_primer_to_fastp(tmp_path, monkeypatch):
    nf = os.environ.get('ASPIRE_TEST_NEXTFLOW')
    if not nf or not shutil.which('fastp'):
        pytest.skip('Set ASPIRE_TEST_NEXTFLOW and put fastp on PATH')
    raw, settings = fixture(tmp_path, monkeypatch)
    source = (ROOT/'asv_pipeline.nf').read_text()
    processes = source.split('process PRIMER_TRIM {', 1)[1].split('process MERGE_READS {', 1)[0]
    # Exercise the actual dependency gate and process definitions.
    start = source.index('    def rawReadsForAsv = raw_reads_input')
    end = source.index('FASTP_QC(rawReadsForAsv)', start) + len('FASTP_QC(rawReadsForAsv)')
    graph = source[start:end]
    script = tmp_path/'smoke.nf'
    script.write_text('nextflow.enable.dsl=2\n' + f'''
def sampleThreads = 1
def primerEnvPath = '{ROOT}/processes/primer_trimming/env.yml'
def condaEnvPath = '{ROOT}/env.yml'
def outputDir = '{tmp_path}/published'
def dirMap = [fastp: "${{outputDir}}/fastp"]
def fastpTrimValues = [front_r1:0,tail_r1:0,front_r2:0,tail_r2:0]
def primerTrimmingEnabled = true
def primerRunner = '{RUNNER}'
def primerUtility = '{UTILITY}'
def primerSettingsJson = new File('{settings}').text
workflow {{
    raw_reads_input = Channel.of(tuple([sample_id:'manifest-01_A',paired:true],file('{raw}/original_R1.fastq'),file('{raw}/original_R2.fastq')))
''' + graph + '\n}\nprocess PRIMER_TRIM {' + processes)
    result = subprocess.run([nf, 'run', str(script)], cwd=tmp_path,
        env=dict(os.environ, NXF_OFFLINE='true', NXF_SYNTAX_PARSER='v1'),
        capture_output=True, text=True, timeout=180)
    assert result.returncode == 0, result.stdout + result.stderr
    output = tmp_path/'published'
    cohort = json.loads((output/'primer_trimming/primer_cohort.json').read_text())
    assert cohort['family'] == 'V4'
    assert (output/'primer_trimming/manifest-01_A/summary.tsv').exists()
    assert (output/'primer_trimmed/manifest-01_A/trimmed/manifest-01_A_R1.fastq.gz').exists()
    fastp = json.loads((output/'fastp/manifest-01_A.fastp.json').read_text())
    assert fastp['summary']['before_filtering']['total_reads'] == 20
    spec = importlib.util.spec_from_file_location('layout', ROOT/'processes/output_layout/organize_outputs.py')
    layout = importlib.util.module_from_spec(spec); spec.loader.exec_module(layout)
    layout.organize(output)
    summary = json.loads((output/'modules/primer_trimming/tables/manifest-01_A.primer_summary.json').read_text())
    assert (output / summary['fastq_r1']).exists()
    assert (output/'modules/primer_trimming/tables/manifest-01_A/summary.tsv').exists()
    assert (output/'modules/primer_trimming/tables/primer_cohort.json').exists()



def test_screen_tolerance_requires_complete_primers(tmp_path):
    spec = importlib.util.spec_from_file_location('trim_tolerant', UTILITY)
    trim = importlib.util.module_from_spec(spec); spec.loader.exec_module(trim)
    import re
    pattern = re.compile(''.join('[' + trim.IUPAC[b] + ']' for b in F))
    one = F[:5] + 'T' + F[6:]
    two = one[:6] + 'T' + one[7:]
    assert trim.primer_position(one + INSERT, F, pattern, 12, 0) is None
    assert trim.primer_position('AA' + one + INSERT, F, pattern, 12, .1) == 2
    assert trim.primer_position(two + INSERT, F, pattern, 0, .1) is None
    assert trim.primer_position(F[:-1], F, pattern, 12, .1) is None


def test_config_guardrails(tmp_path):
    nf = os.environ.get('ASPIRE_TEST_NEXTFLOW')
    if not nf:
        pytest.skip('Set ASPIRE_TEST_NEXTFLOW')
    source = (ROOT/'asv_pipeline.nf').read_text()
    boundary = source[source.index('def primerConfig = '):source.index('\ndef merge', source.index('def primerConfig = '))]
    cases = [({}, {'front_r1':0,'tail_r1':0,'front_r2':0,'tail_r2':0}, False, 0),
             ({'enabled':True}, {'front_r1':19}, False, 1),
             ({'enabled':True}, {'front_r1':0}, True, 1),
             ({'screen_error_rate':.3}, {'front_r1':0}, False, 1),
             ({'min_pair_fraction':0}, {'front_r1':0}, False, 1),
             ({'primers':[{'name':'bad','forward':'AAAA','reverse':R}]}, {'front_r1':0}, False, 1),
             ({'enabled':'false'}, {'front_r1':0}, False, 1)]
    for i,(primer, fastp, single, code) in enumerate(cases):
        config = tmp_path/f'case{i}.json'; config.write_text(json.dumps({'primer_trimming':primer}))
        script = tmp_path/f'case{i}.nf'
        script.write_text('nextflow.enable.dsl=2\n' + f'''
def config = new groovy.json.JsonSlurper().parse(new File('{config}'))
def resolvePath = {{ value -> value }}
def sampleThreads = 1
def fastpTrimValues = new groovy.json.JsonSlurper().parseText('{json.dumps(fastp)}')
def allowSingleEnd = {str(single).lower()}
''' + boundary + '\nworkflow {}\n')
        result = subprocess.run([nf,'run',str(script),'-preview'],cwd=tmp_path,
          env=dict(os.environ,NXF_OFFLINE='true',NXF_SYNTAX_PARSER='v1'),capture_output=True,text=True,timeout=180)
        assert result.returncode == code, result.stdout + result.stderr
        if code:
            assert 'primer_trimming' in result.stdout + result.stderr
