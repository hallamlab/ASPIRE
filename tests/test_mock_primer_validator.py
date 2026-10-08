import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('mock_validator', ROOT/'examples/mock_test/validate_results.py')
v = importlib.util.module_from_spec(spec); spec.loader.exec_module(v)


def test_primer_validator_detects_broken_handoff_and_missing_audit(tmp_path):
    def write(name, value):
        path=tmp_path/name; path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(value))
    write('summary/tables/run_config.yml', {'core':{'primer_trimming':{'enabled':True},'fastp':dict.fromkeys(['trim_front_r1','trim_tail_r1','trim_front_r2','trim_tail_r2'],0)}})
    row=dict(sample_id='S1',family='515F_806R',input_pairs=10,trimmed_pairs=8,unmatched_pairs=2,too_short_pairs=0,screened_pairs=10,supported_screen_pairs=8,fastq_r1='reads/r1.gz',fastq_r2='reads/r2.gz')
    write('modules/primer_trimming/tables/primer_cohort.json',{'family':'515F_806R','samples':[row]})
    write('modules/primer_trimming/tables/S1.primer_summary.json',row)
    write('intermediates/fastp/S1.fastp.json',{'summary':{'before_filtering':{'total_reads':16}}})
    for mate in (1,2): write(f'reads/r{mate}.gz',{})
    audit=v.Audit(); v.validate_primer_trimming(tmp_path,{'S1'},audit); assert not audit.failures
    write('intermediates/fastp/S1.fastp.json',{'summary':{'before_filtering':{'total_reads':20}}})
    audit=v.Audit(); v.validate_primer_trimming(tmp_path,{'S1'},audit); assert len(audit.failures)==1
    (tmp_path/'modules/primer_trimming/tables/primer_cohort.json').unlink()
    audit=v.Audit(); v.validate_primer_trimming(tmp_path,{'S1'},audit); assert len(audit.failures)==1
