#!/usr/bin/env python3
"""Run the standalone primer trimmer for one manifest sample, or check a cohort."""
import argparse
import csv
import importlib.util
import json
from pathlib import Path
import re
import shutil


def check_cohort(paths, output):
    rows = [json.loads(Path(p).read_text()) for p in paths]
    if not rows or len({r['sample_id'] for r in rows}) != len(rows):
        raise ValueError('Primer cohort requires unique, nonempty sample summaries')
    families = {r['family'] for r in rows}
    if len(families) != 1:
        raise ValueError('Different amplicon families across samples; split into separate ASPIRE runs: ' + ', '.join(sorted(families)))
    Path(output).write_text(json.dumps({'family': next(iter(families)), 'samples': rows}, indent=2) + '\n')


def run_sample(sample, r1, r2, settings, utility):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', sample):
        raise ValueError('Primer trimming sample IDs must contain only letters, digits, ._-')
    spec = importlib.util.spec_from_file_location('primer_utility', utility)
    trim = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(trim)
    config = json.loads(Path(settings).read_text())
    inputs = Path('primer_inputs'); inputs.mkdir()
    for mate, source in enumerate([r1, r2], 1):
        source = Path(source).resolve()
        suffix = '.fastq.gz' if source.name.endswith('.gz') else '.fastq'
        (inputs / f'{sample}_R{mate}{suffix}').symlink_to(source)
    catalogue = Path('primers.json')
    catalogue.write_text(json.dumps(config.pop('primers')))
    args = ['--fastq-dir', str(inputs), '--output-dir', 'audit', '--primers-json', str(catalogue)]
    for key, value in config.items():
        args += ['--' + key.replace('_', '-'), str(value)]
    status = trim.main(args)
    with Path('audit/summary.tsv').open() as handle:
        row = next(csv.DictReader(handle, delimiter='\t'))
    if status or row['status'] != 'trimmed':
        raise ValueError(f"{sample}: {row['status']}: {row['message']}; inspect audit/ reports. Split mixed-family libraries before ASV inference.")
    families = {s.split(':')[0] for s in row['selected_families'].split(';')}
    if len(families) != 1:
        raise ValueError(f'{sample}: mixed amplicon families')
    reads = Path('reads'); reads.mkdir()
    for category in ['trimmed', 'unmatched', 'too_short']:
        shutil.move(str(Path('audit') / category), reads / category)
    # Record portable public paths, separate from the utility's task-local manifests.
    summary = {'sample_id': sample, 'family': next(iter(families)),
               'selected_families': row['selected_families']}
    for key in ['screened_pairs', 'supported_screen_pairs', 'input_pairs', 'trimmed_pairs', 'unmatched_pairs', 'too_short_pairs']:
        summary[key] = int(row[key])
    for mate in (1, 2):
        summary[f'fastq_r{mate}'] = f'intermediates/primer_trimmed/{sample}/trimmed/{sample}_R{mate}.fastq.gz'
    Path(f'{sample}.primer_summary.json').write_text(json.dumps(summary, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='mode', required=True)
    run = sub.add_parser('run')
    for name in ['sample', 'r1', 'r2', 'settings', 'utility']:
        run.add_argument('--' + name, required=True)
    check = sub.add_parser('check')
    check.add_argument('--summaries', nargs='+', required=True)
    check.add_argument('--output', default='primer_cohort.json')
    a = parser.parse_args()
    if a.mode == 'check':
        check_cohort(a.summaries, a.output)
    else:
        run_sample(a.sample, a.r1, a.r2, a.settings, a.utility)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError) as exc:
        raise SystemExit(str(exc))
