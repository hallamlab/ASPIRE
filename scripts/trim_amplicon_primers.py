#!/usr/bin/env python3
"""Standalone paired FASTQ primer-family screening and Cutadapt trimming.

No ASPIRE workflow/configuration is modified. Python standard library plus
Cutadapt >=4 is required. See scripts/PRIMER_TRIMMING.md.
"""
import argparse
import csv
import gzip
import itertools
import json
import logging
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

IUPAC = dict(zip('ACGTRYSWKMBDHVN', ['A','C','G','T','AG','CT','CG','AT','GT','AC','CGT','AGT','ACT','ACG','ACGT']))
COMPLEMENT = str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN')
DEFAULTS = [
    {'name': '515F_806R', 'forward': 'GTGYCAGCMGCCGCGGTAA', 'reverse': 'GGACTACNVGGGTWTCTAAT'},
    {'name': '515F_926R', 'forward': 'GTGYCAGCMGCCGCGGTAA', 'reverse': 'CCGYCAATTYMTTTRAGTTT'},
]
LOG = logging.getLogger('primer_trim')


def discover(directory, recursive=False, excluded=None):
    pairs = {}
    unsupported = []
    for p in sorted(directory.rglob('*') if recursive else directory.iterdir()):
        if excluded and (p == excluded or excluded in p.parents):
            continue
        if not p.is_file() or not re.search(r'\.(?:fastq|fq)(?:\.gz)?$', p.name, re.I):
            continue
        stem = re.sub(r'\.(?:fastq|fq)(?:\.gz)?$', '', p.name, flags=re.I)
        m = re.fullmatch(r'(.+?)[._]R([12])(?:[._](\d+))?', stem, re.I)
        if not m:
            m = re.fullmatch(r'(.+?)[._]([12])', stem)
        if not m:
            unsupported.append(str(p))
            continue
        sample = m[1] + (('_' + m[3]) if len(m.groups()) == 3 and m[3] else '')
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', sample):
            raise ValueError(f'Unsupported sample filename: {p.name}; use letters/numbers/._-')
        mate = int(m[2])
        if mate in pairs.setdefault(sample, {}):
            raise ValueError(f'Duplicate sample/mate {sample}/R{mate}; rename lanes or directories explicitly')
        pairs[sample][mate] = p.resolve()
    if unsupported:
        raise ValueError('Cannot assign FASTQ mate from filename: ' + ', '.join(unsupported[:8]))
    missing = [s for s, mates in pairs.items() if set(mates) != {1, 2}]
    if missing:
        raise ValueError('Missing paired FASTQ for: ' + ', '.join(missing))
    if not pairs:
        raise ValueError('No paired FASTQ files found')
    return pairs


def fastq(path):
    opener = gzip.open if path.name.lower().endswith('.gz') else open
    with opener(path, 'rt') as f:
        while True:
            header = f.readline()
            if not header:
                return
            seq, plus, quality = f.readline(), f.readline(), f.readline()
            seq, quality = seq.rstrip('\r\n'), quality.rstrip('\r\n')
            if not header.startswith('@') or not plus.startswith('+') or len(seq) != len(quality) or not quality:
                raise ValueError(f'Malformed four-line FASTQ record in {path}: {header.strip()}')
            yield header.rstrip('\r\n'), seq.upper(), plus.rstrip('\r\n'), quality


def read_id(header):
    return re.sub(r'/[12]$', '', header.split()[0])


def paired_records(r1, r2):
    for a, b in itertools.zip_longest(fastq(r1), fastq(r2)):
        if a is None or b is None:
            raise ValueError('FASTQ mates contain different numbers of records')
        if read_id(a[0]) != read_id(b[0]):
            raise ValueError(f'FASTQ mates out of sync: {a[0]} versus {b[0]}')
        yield a, b


def load_primers(path):
    families = json.loads(path.read_text()) if path else DEFAULTS
    if not isinstance(families, list) or not families:
        raise ValueError('Primer JSON must be a nonempty list of name/forward/reverse objects')
    names = set()
    for f in families:
        if not isinstance(f, dict) or set(f) != {'name', 'forward', 'reverse'}:
            raise ValueError('Each primer entry must contain exactly name, forward and reverse')
        if not isinstance(f['name'], str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', f['name']) or f['name'] in names:
            raise ValueError('Primer names must be unique and contain only letters, digits, _ or -')
        names.add(f['name'])
        for side in ['forward', 'reverse']:
            if not isinstance(f[side], str):
                raise ValueError(f'Primer sequence must be text: {f["name"]}/{side}')
            seq = f[side].upper()
            if not 8 <= len(seq) <= 100 or set(seq) - set(IUPAC):
                raise ValueError(f'Invalid IUPAC primer: {f["name"]}/{side}')
            f[side] = seq
    return families


def screen(r1, r2, families, sample_reads, max_prefix):
    patterns = {s: re.compile(''.join('['+IUPAC[b]+']' for b in s))
                for f in families for s in [f['forward'], f['reverse']]}
    counts, offsets = Counter(), {}
    n = supported = ambiguous = 0
    for reads in itertools.islice(paired_records(r1, r2), sample_reads):
        n += 1
        hits = []
        for i, f in enumerate(families):
            for orientation in ['forward', 'swapped']:
                primers = [f['forward'], f['reverse']]
                if orientation == 'swapped':
                    primers.reverse()
                positions = []
                for read, seq in zip(reads, primers):
                    match = patterns[seq].search(read[1][:max_prefix + len(seq)])
                    positions.append(match.start() if match else None)
                if all(pos is not None and pos <= max_prefix for pos in positions):
                    key = (i, orientation)
                    hits.append(key)
                    counts[key] += 1
                    offsets.setdefault(key, [Counter(), Counter()])
                    for j, pos in enumerate(positions):
                        offsets[key][j][pos] += 1
        supported += bool(hits)
        ambiguous += len(hits) > 1
    return n, supported, ambiguous, counts, offsets


def adapter_specs(selected, families, mate, max_prefix, end_overlap):
    specs = []
    for index, orientation in selected:
        f = families[index]
        primers = [f['forward'], f['reverse']]
        if orientation == 'swapped':
            primers.reverse()
        front = primers[mate-1]
        opposite = primers[2-mate].translate(COMPLEMENT)[::-1]
        for offset in range(max_prefix+1):
            name = f'{f["name"]}_{orientation}_p{offset}'
            # Anchored full 5' primer required; 3' read-through primer optional.
            specs.append(f'{name}=^'+('N'*offset)+front+'...'+opposite+f';optional;min_overlap={min(end_overlap,len(opposite))}')
    return specs


def write_tsv(path, rows, fields):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--fastq-dir', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True, help='New or empty output directory')
    p.add_argument('--recursive', action='store_true')
    p.add_argument('--primers-json', type=Path, help='Optional candidate primer-pair catalogue')
    p.add_argument('--cutadapt', default='cutadapt', help='Cutadapt executable (>=4)')
    p.add_argument('--threads', type=int, default=4, help='Cutadapt threads per sample; samples run sequentially')
    p.add_argument('--sample-reads', type=int, default=5000, help='First N read pairs screened per sample')
    p.add_argument('--max-prefix', type=int, default=12, help='Maximum bases allowed before the 5-prime primer')
    p.add_argument('--min-pair-fraction', type=float, default=.5, help='Minimum screened fraction with a supported primer pair')
    p.add_argument('--min-family-fraction', type=float, default=.01, help='Minimum screened fraction for a family/orientation to be used')
    p.add_argument('--min-family-reads', type=int, default=5)
    p.add_argument('--error-rate', type=float, default=.1, help='Cutadapt mismatch rate; screening uses exact IUPAC matches')
    p.add_argument('--end-overlap', type=int, default=12, help='Minimum matching bases for optional 3-prime read-through trimming')
    p.add_argument('--minimum-length', type=int, default=1, help='Shorter processed pairs go to too_short/')
    p.add_argument('--log-file', type=Path, help='Additional progress log (default OUTPUT/run.log)')
    args = p.parse_args(argv)
    if any(x < 1 for x in [args.threads,args.sample_reads,args.min_family_reads,args.end_overlap,args.minimum_length]):
        p.error('Threads, read counts, overlap and minimum length must be positive')
    if not 0 <= args.max_prefix <= 30 or not 0 <= args.error_rate <= .25:
        p.error('max-prefix must be 0–30 and error-rate 0–0.25')
    if not 0 < args.min_pair_fraction <= 1 or not 0 < args.min_family_fraction <= 1:
        p.error('Fractions must be in (0,1]')
    return args


def main(argv=None):
    args = parse_args(argv)
    args.fastq_dir = args.fastq_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    out = args.output_dir
    if out == args.fastq_dir:
        raise ValueError('Input and output directories must differ')
    if out.exists() and any(out.iterdir()):
        raise ValueError(f'Output is not empty; choose a new directory: {out}')
    families = load_primers(args.primers_json)
    pairs = discover(args.fastq_dir,args.recursive,out)
    exe = shutil.which(args.cutadapt)
    if not exe:
        raise ValueError('Cutadapt not found. Install: mamba create -n aspire_primer_qc -c conda-forge -c bioconda python=3.11 cutadapt=5.2')
    version = subprocess.check_output([exe,'--version'],text=True).strip()
    if int(version.split('.')[0]) < 4:
        raise ValueError('Cutadapt >=4 is required')
    out.mkdir(parents=True,exist_ok=True)
    for sub in ['trimmed','unmatched','too_short','reports']:
        (out/sub).mkdir()
    log = args.log_file or out/'run.log'
    log.parent.mkdir(parents=True,exist_ok=True)
    logging.basicConfig(level=logging.INFO,format='%(asctime)s %(message)s',handlers=[logging.StreamHandler(sys.stdout),logging.FileHandler(log)],force=True)
    config = {k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()}
    config.update(cutadapt_version=version,primer_families=families,screening='exact IUPAC matches in the first N pairs; family counts may overlap',variant_caveat='Family-level compatibility does not identify exact degenerate oligo formulation')
    (out/'run_config.json').write_text(json.dumps(config,indent=2)+'\n')
    summaries, detection, inputs = [], [], []
    for number,(sample,mates) in enumerate(pairs.items(),1):
        LOG.info('[%d/%d] %s: screening',number,len(pairs),sample)
        for mate,path in mates.items():
            stat=path.stat();inputs.append(dict(sample_id=sample,mate=mate,path=str(path),bytes=stat.st_size,mtime_ns=stat.st_mtime_ns))
        row = dict(sample_id=sample,status='failed',screened_pairs=0,supported_screen_pairs=0,ambiguous_screen_pairs=0,selected_families='',input_pairs='',trimmed_pairs='',unmatched_pairs='',too_short_pairs='',fastq_r1='',fastq_r2='',message='')
        try:
            n, supported, ambiguous, counts, offsets = screen(mates[1],mates[2],families,args.sample_reads,args.max_prefix)
            selected = [key for key,c in sorted(counts.items()) if c >= args.min_family_reads and c/max(n,1) >= args.min_family_fraction]
            row.update(screened_pairs=n,supported_screen_pairs=supported,ambiguous_screen_pairs=ambiguous,selected_families=';'.join(families[i]['name']+':'+o for i,o in selected))
            for i,f in enumerate(families):
                for orientation in ['forward','swapped']:
                    key=(i,orientation);ofs=offsets.get(key,[{},{}])
                    detection.append(dict(sample_id=sample,family=f['name'],orientation=orientation,screened_pairs=n,matching_pairs=counts[key],fraction=counts[key]/max(n,1),selected=key in selected,R1_offsets=json.dumps(ofs[0]),R2_offsets=json.dumps(ofs[1])))
            if not n or supported/n < args.min_pair_fraction or not selected:
                row.update(status='unsupported',message='Insufficient candidate-primer support; original inputs unchanged; no trimmed FASTQs produced')
                LOG.warning('%s: %s (%d/%d screened pairs)',sample,row['message'],supported,n)
            else:
                with tempfile.TemporaryDirectory(prefix='.'+sample+'-',dir=out) as tmp:
                    tmp=Path(tmp)
                    command=[exe,'-j',str(args.threads),'--no-indels','-e',str(args.error_rate),'--pair-filter=any','-m',str(args.minimum_length)]
                    for mate,flag in [(1,'-a'),(2,'-A')]:
                        for spec in adapter_specs(selected,families,mate,args.max_prefix,args.end_overlap):command += [flag,spec]
                    targets={}
                    for category,flags in [('trimmed',('-o','-p')),('unmatched',('--untrimmed-output','--untrimmed-paired-output')),('too_short',('--too-short-output','--too-short-paired-output'))]:
                        for mate,flag in enumerate(flags,1):
                            name=f'{sample}_R{mate}.fastq.gz';temp=tmp/(category+'_'+name)
                            targets[temp]=out/category/name;command += [flag,str(temp)]
                    report=out/'reports'/f'{sample}.cutadapt.json'
                    command += ['--json',str(report),str(mates[1]),str(mates[2])]
                    (out/'reports'/f'{sample}.command.json').write_text(json.dumps(command,indent=2)+'\n')
                    LOG.info('%s: Cutadapt using %s',sample,row['selected_families'])
                    with (out/'reports'/f'{sample}.cutadapt.log').open('w') as f:
                        subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=True)
                    stats=json.loads(report.read_text())['read_counts'];filters=stats['filtered']
                    total=stats['input'];kept=stats['output'];unknown=filters.get('discard_untrimmed') or 0;short=filters.get('too_short') or 0
                    if kept+unknown+short != total:
                        raise ValueError('Cutadapt read-pair accounting failed; inspect report')
                    for source,target in targets.items():source.replace(target)
                    mixed=len({i for i,_ in selected})>1
                    row.update(status='trimmed_mixed_families' if mixed else 'trimmed',input_pairs=total,trimmed_pairs=kept,unmatched_pairs=unknown,too_short_pairs=short,fastq_r1=str(out/'trimmed'/f'{sample}_R1.fastq.gz'),fastq_r2=str(out/'trimmed'/f'{sample}_R2.fastq.gz'),message='Mixed amplicon families: do not pool unlike regions for ASV inference' if mixed else '')
                    if kept==0:
                        row.update(status='no_pairs_retained',message='No pairs passed trimming; inspect reports')
                    LOG.info('%s: %d/%d pairs retained; %d unmatched; %d too short',sample,kept,total,unknown,short)
        except (ValueError,OSError,subprocess.SubprocessError,KeyError) as exc:
            row.update(status='failed',message=str(exc));LOG.error('%s: %s',sample,exc)
        summaries.append(row)
        write_tsv(out/'summary.tsv',summaries,list(row))
        if detection:write_tsv(out/'primer_detection.tsv',detection,list(detection[0]))
        write_tsv(out/'input_manifest.tsv',inputs,list(inputs[0]))
    retained=[{k:r[k] for k in ['sample_id','fastq_r1','fastq_r2']} for r in summaries if r['status'].startswith('trimmed')]
    write_tsv(out/'trimmed_manifest.tsv',retained,['sample_id','fastq_r1','fastq_r2'])
    problems=sum(not r['status'].startswith('trimmed') for r in summaries)
    LOG.info('Complete: %d samples; %d require review. Results: %s',len(summaries),problems,out)
    return 2 if problems else 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError,OSError,subprocess.SubprocessError) as exc:
        print(f'error: {exc}',file=sys.stderr)
        sys.exit(2)
