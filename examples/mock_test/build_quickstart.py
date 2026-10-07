#!/usr/bin/env python3
"""Build the small, deterministic ASPIRE installation fixture (standard library only)."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import random
from pathlib import Path

SEED = 20261007
FORWARD_PRIMER = 'GTGYCAGCMGCCGCGGTAA'.replace('Y', 'C').replace('M', 'A')
REVERSE_PRIMER = 'GGACTACNVGGGTWTCTAAT'.replace('N', 'A').replace('V', 'A').replace('W', 'A')
# The pipeline trims the 19/20 primer bases from these error-free reads.
FORWARD_PREFIX = FORWARD_PRIMER
REVERSE_PREFIX = REVERSE_PRIMER


def table(path, rows, columns=None):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns or list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def allocate(weights, total):
    values = [w / sum(weights) * total for w in weights]
    counts = [int(v) for v in values]
    for i in sorted(range(len(values)), key=lambda i: values[i] - counts[i], reverse=True)[:total-sum(counts)]:
        counts[i] += 1
    return counts


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def build(destination: Path):
    if destination.exists():
        raise SystemExit(f'Destination already exists; choose a new directory: {destination}')
    seed_dir = Path(__file__).resolve().parents[1] / 'quickstart'
    features = list(csv.DictReader((seed_dir / 'features.tsv').open(), delimiter='\t'))
    rng = random.Random(SEED)
    destination.mkdir(parents=True)
    (destination / 'fastq').mkdir()
    (destination / 'references').mkdir()
    metadata, chemistry, manifest, counts, controls = [], [], [], {}, []
    for patient in range(12):
        case = 'Control' if patient < 6 else 'Cancer'
        participant = f'P{patient+1:02d}'
        sites = [('Airways', 'Healthy' if case == 'Control' else 'TumorSide'), ('Oral', ''), ('Skin', '')]
        if case == 'Cancer':
            sites.insert(1, ('Airways', 'Contralateral'))
        for site, status in sites:
            sample = f'{participant}_{site}' + (f'_{status}' if status else '')
            row = dict(sample_id=sample, Sample=sample, Participant_ID=participant, Case=case,
                       lung_status=status, Type_Group=site, body_site=site,
                       study_role='biological_control' if site == 'Skin' else 'comparison',
                       synthetic_clinical=True, source_sample=participant,
                       batch=f'plate_{patient%2+1}', DNA_conc=10+patient,
                       is_negative_control=False, is_positive_control=False)
            metadata.append(row)
            weights = [0.] * len(features)
            if site == 'Skin':
                for i in range(24,28): weights[i] = rng.uniform(.7,1.3)
                total = 1800
            else:
                latent = [rng.gauss(0,1.2) for _ in range(3)]
                for i in range(24):
                    if i < 15:
                        weights[i] = math.exp(latent[i//5]+rng.gauss(0,.08))
                    elif i < 18:
                        weights[i] = (8 if site == 'Airways' else .2)*math.exp(rng.gauss(0,.2))
                    elif i < 21:
                        weights[i] = (8 if site == 'Oral' else .2)*math.exp(rng.gauss(0,.2))
                    else:
                        weights[i] = math.exp(rng.gauss(0,.6))*(2 if case == 'Cancer' else 1)
                # Sparse carry-over exposes both control-removal arms.
                if patient == 0:
                    for i in range(24,31): weights[i] = .12
                weights[31] = .4  # mitochondrial fixture in every biological sample
                total = 6000
            count = allocate(weights,total)
            counts[sample] = count
            if site != 'Skin':
                chemistry.append(dict(sample_id=sample,
                    acetone=round(1000*count[0]/total,8),
                    ethanol=round(1000*count[5]/total,8),
                    acetate=round(1000*count[10]/total,8)))
    for n, total in enumerate([800,1600,6000,8000],1):
        sample=f'Blank_{n}'
        metadata.append(dict(sample_id=sample,Sample=sample,Participant_ID=sample,Case='Control',
            lung_status='',Type_Group='Control',body_site='Control',study_role='technical_control',
            synthetic_clinical=True,source_sample=sample,batch=f'plate_{n%2+1}',DNA_conc=.02*n,
            is_negative_control=True,is_positive_control=False))
        counts[sample]=allocate([1. if 28<=i<=30 else 0. for i in range(len(features))],total)
        controls.append(dict(sample_id=sample,control_type='negative',total_reads=total,DNA_conc=.02*n))
    for row in metadata:
        sample=row['sample_id']
        paths=[f'fastq/{sample}_R{read}.fastq.gz' for read in (1,2)]
        handles=[gzip.GzipFile(filename=str(destination/p),mode='wb',mtime=0) for p in paths]
        try:
            serial=0
            for feature,count in zip(features,counts[sample]):
                sequence=feature['v4_sequence']
                reads=[(FORWARD_PREFIX+sequence)[:250],(REVERSE_PREFIX+reverse_complement(sequence))[:250]]
                for _ in range(count):
                    serial+=1
                    for read,handle in zip(reads,handles):
                        handle.write(f'@{sample}:{serial}\n{read}\n+\n{"I"*len(read)}\n'.encode())
        finally:
            for handle in handles: handle.close()
        manifest.append(dict(sample_id=sample,fastq_r1=paths[0],fastq_r2=paths[1]))
    table(destination/'sample_metadata.tsv',metadata)
    table(destination/'chemistry.tsv',chemistry)
    table(destination/'fastq_manifest.tsv',manifest)
    table(destination/'ground_truth_feature_registry.tsv',features)
    table(destination/'asv_counts.tsv',[dict(ASV_ID=f['ASV_ID'],**{s:c[i] for s,c in counts.items()}) for i,f in enumerate(features)])
    table(destination/'ground_truth_group_effects.tsv',[dict(ASV_ID=features[i]['ASV_ID'],group_column='Type_Group',group=site) for site,inds in [('Airways',range(15,18)),('Oral',range(18,21))] for i in inds])
    table(destination/'ground_truth_network_modules.tsv',[dict(ASV_ID=f['ASV_ID'],module=f['module']) for f in features if f['module']])
    table(destination/'ground_truth_asv_chem.tsv',[dict(ASV_ID=features[i]['ASV_ID'],compound=compound,direction='positive') for i,compound in [(0,'acetone'),(5,'ethanol'),(10,'acetate')]])
    table(destination/'ground_truth_extraction_controls.tsv',controls)
    table(destination/'ground_truth_reference_filters.tsv',[dict(ASV_ID=f['ASV_ID'],expected_filter='mitochondrial' if f['fixture_role']=='mitochondrial' else 'contaminant') for f in features[28:]])
    for name,items in [('contaminants',features[28:31]),('mitochondria',features[31:])]:
        (destination/'references'/f'{name}.fasta').write_text(''.join(f">{f['ASV_ID']}\n{f['v4_sequence']}\n" for f in items))
    (destination/'manifest.json').write_text(json.dumps(dict(name='aspire-quickstart',seed=SEED,
        patients=12,libraries=len(metadata),read_pairs=sum(map(sum,counts.values())),features=len(features),
        description='Deterministic installation fixture with strong planted signals and error-free V4 reads.',
        source=json.loads((seed_dir/'source.json').read_text())),indent=2)+'\n')
    table(destination/'checksums.tsv',[dict(relative_path=str(p.relative_to(destination)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(destination.rglob('*')) if p.is_file()])
    print(f'Created {destination}: {len(metadata)} libraries, {sum(map(sum,counts.values())):,} read pairs')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    build(args.output.expanduser().resolve())


if __name__=='__main__':
    main()
