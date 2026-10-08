#!/usr/bin/env python3
"""Generate the current TECH/BIO workflow SVG/PDF figures from one source.

Requires CairoSVG from docs/diagram-requirements.txt for vector PDF export.
"""
from pathlib import Path
from xml.etree import ElementTree as ET
import cairosvg

ROOT = Path(__file__).resolve().parents[1]
NS = 'http://www.w3.org/2000/svg'
ET.register_namespace('', NS)

FULL = [
    ('Optional primer removal', 'PRIMER_TRIM (Cutadapt) → PRIMER_TRIM_CHECK',
     'Detect paired primers; audit discarded pairs; require one amplicon family.'),
    ('Read quality control', 'FASTP_QC → MERGE_READS → FILTER_READS → RELABEL_FILTERED',
     'All biological samples and controls share read processing.'),
    ('ASV construction', 'CONCAT_FASTAS → DEREPLICATE → DENOISE → CHIMERA_CHECK → CREATE_COUNT_MATRIX',
     'Preserve original ASV counts and sequences, including control-only ASVs.'),
    ('Alignment and full taxonomy', 'SINA_TRIM → TAXONOMY',
     'Assign taxonomy to all inferred ASVs before biological feature filtering.'),
    ('CONTROL_DECONTAM · optional, metadata-defined', 'Biological samples: ≥5,000 post-QC reads; nonzero controls: no depth cutoff',
     'TECH: biological + blanks. BIO: same biological cohort + sample controls.',
     'Independent prevalence scores < configured cutoffs → union ASV removal.',
     'Retained counts unchanged. Positive controls: QC only. Unused arms can be disabled.'),
    ('Reference screening', 'PREPARE_BLAST_DATABASES → MITOMASTER → MITO_DECONTAM',
     'Local host / mitochondrial evidence; external MITOMASTER lookup is optional.'),
    ('FILTER_ASVS · combined biological feature filtering', 'Initial abundance / nonzero → group and non-target → abundance / taxonomy → final nonzero',
     'Mock: ≥0.1% RA in any biological sample AND nonzero counts in ≥5% of biological samples.',
     'No second depth cutoff after CONTROL_DECONTAM. Final counts and FASTA agree.'),
    ('PLOT_METADATA · final biological tables', 'ASV_target.tsv → synchronized long / wide tables and metadata',
     'CAMI: Airways and Oral only; Skin BIO controls and TECH blanks excluded.'),
    ('Optional analysis-table preparation', 'Group diagnostics / validated labels → batch correction and count selection',
     'Diversity retains all final biological samples; configured study-group exclusions apply to other analyses.'),
    ('ANALYSIS_COHORT · optional study-group selection', 'Synchronize metadata, counts, long ASV tables and CLR sample rows; publish a selection audit',
     'Diagnostic cohort selection precedes group diagnostics; downstream selection follows batch preparation.'),
    ('Optional analysis branches · dependencies govern execution', 'Diversity / ordination · indicators · patient / paired contrasts · group-effect power',
     'ASV–VOC associations · SPIEC-EASI → modules / topology → annotated networks',
     'Final filtered FASTA + genome references → ASV–MAG links and network overlays.'),
    ('Publish results and provenance', 'Module tables / plots → master summaries → HTML report and checksum inventories',
     'GENERAL_STATS runs as a parallel QC branch; audit intermediates remain available.'),
]
BRIEF = [
    ('Reads → ASVs → full taxonomy', 'Shared QC for biological samples, BIO controls and TECH blanks',
     'Optional Cutadapt → fastp / VSEARCH → ASV count matrix → SINA / QIIME 2 taxonomy'),
    ('Optional TECH / BIO control decontamination', '≥5,000 post-QC reads for biological samples; retain all nonzero controls',
     'Independent prevalence tests → union ASV removal; retained counts unchanged',
     'Positive controls stay in QC. Either control arm may be disabled.'),
    ('Reference screening → FILTER_ASVS', 'Host / mitochondrial evidence → biological abundance and taxonomy filters',
     'Mock: ≥0.1% RA in any biological sample AND nonzero counts in ≥5% of biological samples',
     'This percentage is separate from decontam’s prevalence-score cutoffs.'),
    ('Metadata → selected analysis tables', 'PLOT_METADATA receives the final filtered biological count table',
     'Optional cohort selection for other analyses; diversity retains all final biological samples'),
    ('Optional analyses and integration', 'Diversity · indicators · patient contrasts · power · participant-level ASV–VOC tests',
     'SPIEC-EASI networks / topology / modules · filtered ASV-to-genome linkage'),
    ('Results, QC and provenance', 'Tables · editable figures · removal audits · HTML report · logs and checksums',
     'CAMI cohort: Airways / Oral biological; Skin BIO controls; Control TECH blanks'),
]


def element(parent, tag, text=None, **attrs):
    child = ET.SubElement(parent, f'{{{NS}}}{tag}', {k.replace('_', '-'): str(v) for k, v in attrs.items()})
    child.text = text
    return child


def build(name, rows):
    width = 1400
    height = 155 + sum(80 + 29 * (len(row) - 1) for row in rows) + 55
    svg = ET.Element(f'{{{NS}}}svg', dict(width=str(width), height=str(height), viewBox=f'0 0 {width} {height}', role='img', **{'aria-labelledby': 'title desc'}))
    element(svg, 'title', 'ASPIRE: TECH/BIO control decontamination workflow', id='title')
    element(svg, 'desc', 'Full taxonomy precedes control prevalence filtering, reference screening, combined ASV filtering and metadata tables. Optional analyses follow the selected biological tables.', id='desc')
    element(svg, 'rect', width=width, height=height, fill='white')
    element(svg, 'style', 'text{font-family:Arial,Helvetica,sans-serif;fill:#172b3a}.heading{font-weight:bold}')
    defs = element(svg, 'defs')
    marker = element(defs, 'marker', id='arrow', viewBox='0 0 10 10', refX=9, refY=5, markerWidth=7, markerHeight=7, orient='auto')
    element(marker, 'path', d='M 0 0 L 10 5 L 0 10 z', fill='#526d82')
    element(svg, 'text', 'ASPIRE · TECH/BIO workflow', x=width//2, y=48, font_size=32, text_anchor='middle', **{'class': 'heading'})
    element(svg, 'text', 'YAML configuration · Nextflow / Mamba · resumable execution · audited outputs', x=width//2, y=83, font_size=22, text_anchor='middle')
    y=115
    for i, row in enumerate(rows):
        h=57 + 29 * (len(row)-1)
        color = '#e7f3ed' if ('control decontamination' in row[0] or 'CONTROL_DECONTAM' in row[0]) else '#edf3f8'
        element(svg, 'rect', x=35, y=y, width=width-70, height=h, rx=10, fill=color, stroke='#70889a')
        element(svg, 'text', row[0], x=58, y=y+33, font_size=24, **{'class':'heading'})
        for j, line in enumerate(row[1:]):
            element(svg, 'text', line, x=58, y=y+65+j*29, font_size=20)
        y+=h
        if i<len(rows)-1:
            element(svg, 'path', d=f'M {width//2} {y+3} V {y+21}', stroke='#526d82', stroke_width=2, marker_end='url(#arrow)')
        y+=23
    element(svg, 'text', 'Arrows show the main data path; enabled analyses may run concurrently. Control decontamination can be bypassed.', x=width//2, y=y+19, font_size=19, text_anchor='middle')
    target=ROOT/'docs/assets'/f'{name}.svg'
    ET.ElementTree(svg).write(target, encoding='unicode')
    cairosvg.svg2pdf(url=str(target), write_to=str(target.with_suffix('.pdf')))
    print(f'Built {target.relative_to(ROOT)} and PDF')


if __name__ == '__main__':
    build('workflow',FULL)
    build('workflow-brief',BRIEF)
