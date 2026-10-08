#!/usr/bin/env python3
"""Build ASPIRE's MP-style nodal SVG/PDF workflows.

The visual contract is docs/WORKFLOW_STYLE.md. Keep scientific content in ROWS;
keep the numbered module spine, node grammar and palette in the renderer.
Requires CairoSVG from docs/diagram-requirements.txt.
"""
from pathlib import Path
from xml.etree import ElementTree as ET
import cairosvg

ROOT = Path(__file__).resolve().parents[1]
NS = 'http://www.w3.org/2000/svg'
ET.register_namespace('', NS)
COLORS = {'data': ('#F5F5F5', '#666666'), 'input': ('#DAE8FC', '#6C8EBF'),
          'output': ('#D5E8D4', '#82B366')}

# title, input, compute nodes, output, explanatory annotation
FULL = [
    ('Primer removal\n& read QC', 'Paired reads\nand manifest',
     ['Primer trimming\nCutadapt · optional', 'Read QC\nfastp', 'Merge / filter\nVSEARCH'], 'Processed\nreads',
     'Paired-primer audit and amplicon-family check precede shared QC for biological samples and controls.'),
    ('ASV\nconstruction', 'Processed\nreads',
     ['Dereplicate', 'Denoise', 'Chimera check\nand counts'], 'Raw ASVs\nand counts',
     'Preserve original counts and control-only ASVs for the control prevalence tests.'),
    ('Alignment\n& taxonomy', 'All inferred ASVs\nand references',
     ['Align sequences\nSINA', 'Assign taxonomy\nQIIME 2'], 'Full ASV\ntaxonomy',
     'Taxonomy covers the raw ASV set before biological feature filtering.'),
    ('Control\ndecontamination', 'Raw counts and\nsample classes',
     ['Biological\ndepth QC', 'TECH / BIO\nprevalence tests', 'Union of\nflagged ASVs'], 'Decontaminated\nbiological counts',
     'Optional · Biological depth cutoff is configurable; nonzero TECH / BIO controls retain their original counts.'),
    ('Non-target\nscreening', 'ASVs and\nreference sequences',
     ['MITOMASTER\noptional API', 'Local BLAST\nreference screen'], 'Non-target\nevidence',
     'Configured API retry and failure policies publish availability audits alongside local reference results.'),
    ('Biological\nfeature filtering', 'Biological counts\nand taxonomy',
     ['Abundance and\nprevalence', 'Host / unassigned\ntaxon exclusions'], 'Final microbial\ncounts and FASTA',
     'Mock: ≥0.1% relative abundance in any biological sample AND nonzero counts in ≥5% of biological samples.'),
    ('Metadata\n& preparation', 'Final counts\nand metadata',
     ['Metadata-linked\nlong / wide tables', 'Group diagnostics\nselected cohort', 'Batch preparation\nfull cohort'], 'Synchronized\nanalysis tables',
     'Metadata plots retain all final biological samples; batch correction is optional.'),
    ('Diversity\n& cohort selection', 'Prepared\nbiological tables',
     ['Full cohort\ndiversity', 'Select study groups\nfor other analyses'], 'Diversity and\nselected tables',
     'Parallel uses of prepared tables: excluded groups remain in diversity; sample-selection audits accompany other analyses.'),
    ('Associations\n& study design', 'Selected tables\nand optional VOCs',
     ['Indicators and\npatient contrasts', 'ASV–VOC\nassociations', 'Patient-count\npower simulations'], 'Association tables\nand figures',
     'Optional branches share synchronized inputs; VOC clustermaps cluster VOC rows and ASV columns.'),
    ('Networks\n& genome links', 'Selected tables\nand optional genomes',
     ['SPIEC-EASI\ninference', 'Modules and\nnull topology', 'ASV–MAG links\nand overlays'], 'Networks and\nlinkage tables',
     'Genome linkage uses final filtered ASV sequences and supplied references; enabled branches follow their dependencies.'),
    ('Reports\n& provenance', 'Module outputs\nand QC audits',
     ['Master summaries', 'All-library Sankey\nand HTML report', 'Logs and\nchecksum inventory'], 'Published results\nand run record',
     'The Sankey includes TECH / BIO controls and explicit dropout bands; kept nodes remain above removed nodes.'),
]
BRIEF = [
    ('Reads\n& ASVs', 'Reads, metadata\nand references',
     ['Optional Cutadapt\nthen shared QC', 'ASV construction', 'Full taxonomy'], 'Raw ASVs\nand taxonomy',
     'Biological samples and controls share primer processing and read QC.'),
    ('Control\ndecontamination', 'Raw counts and\nsample classes',
     ['Biological-only\ndepth QC', 'TECH / BIO\nprevalence', 'Union ASV\nremoval'], 'Biological\ncounts',
     'Optional · Each control arm is configurable; nonzero controls are exempt from the biological depth cutoff.'),
    ('Microbial\nfeature selection', 'Counts, taxonomy\nand references',
     ['Non-target\nscreening', 'Abundance /\nprevalence filters', 'Taxonomy\nexclusions'], 'Final microbial\ncounts and FASTA',
     'The mock uses 0.1% abundance in at least one biological sample and 5% nonzero biological prevalence.'),
    ('Metadata\n& cohorts', 'Final counts\nand metadata',
     ['Linked metadata\nand preparation', 'Full-cohort\ndiversity', 'Study groups for\nother analyses'], 'Diversity and\nselected tables',
     'Configured group exclusions apply to other analyses; metadata plots and diversity retain those groups.'),
    ('Analyses\n& integration', 'Selected tables;\nVOCs / genomes',
     ['Indicators /\npatient tests', 'Associations\nand power', 'Networks /\ngenome links'], 'Analysis tables\nand figures',
     'Optional analyses run according to their input dependencies and available study data.'),
    ('Results\n& provenance', 'Module outputs\nand QC audits',
     ['All-library\nSankey', 'Integrated\nHTML report', 'Logs and\nchecksums'], 'Auditable\nrun outputs',
     'Explicit control-dropout bands accompany the retained-sample and retained-read accounting.'),
]


def element(parent, tag, text=None, **attrs):
    child = ET.SubElement(parent, f'{{{NS}}}{tag}', {k.replace('_', '-'): str(v) for k, v in attrs.items()})
    child.text = text
    return child


def build(name, rows):
    width, gap = 1500, 175
    height = 335 + gap * len(rows) + (60 if rows is FULL else 0)
    svg = ET.Element(f'{{{NS}}}svg', dict(width=str(width), height=str(height), viewBox=f'0 0 {width} {height}', role='img', **{'aria-labelledby': 'title desc', 'data-workflow-style': 'mp-nodal-v1'}))
    element(svg, 'title', 'ASPIRE: amplicon analysis workflow', id='title')
    element(svg, 'desc', 'Numbered conceptual modules use the MetaPathways nodal style. Shared read processing and full taxonomy precede control decontamination and final feature filtering. Diversity keeps the full biological cohort; selected groups feed other analyses. Module numbers are not a serial execution schedule.', id='desc')
    element(svg, 'rect', width=width, height=height, fill='#FFFFFF')
    element(svg, 'style', 'text{font-family:"Times New Roman",Times,serif;fill:#111111}.wire{fill:none;stroke:#111111;stroke-width:1.5;marker-end:url(#arrow)}')
    defs = element(svg, 'defs')
    marker = element(defs, 'marker', id='arrow', viewBox='0 0 10 10', refX=10, refY=5, markerWidth=7, markerHeight=7, orient='auto')
    element(marker, 'path', d='M0 0 L10 5 L0 10z', fill='#111111')

    def text(x, y, value, size=22, anchor='middle'):
        parent = element(svg, 'text', x=x, y=y, font_size=size, text_anchor=anchor)
        for i, line in enumerate(value.split('\n')):
            element(parent, 'tspan', line, x=x, dy=0 if i == 0 else size * 1.15)

    def panel(x, y, w, h, fill, stroke):
        element(svg, 'rect', x=x, y=y, width=w, height=h, rx=16, fill=fill, stroke=stroke, stroke_width=1.5)

    def node(x, y, kind):
        fill, stroke = COLORS.get(kind, COLORS['data'])
        if kind == 'module':
            element(svg, 'rect', x=x-13, y=y-13, width=26, height=26, fill='#F5F5F5', stroke='#111111', stroke_width=3, **{'data-node': kind})
        elif kind == 'compute':
            element(svg, 'path', d=f'M{x} {y-13} L{x+13} {y} L{x} {y+13} L{x-13} {y}Z', fill=fill, stroke=stroke, stroke_width=3, **{'data-node': kind})
        else:
            element(svg, 'circle', cx=x, cy=y, r=12, fill=fill, stroke=stroke, stroke_width=3, **{'data-node': kind})

    def wire(x1, y1, x2, y2, arrow=True):
        attrs = {'class': 'wire'} if arrow else dict(fill='none', stroke='#111111', stroke_width=1.5)
        element(svg, 'path', d=f'M{x1} {y1} L{x2} {y2}', **attrs)

    panel(20, 20, 1460, 80, '#CCCCCC', '#666666')
    text(750, 54, 'ASPIRE • Amplicon analysis with Nextflow', 29)
    text(750, 83, 'Mamba environments • Resource controls • Resumable execution • Audited outputs', 21)
    panel(20, 120, 770, 130, *COLORS['input'])
    text(42, 152, 'Inputs', 25, 'start')
    text(42, 181, 'Paired amplicon reads, sample manifest, metadata and references\nConfigured biological / TECH / BIO sample classes\nOptional: VOC measurements and genome references', 21, 'start')
    panel(820, 120, 660, 100, '#CCCCCC', '#666666')
    for x, kind in zip([880, 1010, 1140, 1270, 1400], ['module', 'compute', 'data', 'input', 'output']):
        text(x, 150, kind.title(), 21)
        node(x, 183, kind)

    for i, (title, source, steps, output, note) in enumerate(rows):
        y = 345 + i * gap + (60 if rows is FULL and i > 7 else 0)
        text(140, y-8, title, 27)
        node(290, y, 'module')
        text(290, y+7, str(i+1), 20)
        node(430, y, 'input')
        text(430, y-66, source, 21)
        node(1370, y, 'output')
        text(1370, y-66, output, 21)
        wire(304.5, y, 416.5, y)
        xs = [650, 920, 1160] if len(steps) == 3 else [730, 1070]
        # Diversity and selected-cohort analysis are parallel consumers.
        if title == 'Diversity\n& cohort selection':
            for yy, label in zip([y-25, y+35], steps):
                node(900, yy, 'compute')
                text(900, yy-48 if yy < y else yy+30, label, 21)
                element(svg, 'path', d=f'M443.5 {y} H550 V{yy} H885', **{'class': 'wire'})
                element(svg, 'path', d=f'M915 {yy} H1260 V{y} H1356.5', **{'class': 'wire'})
            text(900, y+107, 'Full cohort and selected cohort are exported separately; sample-selection audits record exclusions.', 19)
        else:
            last = 443.5
            for x, label in zip(xs, steps):
                node(x, y, 'compute')
                text(x, y-66, label, 21)
                wire(last, y, x-15, y)
                last = x+15
            wire(last, y, 1356.5, y)
            text(900, y+49, note, 18)
        if i < len(rows)-1:
            wire(290, y+14.5, 290, y+gap+(60 if rows is FULL and i == 7 else 0)-14.5)
    text(750, height-20, 'Numbered modules group related operations; arrows summarize flow and enabled branches follow their dependencies.', 19)
    target = ROOT / 'docs/assets' / f'{name}.svg'
    ET.ElementTree(svg).write(target, encoding='unicode')
    cairosvg.svg2pdf(url=str(target), write_to=str(target.with_suffix('.pdf')))
    print(f'Built {target.relative_to(ROOT)} and PDF')


if __name__ == '__main__':
    build('workflow', FULL)
    build('workflow-brief', BRIEF)
