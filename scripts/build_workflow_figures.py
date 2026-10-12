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

# title, input, (process label, tool label) pairs, output
FULL = [('Primer removal\n& read QC',
  'Paired reads\nand manifest',
  [('Primer trimming\n(optional)', 'Cutadapt'), ('Read QC', 'fastp'), ('Merge / filter', 'VSEARCH')],
  'Processed\nreads'),
 ('ASV\nconstruction',
  'Processed\nreads',
  [('Dereplicate', 'VSEARCH'),
   ('Denoise', 'VSEARCH'),
   ('Chimera check\nand counts', 'VSEARCH')],
  'Raw ASVs\nand counts'),
 ('Alignment\n& taxonomy',
  'All inferred ASVs\nand references',
  [('Align sequences', 'SINA'), ('Assign taxonomy', 'QIIME 2 / VSEARCH')],
  'Full ASV\ntaxonomy'),
 ('Control\ndecontamination',
  'Raw counts and\nsample classes',
  [('Biological\ndepth QC', 'ASPIRE / pandas'),
   ('TECH / BIO\nprevalence tests', 'decontam'),
   ('Union of\nflagged ASVs', 'ASPIRE / pandas')],
  'Decontaminated\nbiological counts'),
 ('Non-target\nscreening',
  'ASVs and\nreference sequences',
  [('Mitochondrial lookup\n(optional API)', 'MITOMASTER'), ('Local reference\nscreen', 'BLAST+')],
  'Non-target\nevidence'),
 ('Biological\nfeature filtering',
  'Biological counts\nand taxonomy',
  [('Abundance and\nprevalence', 'ASPIRE / pandas'),
   ('Host / unassigned\ntaxon exclusions', 'ASPIRE / pandas')],
  'Final microbial\ncounts and FASTA'),
 ('Metadata\n& preparation',
  'Final counts\nand metadata',
  [('Metadata-linked\nlong / wide tables', 'pandas'),
   ('Group diagnostics\nselected cohort', 'scikit-learn'),
   ('Batch correction\n(optional)', 'ConQuR')],
  'Synchronized\nanalysis tables'),
 ('Diversity\n& cohort selection',
  'Prepared\nbiological tables',
  [('Full cohort\ndiversity', 'scikit-bio / vegan'),
   ('Select study groups\nfor other analyses', 'ASPIRE / pandas')],
  'Diversity and\nselected tables'),
 ('Associations\n& study design',
  'Selected tables\nand optional VOCs',
  [('Indicators and\npatient contrasts', 'indicspecies / vegan'),
   ('ASV–VOC\nassociations', 'SciPy / statsmodels'),
   ('Patient-count\npower simulations', 'NumPy / SciPy')],
  'Association tables\nand figures'),
 ('Networks\n& genome links',
  'Selected tables\nand optional genomes',
  [('Network\ninference', 'SPIEC-EASI'),
   ('Modules and\nnull topology', 'igraph / NetworkX'),
   ('ASV–MAG links\nand overlays', 'BLAST+ / Biopython')],
  'Networks and\nlinkage tables'),
 ('Reports\n& provenance',
  'Module outputs\nand QC audits',
  [('Master summaries', 'ASPIRE / pandas'),
   ('All-library Sankey\nand HTML report', 'ASPIRE'),
   ('Logs and\nchecksum inventory', 'Nextflow / ASPIRE')],
  'Published results\nand run record')]
BRIEF = [('Reads\n& ASVs',
  'Reads, metadata\nand references',
  [('Primer trimming\nand read QC', 'Cutadapt / fastp'),
   ('ASV construction', 'VSEARCH'),
   ('Full taxonomy', 'SINA / QIIME 2')],
  'Raw ASVs\nand taxonomy'),
 ('Control\ndecontamination',
  'Raw counts and\nsample classes',
  [('Biological-only\ndepth QC', 'ASPIRE / pandas'),
   ('TECH / BIO\nprevalence', 'decontam'),
   ('Union ASV\nremoval', 'ASPIRE / pandas')],
  'Biological\ncounts'),
 ('Microbial\nfeature selection',
  'Counts, taxonomy\nand references',
  [('Non-target\nscreening', 'MITOMASTER / BLAST+'),
   ('Abundance /\nprevalence filters', 'ASPIRE / pandas'),
   ('Taxonomy\nexclusions', 'ASPIRE / pandas')],
  'Final microbial\ncounts and FASTA'),
 ('Metadata\n& cohorts',
  'Final counts\nand metadata',
  [('Linked metadata\nand preparation', 'pandas / ConQuR'),
   ('Full-cohort\ndiversity', 'scikit-bio / vegan'),
   ('Study groups for\nother analyses', 'ASPIRE / pandas')],
  'Diversity and\nselected tables'),
 ('Analyses\n& integration',
  'Selected tables;\nVOCs / genomes',
  [('Indicators /\npatient tests', 'indicspecies / vegan'),
   ('Associations\nand power', 'SciPy / statsmodels'),
   ('Networks /\ngenome links', 'SPIEC-EASI / BLAST+')],
  'Analysis tables\nand figures'),
 ('Results\n& provenance',
  'Module outputs\nand QC audits',
  [('All-library\nSankey', 'ASPIRE'),
   ('Integrated\nHTML report', 'ASPIRE'),
   ('Logs and\nchecksums', 'Nextflow / ASPIRE')],
  'Auditable\nrun outputs')]

def element(parent, tag, text=None, **attrs):
    child = ET.SubElement(parent, f'{{{NS}}}{tag}', {k.replace('_', '-'): str(v) for k, v in attrs.items()})
    child.text = text
    return child


def build(name, rows):
    width, gap = 1500, 175
    height = 335 + gap * len(rows) + (110 if rows is FULL else 0)
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

    for i, (title, source, steps, output) in enumerate(rows):
        y = 345 + i * gap + (110 if rows is FULL and i > 7 else 0)
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
            for yy, (label, tool) in zip([y-25, y+110], steps):
                node(900, yy, 'compute')
                text(900, yy-62, label, 21)
                text(900, yy+37, tool, 19)
                element(svg, 'path', d=f'M443.5 {y} H550 V{yy} H885', **{'class': 'wire'})
                element(svg, 'path', d=f'M915 {yy} H1260 V{y} H1356.5', **{'class': 'wire'})
        else:
            last = 443.5
            for x, (label, tool) in zip(xs, steps):
                node(x, y, 'compute')
                text(x, y-66, label, 21)
                text(x, y+37, tool, 19)
                wire(last, y, x-15, y)
                last = x+15
            wire(last, y, 1356.5, y)
        if i < len(rows)-1:
            wire(290, y+14.5, 290, y+gap+(110 if rows is FULL and i == 7 else 0)-14.5)
    text(750, height-20, 'Numbered modules group related operations; arrows summarize flow and enabled branches follow their dependencies.', 19)
    target = ROOT / 'docs/assets' / f'{name}.svg'
    ET.ElementTree(svg).write(target, encoding='unicode')
    cairosvg.svg2pdf(url=str(target), write_to=str(target.with_suffix('.pdf')))
    print(f'Built {target.relative_to(ROOT)} and PDF')


if __name__ == '__main__':
    build('workflow', FULL)
    build('workflow-brief', BRIEF)
