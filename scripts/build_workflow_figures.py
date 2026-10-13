#!/usr/bin/env python3
"""Build ASPIRE's MP-style nodal SVG/PDF workflows.

The visual contract is docs/WORKFLOW_STYLE.md. Keep scientific content in FULL / BRIEF;
keep the numbered module spine, node grammar and palette in the renderer.
Requires CairoSVG from docs/diagram-requirements.txt.
"""
from pathlib import Path
import json
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
  [('Biological\ndepth QC', ''),
   ('TECH / BIO\nprevalence tests', 'decontam'),
   ('Union of\nflagged ASVs', '')],
  'Decontaminated\nbiological counts'),
 ('Non-target\nscreening',
  'ASVs and\nreference sequences',
  [('Mitochondrial lookup\n(optional API)', 'MITOMASTER'), ('Local reference\nscreen', 'BLAST+')],
  'Non-target\nevidence'),
 ('Biological\nfeature filtering',
  'Biological counts\nand taxonomy',
  [('Abundance and\nprevalence', ''),
   ('Host / unassigned\ntaxon exclusions', '')],
  'Final microbial\ncounts and FASTA'),
 ('Metadata\n& preparation',
  'Final counts\nand metadata',
  [('Metadata-linked\nlong / wide tables', ''),
   ('Group diagnostics\nselected cohort', 'scikit-learn'),
   ('Batch correction\n(optional)', 'ConQuR / UMAP / HDBSCAN')],
  'Synchronized\nanalysis tables'),
 ('Diversity\n& cohort selection',
  'Prepared\nbiological tables',
  [('Full cohort\ndiversity', 'scikit-bio / vegan'),
   ('Select study groups\nfor other analyses', '')],
  'Diversity and\nselected tables'),
 ('Associations\n& study design',
  'Selected tables\nand optional VOCs',
  [('Indicators and\npatient contrasts', 'indicspecies / vegan'),
   ('ASV–VOC\nassociations', ''),
   ('Patient-count\npower simulations', '')],
  'Association tables\nand figures'),
 ('Networks\n& genome links',
  'Selected tables\nand optional genomes',
  [('Network\ninference', 'SPIEC-EASI'),
   ('Modules and\nnull topology', 'igraph / NetworkX'),
   ('ASV–MAG links\nand overlays', 'BLAST+')],
  'Networks and\nlinkage tables'),
 ('Reports\n& provenance',
  'Module outputs\nand QC audits',
  [('Read statistics /\nmaster summaries', 'SeqKit'),
   ('All-library Sankey\nand HTML report', ''),
   ('Logs and\nchecksum inventory', 'Nextflow')],
  'Published results\nand run record')]
BRIEF = [('FASTQ → ASVs',
  'Paired FASTQs',
  [('Primer removal /\nread QC', 'Cutadapt / fastp'),
   ('Merge / filter', 'VSEARCH'),
   ('Denoise / remove\nchimeras / count', 'VSEARCH')],
  'Initial ASVs\n+ counts'),
 ('ASV filtering',
  'Initial ASVs\n+ counts',
  [('Taxonomy', 'SINA / QIIME 2'),
   ('Depth QC /\ncontrol decontamination', 'decontam (TECH / BIO)'),
   ('Non-target removal /\nabundance / prevalence', 'MITOMASTER / BLAST+')],
  'Filtered ASVs\n+ counts'),
 ('Metadata &\ntransformations',
  'ASVs + metadata',
  [('Link sample\nmetadata', ''),
   ('Transform counts /\nbatch correction', 'ConQuR\nUMAP / HDBSCAN'),
   ('Group diagnostics /\ncohort selection', 'scikit-learn')],
  'Counts / CLR\n+ metadata'),
 ('Community\nanalysis',
  'Prepared tables',
  [('Diversity / ordination', 'scikit-bio / vegan'),
   ('Indicators / ASV clustermaps\ncontrasts / power', 'indicspecies / vegan'),
   ('Networks / genome links', 'SPIEC-EASI / BLAST+\nigraph / NetworkX')],
  'Community\nresults'),
 ('ASV–VOC\nintegration',
  'ASVs / VOCs / ISA\n+ sample metadata',
  [('Match samples\nand measurements', ''),
   ('Exploratory sample-level\nSpearman correlations', ''),
   ('Patient-level inference\npermutations / FDR', ''),
   ('ISA overlays /\nVOC–ASV clustermaps', '')],
  'Association tables\nand heatmaps'),
 ('Summaries &\nreports',
  'Results + audits',
  [('Read statistics /\nresult summaries', 'SeqKit'),
   ('Sankey / report', ''),
   ('Logs / checksums', 'Nextflow')],
  'Tables / figures\n+ provenance')]

def element(parent, tag, text=None, **attrs):
    child = ET.SubElement(parent, f'{{{NS}}}{tag}', {k.replace('_', '-'): str(v) for k, v in attrs.items()})
    child.text = text
    return child


def build(name, rows):
    width, gap = 1500, 175
    height = 335 + gap * len(rows) + (110 if rows is FULL else 400)
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
    text(750, 83, 'Mamba environments • Resource controls • Resumable execution • Audited outputs' if rows is FULL else 'Reads → filtered ASVs → prepared tables → analyses → reports', 21)
    panel(20, 120, 770, 130, *COLORS['input'])
    text(42, 152, 'Inputs', 25, 'start')
    text(42, 181, 'Paired amplicon reads, sample manifest, metadata and references\nConfigured biological / TECH / BIO sample classes\nOptional: VOC measurements and genome references' if rows is FULL else 'FASTQs • Manifest • Metadata • References\nOptional: controls • VOCs • Genomes', 21, 'start')
    panel(820, 120, 660, 100, '#CCCCCC', '#666666')
    for x, kind in zip([880, 1010, 1140, 1270, 1400], ['module', 'compute', 'data', 'input', 'output']):
        text(x, 150, kind.title(), 21)
        node(x, 183, kind)

    for i, (title, source, steps, output) in enumerate(rows):
        y = 345 + i * gap + (110 if rows is FULL and i > 7 else 0)
        if rows is BRIEF:
            y += {3: 120, 4: 320, 5: 400}.get(i, 0)
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
        if rows is BRIEF and i == 3:
            # Independent analyses consume the appropriate prepared cohort.
            for offset, (label, tool), cohort in zip(
                    [-120, 0, 120], steps, ['Full biological cohort', 'Selected cohort', 'Selected cohort']):
                yy = y + offset
                node(920, yy, 'compute')
                text(920, yy-(58 if '\n' in label else 35), label, 21)
                if tool:
                    text(920, yy+37, tool, 19)
                text(665, yy-12, cohort, 19)
                element(svg, 'path', d=f'M443.5 {y} H535 V{yy} H905', **{'class': 'wire'})
                element(svg, 'path', d=f'M935 {yy} H1260 V{y} H1356.5', **{'class': 'wire'})
        elif rows is BRIEF and i == 4:
            # Sample-level exploration and patient-level inference are separate
            # analyses of matched inputs, not successive statistical filters.
            for x, yy, (label, tool) in zip([650, 920, 920, 1160], [y, y-80, y+80, y], steps):
                node(x, yy, 'compute')
                text(x, yy-62, label, 21)
                if tool:
                    text(x, yy+37, tool, 19)
            wire(443.5, y, 635, y)
            for yy in [y-80, y+80]:
                element(svg, 'path', d=f'M665 {y} H780 V{yy} H905', **{'class': 'wire'})
                element(svg, 'path', d=f'M935 {yy} H1040 V{y} H1145', **{'class': 'wire'})
            wire(1175, y, 1356.5, y)
        elif title == 'Diversity\n& cohort selection':
            for yy, (label, tool) in zip([y-25, y+110], steps):
                node(900, yy, 'compute')
                text(900, yy-62, label, 21)
                if tool:
                    text(900, yy+37, tool, 19)
                element(svg, 'path', d=f'M443.5 {y} H550 V{yy} H885', **{'class': 'wire'})
                element(svg, 'path', d=f'M915 {yy} H1260 V{y} H1356.5', **{'class': 'wire'})
        else:
            last = 443.5
            for x, (label, tool) in zip(xs, steps):
                node(x, y, 'compute')
                text(x, y-66, label, 21)
                if tool:
                    text(x, y+37, tool, 19)
                wire(last, y, x-15, y)
                last = x+15
            wire(last, y, 1356.5, y)
        if i < len(rows)-1:
            extra = (110 if i == 7 else 0) if rows is FULL else {2: 120, 3: 200, 4: 80}.get(i, 0)
            wire(290, y+14.5, 290, y+gap+extra-14.5)
    text(750, height-20, 'Numbered modules group related operations; arrows summarize flow and enabled branches follow their dependencies.' if rows is FULL else 'Optional steps follow configuration. CLR: centered log-ratio • ISA: indicator species analysis.', 19)
    target = ROOT / 'docs/assets' / f'{name}.svg'
    ET.ElementTree(svg).write(target, encoding='unicode')
    cairosvg.svg2pdf(url=str(target), write_to=str(target.with_suffix('.pdf')))
    # Keep the documentation preview at the shared scale without rebuilding
    # the independent Mermaid diagrams.
    canvas_width = json.loads((ROOT / 'docs/diagrams/figures.json').read_text())['canvas_width']
    preview = ET.Element(f'{{{NS}}}svg', {
        'viewBox': f'0 0 {canvas_width} {height + 40}', 'width': str(canvas_width),
        'height': str(height + 40), 'role': 'img',
    })
    element(preview, 'title', name.replace('-', ' ').title())
    element(preview, 'rect', width=canvas_width, height=height + 40, fill='#ffffff')
    svg.set('x', f'{(canvas_width - width) / 2:g}')
    svg.set('y', '20')
    preview.append(svg)
    (ROOT / 'docs/assets/diagrams' / target.name).write_text(ET.tostring(preview, encoding='unicode') + '\n')
    print(f'Built {target.relative_to(ROOT)}, PDF and preview')


if __name__ == '__main__':
    build('workflow', FULL)
    build('workflow-brief', BRIEF)
