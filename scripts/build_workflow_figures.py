#!/usr/bin/env python3
"""Build the detailed workflow and reproduce the approved brief artwork exactly.

Requires CairoSVG and locally installed Trebuchet MS for the detailed PDF.
The SVG source preserves the Illustrator-exported brief lettering as outlines.
"""
from pathlib import Path
import shutil
import subprocess
from xml.etree import ElementTree as ET
import cairosvg

ROOT = Path(__file__).resolve().parents[1]
NS = 'http://www.w3.org/2000/svg'
ET.register_namespace('', NS)

# Stage label, software / methods, additional operations or resulting evidence.
SHARED = [
 ('Read QC', 'Cutadapt / fastp / VSEARCH / SeqKit',
  'Optional paired-primer removal and amplicon-family audit\nRead QC, paired-read merging, filtering and relabelling'),
 ('Denoise and remove chimeras', 'VSEARCH: dereplication / UNOISE / UCHIME',
  'Dereplicate and denoise reads; remove chimeras\nConstruct ASV sequences and the sample-by-ASV count matrix'),
 ('Alignment and taxonomic labelling', 'SINA / QIIME 2 / SILVA',
  'Align inferred sequences and assign full ASV taxonomy\nRetain taxonomy confidence and reference provenance'),
 ('Contamination and host removal', 'decontam / MITOMASTER / BLAST+ / SeqKit',
  'Biological depth QC; TECH / BIO prevalence tests when enabled\nCombine flagged ASVs; screen mitochondrial and contaminant references'),
 ('Abundance / prevalence filtering\nand ambiguous taxon removal', 'Combined count-table and FASTA filtering',
  'Apply configured abundance, prevalence and taxonomy rules\nRemove zero-count features; synchronize final counts and FASTA; audit exclusions'),
 ('Batch correction with metadata', 'ConQuR (optional) / metadata preparation',
  'Join sample metadata; construct long / wide tables and CLR representations\nCohort diagnostics and selection; synchronize corrected and raw analytical inputs'),
]
ASV = [
 ('Diversity estimation', 'Shannon / Bray–Curtis / Jaccard / PERMANOVA',
  'scikit-bio / vegan; ordination and collectors curves\nFull biological cohort retained; matched subsets where configured'),
 ('Clustering and outlier detection', 'UMAP / CLR / HDBSCAN',
  'Embedding, group diagnostics and outlier summaries\nSelected cohort feeds subsequent study-specific analyses'),
 ('Indicator species analysis', 'indicspecies; permutation tests / FDR',
  'Group and multigroup indicators; taxonomy-linked plots\nIndicator membership supplies ASV–VOC overlays'),
 ('Taxon clustering and mapping', 'Hierarchical clustering / heatmaps / vegan',
  'ASV clustermaps, taxonomy / group tests and paired contrasts\nPatient-count power simulations and study-design summaries'),
 ('Correlation network analysis', 'SPIEC-EASI / igraph: Leiden / NetworkX',
  'Filtered network inference, consensus modules and null topology\nGraph exports; optional genome-linked network overlays'),
]
VOC = [
 ('Correlation analysis', 'Spearman / p-values / BH-FDR',
  'Match ASVs, processed VOCs and sample metadata\nSample-type subsets and sample-level association tables'),
 ('Clustering and mapping', 'Hierarchical clustering / heatmaps',
  'VOC–ASV and VOC-profile clustermaps\nIndicator overlays and focused ISA–VOC comparisons'),
 ('Sample–VOC analysis', 'Z-score summaries',
  'Sample-level VOC patterns and annotated heatmaps\nCase summaries display mean and standard deviation'),
 ('Participant–VOC analysis', 'Mann–Whitney / permutations / BH-FDR',
  'Aggregate within participants; RA / CLR concordance\nPatient-level inference and cancer / control comparisons'),
]


def element(parent, tag, text=None, **attrs):
    child = ET.SubElement(parent, f'{{{NS}}}{tag}',
                          {k.replace('_', '-'): str(v) for k, v in attrs.items()})
    child.text = text
    return child


def build_detailed():
    width, height = 2300, 2860
    svg = ET.Element(f'{{{NS}}}svg', width=str(width), height=str(height),
                     viewBox=f'0 100 {width} {height}', role='img',
                     **{'aria-labelledby': 'title desc', 'data-workflow-style': 'aspire-approved-branched'})
    element(svg, 'title', 'ASPIRE detailed workflow', id='title')
    element(svg, 'desc', 'Shared read processing, taxonomy, contamination removal, filtering and metadata preparation split into ASV and optional ASV–VOC analyses. Genome linkage consumes filtered sequences independently of network inference. Results converge on summary reporting. Analysis spines group operations rather than impose a serial execution schedule.', id='desc')
    element(svg, 'rect', x=0, y=100, width=width, height=height, fill='#FFFFFF')
    element(svg, 'style', 'text{font-family:"Trebuchet MS",sans-serif;fill:#111111}.wire{fill:none;stroke:#111111;stroke-width:1.5;marker-end:url(#arrow)}.guide{fill:none;stroke:#888;stroke-width:1;stroke-dasharray:4 4}')
    defs = element(svg, 'defs')
    marker = element(defs, 'marker', id='arrow', viewBox='0 0 10 10', refX=10, refY=5,
                     markerWidth=7, markerHeight=7, orient='auto')
    element(marker, 'path', d='M0 0 L10 5 L0 10z', fill='#111111')

    def text(x, y, value, size=22, anchor='start'):
        parent = element(svg, 'text', x=x, y=y, font_size=size, text_anchor=anchor)
        for i, line in enumerate(value.split('\n')):
            element(parent, 'tspan', line, x=x, dy=0 if i == 0 else size*1.2)

    def node(x, y, kind):
        if kind == 'compute':
            element(svg, 'path', d=f'M{x} {y-14} L{x+14} {y} L{x} {y+14} L{x-14} {y}Z',
                    fill='#FFFFFF', stroke='#010101', stroke_width=1)
        else:
            fill, stroke = {'input': ('#DAE8F6', '#6D8EBF'),
                           'output': ('#D5E7D4', '#81B366')}[kind]
            element(svg, 'circle', cx=x, cy=y, r=13, fill=fill, stroke=stroke, stroke_width=3)

    def wire(path, arrow=True, dashed=False):
        attrs = {'class': 'wire'} if arrow else dict(fill='none', stroke='#111111', stroke_width=1.5)
        if dashed:
            attrs['stroke-dasharray'] = '5 5'
        element(svg, 'path', d=path, **attrs)

    def label_stage(x, y, label_x, stage, anchor):
        title, method, details = stage
        node(x, y, 'compute')
        title_lines = len(title.split('\n'))
        text(label_x, y-20-(title_lines-1)*28, title, 25, 'middle')
        text(label_x, y+18, method, 20, 'middle')
        text(label_x, y+53, details, 19, 'middle')
        edge = x-23 if anchor == 'end' else x+23
        element(svg, 'path', d=f'M{600 if anchor == "end" else 1710} {y} H{edge}', **{'class': 'guide'})

    node(980, 160, 'input'); text(840, 143, 'Sequencing data\nPaired FASTQs', 23, 'middle')
    node(1320, 160, 'input'); text(1510, 143, 'Auxiliary data\nMetadata and sample manifest', 23, 'middle')
    wire('M994 160 H1306', False)
    wire('M1150 160 V286')
    # Plain vertical legend and right-hand labels, matching the approved brief.
    for y, kind, name in [(460,'input','Input'),(560,'compute','Process'),(660,'output','Output')]:
        node(1450, y, kind)
        text(1505, y+7, name, 23)

    for i, stage in enumerate(SHARED):
        y = 300 + i*160
        title, method, details = stage
        node(1150, y, 'compute')
        text(760, y-20-(len(title.split('\n'))-1)*30, title, 26, 'middle')
        text(760, y+18, method, 20, 'middle')
        text(760, y+50, details, 19, 'middle')
        element(svg, 'path', d=f'M1105 {y} H1127', **{'class': 'guide'})
        wire(f'M1150 {y+14} V{y+146 if i<5 else 1236}')
    node(1150, 1250, 'output')
    text(1500, 1240, 'Filtered and corrected ASVs / count data', 24, 'middle')
    # Filtered ASV sequences also support genome linkage independently of networks.
    wire('M1136 1250 H1100 V2416')
    wire('M1150 1264 V1320 H650 V1416')
    wire('M1150 1320 H1660 V1416')
    text(600, 1360, 'ASV analysis', 28, 'end')
    text(1710, 1360, 'ASV–VOC analysis', 28)
    node(2070, 1320, 'input')
    text(2070, 1250, 'Processed VOC data\n(optional)', 22, 'middle')
    wire('M2056 1320 H1674')

    for i, stage in enumerate(ASV):
        y=1430+i*220
        label_stage(650,y,315,stage,'end')
        if i<4:wire(f'M650 {y+14} V{y+206}')
    for i, stage in enumerate(VOC):
        y=1430+i*220
        label_stage(1660,y,1980,stage,'start')
        if i<3:wire(f'M1660 {y+14} V{y+206}')
    # ISA evidence is a real VOC input; the long spines otherwise show conceptual grouping.
    wire('M664 1870 H1090', arrow=False, dashed=True)
    wire('M1110 1870 H1510 V1650 H1646', dashed=True)
    text(1400, 1845, 'ISA overlays', 19, 'middle')

    node(1370,2320,'input')
    text(1370,2255,'Genome references\n(optional)',22,'middle')
    wire('M1370 2334 V2370 H1100',False)
    node(1100,2430,'compute')
    text(1375,2410,'Genome linkage',25,'middle')
    text(1375,2450,'BLAST+ / genome-derived 16S references',19,'middle')
    text(1375,2480,'Identity / query-coverage and taxonomy gates\nUnique / ambiguous matches; genome QC metadata',18,'middle')
    wire('M1100 2444 V2616')
    node(1100,2630,'output')
    text(1375,2620,'Genome-linkage tables\nNetwork overlays / module anchors',20,'middle')
    wire('M650 2324 V2760 H1150',False)
    wire('M1100 2644 V2760',False)
    wire('M1660 2104 V2760 H1150',False)
    wire('M1150 2760 V2836')
    node(1150,2850,'output')
    text(1150,2900,'Summary reporting',28,'middle')
    text(1150,2940,'(figures and tables)',22,'middle')
    target=ROOT/'docs/assets/workflow.svg'
    ET.ElementTree(svg).write(target,encoding='unicode')
    family=subprocess.check_output(['fc-match','--format=%{family}','Trebuchet MS'],text=True)
    if 'Trebuchet MS' not in family:
        raise RuntimeError('Install Trebuchet MS before exporting the detailed PDF.')
    cairosvg.svg2pdf(url=str(target),write_to=str(target.with_suffix('.pdf')))
    # Use the same complete canvas, without introducing extra display margins.
    shutil.copyfile(target,ROOT/'docs/assets/diagrams/workflow.svg')
    print('Built detailed branched workflow SVG, PDF and preview')


def build_brief():
    source=ROOT/'docs/diagrams/workflow-brief.source.svg'
    target=ROOT/'docs/assets/workflow-brief.svg'
    shutil.copyfile(source,target)
    cairosvg.svg2pdf(url=str(source),write_to=str(target.with_suffix('.pdf')))
    shutil.copyfile(source,ROOT/'docs/assets/diagrams/workflow-brief.svg')
    print('Reproduced approved brief SVG, PDF and preview')


if __name__=='__main__':
    build_brief()
    build_detailed()
