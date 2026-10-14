---
orphan: true
---

# Shared workflow figure style: MP nodal v1

## ASPIRE artwork and typography override

The user's explicit ASPIRE instruction supersedes the shared typography and
brief-layout defaults below: ASPIRE workflow text uses **Trebuchet MS**, the user-approved locally available
humanist sans-serif substitute for Myriad Pro. The brief artwork retains its
original outlined Myriad Pro lettering.
The brief figure is the approved Illustrator SVG at
`docs/diagrams/workflow-brief.source.svg`. The build copies it byte-for-byte,
including its outlined lettering, geometry, colors, connectors and canvas.
Do not recreate it from the legacy BRIEF row definitions or add preview margins.
The detailed workflow expands the brief shared processing chain and its ASV /
ASV–VOC branches. Genome linkage takes filtered sequences independently of
network inference; the analysis spines show conceptual grouping, not job order.
Font files are not distributed by this repository; Trebuchet MS must be
installed locally when rendering the detailed PDF and supplementary diagrams.
This ASPIRE layout override also replaces the numbered-square / horizontal-row
layout below. Detailed figures also match the brief’s unboxed vertical legend,
centered label blocks, white process diamonds and input / output colors. Do not
add a title banner or explanatory footer; place those details in captions.


The MetaPathways publication workflow is the visual reference for workflow
figures across MetaPathways, ASPIRE, SCARAB, BASINS, DECOI and related repositories.
This is the user's standing preference. Content updates must preserve this
visual grammar; a different layout requires an explicit user request.

## Layout and node grammar

- White canvas; a restrained gray title / execution header.
- Inputs and the node legend sit immediately below the header.
- Numbered square module nodes form a connected vertical spine on the left.
  Module labels sit to its left. A number groups related operations; it does
  not assert that independent analyses execute serially.
- Each module reads left to right: circular inputs, diamond compute steps,
  circular data / outputs. Put short process labels above compute diamonds
  and concise software / major-library names directly below each diamond.
  Omit module-wide explanatory notes; keep those details in the accompanying
  documentation. Preserve generous whitespace.
- Show actual forks and joins with clean orthogonal connectors. Label optional
  branches explicitly. Keep cohort-specific paths scientifically accurate.
- Use thin black connectors with small arrowheads and heavier node outlines.
  Avoid lines through text, overlapping labels and truncated exports.
- Do not replace the node graph with stacked paragraph cards, a text table,
  a dashboard, gradients, decorative icons or a generic box flowchart.

## Typography and palette

| Element | Style |
|---|---|
| Labels | Times New Roman, Times, serif; black `#111111` |
| Canvas | White `#FFFFFF` |
| Header / legend | Gray `#CCCCCC`, outline `#666666` |
| Compute / intermediate data | `#F5F5F5`, outline `#666666` |
| Module square | `#F5F5F5`, outline `#111111` |
| Input circle | Blue `#DAE8FC`, outline `#6C8EBF` |
| Output circle | Green `#D5E8D4`, outline `#82B366` |
| Connectors | Black, approximately 1.5 SVG units |
| Node outline | Approximately 3 SVG units |

At publication scale, use approximately 20–22-unit node labels and 26–29-unit
module / title text. Adjust the canvas and whitespace to fit the content;
do not shrink the entire figure to accommodate longer paragraphs. Full and
brief figures use the same grammar. Biological plot palettes are separate
study settings and are not replaced by these workflow colors.

## Documentation progression

1. Put the brief nodal overview on the documentation landing page and link
   directly to the detailed workflow page.
2. Lead that detailed page with a separate, expanded nodal figure in the same
   visual style. Provide a full-size SVG link and a PDF where exported.
3. Place the supporting architecture, dependency and data-flow diagrams below
   the large figure, followed by process explanations and citations.
4. Keep a link back to the brief overview beside the detailed figure’s downloads.

The brief and detailed figures must be distinct levels of detail. Supporting
Mermaid diagrams complement the expanded nodal figure rather than substitute
for it. Apply this progression to new documentation guides as well.

## Sources, exports and review

Maintain an editable SVG or a deterministic generator as the source of truth.
Update the generator when one exists, then rebuild SVGs, PDFs and documentation
previews together. Never copy a different diagram over a generated preview.
Preserve source-to-preview identity and shared preview scale where configured.

Detailed software dependency diagrams may retain Mermaid where they already
serve as technical supplements, using the same serif typography and restrained
palette. They do not replace the primary full or brief MP nodal workflow.

Before completing a workflow edit:

1. Run `python3 scripts/check_workflow_style.py`.
2. Render and inspect both full and brief figures at readable size. Check the
   documentation page and exported figure for clipping, overlap and spacing.
3. Build the documentation with warnings treated as errors where supported.
4. Confirm that module content and optional branches match the pipeline.

The automated check rejects loss of the nodal grammar, serif font or palette.
It complements visual review; it cannot certify layout quality or scientific
correctness. The canonical visual reference is MetaPathways
`docs/assets/workflow-main.svg` (the appnote module-and-node layout).
