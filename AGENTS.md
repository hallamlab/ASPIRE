# Repository instructions

## Workflow figure aesthetic — standing user preference

Read `docs/WORKFLOW_STYLE.md` before changing workflow diagrams, their generators,
exports or documentation embeds. Preserve the MetaPathways nodal publication
style: numbered square module spine, diamond compute nodes, circular input /
data / output nodes, Times serif labels, thin black connectors and the shared
blue / gray / green palette. Update scientific content within that layout.
Do not replace these figures with stacked text cards or generic box flowcharts.
A style change requires an explicit user request.

Rebuild source-linked SVG / PDF / preview artifacts together, run
`python3 scripts/check_workflow_style.py`, and visually inspect the result.
This governs workflow figures; preserve the separately configured study palettes
for scientific plots. Apply this convention to new workflow figures too.

ASPIRE-specific user override (2026-10-14): preserve outlined Myriad Pro in the brief artwork;
use the approved Trebuchet MS substitute for other workflow text.
Reproduce the approved `docs/diagrams/workflow-brief.source.svg` exactly; expand the brief shared-chain / ASV / ASV–VOC structure
in the detailed workflow. See the override in docs/WORKFLOW_STYLE.md.
