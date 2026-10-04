# Workflow, architecture and data flow

```{container} aspire-primary-workflow
[![Complete ASPIRE workflow](assets/workflow.svg)](assets/workflow.svg)
```

[Download SVG](assets/workflow.svg) · [Download PDF](assets/workflow.pdf)

The figure follows the appnote’s overall theme—read/ASV processing, taxonomy and background removal, metadata and quality assessment, community diversity, indicators/networks, and ASV–VOC integration—while expanding those areas into detailed modules and retaining the additional optional branches. Numbered rows
are not a serial execution schedule. Optional analyses run only when enabled
and their dependencies are available. The [process reference](PROCESS_REFERENCE.md)
lists the exact task names, inputs and outputs.

## Conceptual workflow

[![Workflow 1](assets/diagrams/workflow-1.svg)](assets/diagrams/workflow-1.svg)

[Zoom diagram](assets/diagrams/workflow-1.svg) · [Mermaid source](diagrams/workflow-1.mmd)

Three-tier filtering scores contaminants from the raw control-bearing counts,
then filters the microbial analysis tables. Batch correction and proposed group
labels are separately configurable. Network edges and ASV–VOC correlations are
associations, rather than evidence of causal relationships.

## Software architecture

[![Workflow 2](assets/diagrams/workflow-2.svg)](assets/diagrams/workflow-2.svg)

[Zoom diagram](assets/diagrams/workflow-2.svg) · [Mermaid source](diagrams/workflow-2.mmd)

The supported entrypoint is the wrapper. It handles environment setup, resume,
logging and final publication. Independent tasks can overlap when resources
permit. Per-task threads are not a global concurrency limit. See
[expert execution guidance](EXPERT_GUIDE.md) before changing runtime locations
or passing native Nextflow options.

## Data flow and the selected analysis table

[![Workflow 3](assets/diagrams/workflow-3.svg)](assets/diagrams/workflow-3.svg)

[Zoom diagram](assets/diagrams/workflow-3.svg) · [Mermaid source](diagrams/workflow-3.mmd)

The diagrams summarize the analytical table path. Individual descriptive
products can use their own declared metadata or intermediate inputs; consult
the process I/O reference before treating every report panel as the same cohort.
ASV-to-genome linkage uses the filtered ASV sequences and supplied genome
references; network overlays combine those links with the selected analytical
network. Exact file paths and checksums are recorded in the output inventories.

## Brief appnote panel

[![ASPIRE compact workflow overview](assets/diagrams/workflow-brief.svg)](assets/workflow-brief.svg)

[Brief SVG](assets/workflow-brief.svg) · [Brief PDF](assets/workflow-brief.pdf)

The brief panel follows the appnote methods: read/ASV processing, taxonomy and background removal, metadata and quality assessment, community diversity, indicators/networks, and ASV–VOC integration. It highlights the statistical methods; the complete figure includes additional optional branches such as genome linkage. These are conceptual groups: three-tier decontamination runs after metadata-table construction, and SINA trimming follows ASV inference in the code.
