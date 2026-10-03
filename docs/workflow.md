# Workflow, architecture and data flow

[![Complete ASPIRE workflow](assets/workflow.svg)](assets/workflow.svg)

[Download SVG](assets/workflow.svg) · [Download PDF](assets/workflow.pdf)

The figure groups the complete workflow into conceptual modules. Numbered rows
are not a serial execution schedule. Optional analyses run only when enabled
and their dependencies are available. The [process reference](PROCESS_REFERENCE.md)
lists the exact task names, inputs and outputs.

## Conceptual workflow

```mermaid
%%{init: {"theme":"base","themeVariables":{"fontFamily":"Times New Roman, Times, serif","primaryColor":"#F5F5F5","primaryBorderColor":"#666666","lineColor":"#111111"}}}%%
flowchart TD
    I[FASTQs and sample manifest] --> C[Read QC and ASV construction]
    R[Sequence and taxonomy references] --> T[Alignment and taxonomy]
    C --> T
    T --> F[Reference-based non-target filtering]
    F --> M[Metadata-linked microbial tables]
    S[Sample metadata and controls] --> M
    M --> D{Three-tier decontamination enabled?}
    D -->|Yes| DC[Control-based scores and filtering]
    D -->|No| A[Analysis table preparation]
    C -->|Raw counts retaining controls| DC
    S --> DC
    DC --> A
    A --> E[Ecology, indicators and associations]
    A --> N[Networks, modules and topology]
    C -->|ASV sequences| L[Optional ASV-to-genome linkage]
    G[Genome and barrnap references] --> L
    L --> N
    E --> O[Tables, plots and run report]
    N --> O
    classDef input fill:#DAE8FC,stroke:#6C8EBF,color:#111111
    classDef output fill:#D5E8D4,stroke:#82B366,color:#111111
    class I,R,S,G input
    class O output
```

Three-tier filtering scores contaminants from the raw control-bearing counts,
then filters the microbial analysis tables. Batch correction and proposed group
labels are separately configurable. Network edges and ASV–VOC correlations are
associations, rather than evidence of causal relationships.

## Software architecture

```mermaid
%%{init: {"theme":"base","themeVariables":{"fontFamily":"Times New Roman, Times, serif","primaryColor":"#F5F5F5","primaryBorderColor":"#666666","lineColor":"#111111"}}}%%
flowchart TD
    U[User: YAML configuration and manifest] --> W[run_asv_pipeline.sh]
    W --> C[Mamba controller environment: Nextflow and Java]
    C --> N[Nextflow dependency graph]
    N --> E[Per-process Mamba environments]
    E --> P[Tools and ASPIRE Python / R scripts]
    P --> S[Runtime staging and task work]
    S -->|After workflow success| O[Output organizer]
    O --> R[Module tables and plots]
    O --> H[HTML run report and master summaries]
    O --> L[Logs, inventories and checksums]
    S -. Retained cache for resume .-> N
    classDef input fill:#DAE8FC,stroke:#6C8EBF,color:#111111
    classDef output fill:#D5E8D4,stroke:#82B366,color:#111111
    class U input
    class R,H,L output
```

The supported entrypoint is the wrapper. It handles environment setup, resume,
logging and final publication. Independent tasks can overlap when resources
permit. Per-task threads are not a global concurrency limit. See
[expert execution guidance](EXPERT_GUIDE.md) before changing runtime locations
or passing native Nextflow options.

## Data flow and the selected analysis table

```mermaid
%%{init: {"theme":"base","themeVariables":{"fontFamily":"Times New Roman, Times, serif","primaryColor":"#F5F5F5","primaryBorderColor":"#666666","lineColor":"#111111"}}}%%
flowchart TD
    RAW[Raw ASV counts including controls] --> TECH[Technical sample / ASV filters]
    TECH --> TARGET[ASV_target.tsv: non-target-filtered microbial counts]
    TARGET --> BASE[Long ASV annotations and wide counts]
    RAW --> SCORE[Optional contamination scoring]
    SCORE --> FILTER[Apply contamination flags to long / wide tables]
    BASE --> FILTER
    BASE -->|Filtering disabled| SELECT[Selected base tables]
    FILTER --> SELECT
    SELECT --> GROUP[Optional group diagnostics / validated labels]
    GROUP --> BATCH[Optional batch correction and count selection]
    BATCH --> FINAL[Selected counts and synchronized long annotations]
    FINAL --> STATS[Statistical and network modules]
    STATS --> REPORT[Module files and integrated report]
    classDef input fill:#DAE8FC,stroke:#6C8EBF,color:#111111
    classDef output fill:#D5E8D4,stroke:#82B366,color:#111111
    class RAW input
    class FINAL,REPORT output
```

The diagrams summarize the analytical table path. Individual descriptive
products can use their own declared metadata or intermediate inputs; consult
the process I/O reference before treating every report panel as the same cohort.
ASV-to-genome linkage uses the filtered ASV sequences and supplied genome
references; network overlays combine those links with the selected analytical
network. Exact file paths and checksums are recorded in the output inventories.

## Brief appnote panel

[![ASPIRE compact workflow overview](assets/workflow-brief.svg)](assets/workflow-brief.svg)

[Brief SVG](assets/workflow-brief.svg) · [Brief PDF](assets/workflow-brief.pdf)

This smaller figure groups the workflow into sequence processing, filtering, optional analytical branches and reporting. Use the complete figure above for individual tools and the three-tier filtering details.
