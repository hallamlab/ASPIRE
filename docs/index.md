# ASPIRE: Amplicon Sequencing Profiler for Investigating Respiratory Ecosystems

Build amplicon sequence variants (ASVs), assign taxonomy, screen non-targets and contaminants, and connect microbial community profiles to sample metadata, Volatile Organic Compounds (VOCs) and genome references.

**Start here:** [install ASPIRE](installation.md), then run the [public workflow test](test.md). The launcher manages Mamba environments and Nextflow. Use the [study walkthrough](getting-started.md) when you are ready to analyze your own samples.

```{container} aspire-primary-workflow
[![ASPIRE TECH/BIO workflow: full taxonomy, control decontamination, reference and biological feature filtering, metadata and analyses.](assets/workflow-brief.svg)](assets/workflow-brief.svg)
```

[SVG](assets/workflow-brief.svg) · [PDF](assets/workflow-brief.pdf) · [Workflow details and diagrams](workflow.md)

Arrows trace the main data path. Control decontamination precedes biological feature filtering and metadata construction. Optional analyses depend on configuration.

```{toctree}
:maxdepth: 1
:caption: Getting started

installation
test
cami-mock
getting-started
inputs
```

```{toctree}
:maxdepth: 1
:caption: Run and understand analyses

workflow
decontamination
analyses
outputs
VOC_STATISTICS
```

```{toctree}
:maxdepth: 1
:caption: Configuration and advanced use

CONFIGURATION
CONFIG_PARAMETERS
PROCESS_REFERENCE
EXPERT_GUIDE
troubleshooting
SPARK_FUNCTIONALITY
documentation
```

[Source code](https://github.com/hallamlab/ASPIRE) · [Issues and feature requests](https://github.com/hallamlab/ASPIRE/issues)

## Cite ASPIRE

If you use ASPIRE, please cite:

> McLaughlin, R. J., Chen, S., Nag, A., Noonan, A. J. C., Bartolomeu, C., Borden, S. A., Lam, S., Myers, R., & Hallam, S. J. (2026). *ASPIRE: the Amplicon Sequencing Profiler for Investigating Respiratory Ecosystems*. bioRxiv, version 2 (12 August 2026). [DOI: 10.64898/2026.08.05.743000](https://doi.org/10.64898/2026.08.05.743000) · [Read version 2](https://www.biorxiv.org/content/10.64898/2026.08.05.743000v2).
