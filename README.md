# ASPIRE: Amplicon Sequencing Profiler for Investigating Respiratory Ecosystems

ASPIRE is a Nextflow DSL2 workflow for ASV generation, taxonomy assignment, decontamination, metadata-linked ASV summaries, ecological analyses, association analyses linking ASVs with Volatile Organic Compounds (VOCs), network/module analyses, and optional ASV-to-MAG linkage.

**[Full user guide](https://hallamlab-aspire.readthedocs.io/en/latest/index.html)** · [Workflow test](https://hallamlab-aspire.readthedocs.io/en/latest/test.html) · [Issues and feature requests](https://github.com/hallamlab/ASPIRE/issues)

## Quick start

Use a 64-bit Linux system with Git and **Mamba** on `PATH`. The launcher creates its own Nextflow/Java controller and process environments. First-time setup needs internet access for packages and configured references. See [installation](https://hallamlab-aspire.readthedocs.io/en/latest/installation.html) for requirements.

```bash
git clone https://github.com/hallamlab/ASPIRE.git
cd ASPIRE
```

The control workflow runs **full taxonomy → independent TECH/BIO prevalence tests
→ union removal → reference screening → combined ASV filtering → metadata tables**.
Only biological samples face the 5,000-read inclusion cutoff; nonzero controls
are retained for testing. The mock uses **0.1% relative abundance in at least one
retained biological sample** plus **5% nonzero prevalence across biological
samples**, separately from decontam's score thresholds.

### Run the workflow test

Run the bundled **aspire-quickstart** fixture: approximately **2.5 MB**, 46
libraries and 218,000 read pairs. The wrapper builds the data, runs the full
benchmarked workflow and executes its validator:

```bash
./examples/run_quickstart.sh --output "$PWD/aspire-quickstart" --threads 4
```

The separate [**aspire-cami-mock** manuscript demonstration](https://hallamlab-aspire.readthedocs.io/en/latest/cami-mock.html)
contains 50 synthetic patients and 179 libraries. Both use the same full validator.

Success ends with `All mock-run checks passed.` Open
`aspire-quickstart/results/summary/report/ASPIRE_run_report.html` to review accounting, module outputs and execution logs. See [remote report viewing](https://hallamlab-aspire.readthedocs.io/en/latest/outputs.html#view-reports-on-a-remote-server) if running over SSH.

### Run your own study

```bash
cp asv_pipeline_nextflow.yml my_study.yml
```

Edit the input/output paths, references, sample metadata, assay settings and enabled modules in `my_study.yml`. Follow the [study walkthrough](https://hallamlab-aspire.readthedocs.io/en/latest/getting-started.html) and [input guide](https://hallamlab-aspire.readthedocs.io/en/latest/inputs.html); the complete template contains placeholders and study-specific examples.

```bash
./run_asv_pipeline.sh my_study.yml
```

Use the same command to resume. See the [user guide](https://hallamlab-aspire.readthedocs.io/) for resources, configuration, outputs and troubleshooting.

### Optional primer removal

Enable [Cutadapt primer trimming](docs/primer-trimming.md) to detect and remove paired primers before fastp. Set all four fixed fastp clipping values to zero when enabling this module. The primer audit records retained and discarded pairs, and a cohort check requires a single amplicon family before ASV construction.

## Workflow

[![ASPIRE workflow from amplicon reads through ASVs, taxonomy, optional decontamination and analyses, genome links and integrated reports.](docs/assets/workflow-brief.svg?v=mp-tools-20261008)](docs/assets/workflow-brief.svg)

[Vector SVG](docs/assets/workflow-brief.svg) · [PDF](docs/assets/workflow-brief.pdf) · [Workflow and data-flow diagrams](https://hallamlab-aspire.readthedocs.io/en/latest/workflow.html)

For the complete workflow, configuration reference and interpretation, use the **[Read the Docs guide](https://hallamlab-aspire.readthedocs.io/)**. Report bugs and request features through [GitHub issues](https://github.com/hallamlab/ASPIRE/issues).

## Cite ASPIRE

If you use ASPIRE, please cite:

> McLaughlin, R. J., Chen, S., Nag, A., Noonan, A. J. C., Bartolomeu, C., Borden, S. A., Lam, S., Myers, R., & Hallam, S. J. (2026). *ASPIRE: the Amplicon Sequencing Profiler for Investigating Respiratory Ecosystems*. bioRxiv, version 2 (12 August 2026). [DOI: 10.64898/2026.08.05.743000](https://doi.org/10.64898/2026.08.05.743000) · [Read version 2](https://www.biorxiv.org/content/10.64898/2026.08.05.743000v2).
