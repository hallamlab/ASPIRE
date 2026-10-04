# ASPIRE: Amplicon Sequencing Profiler for Investigating Respiratory Ecosystems

ASPIRE is a Nextflow DSL2 workflow for ASV generation, taxonomy assignment, decontamination, metadata-linked ASV summaries, ecological analyses, association analyses linking ASVs with Volatile Organic Compounds (VOCs), network/module analyses, and optional ASV-to-MAG linkage.

**[Full user guide](https://hallamlab-aspire.readthedocs.io/en/latest/index.html)** · [Reviewer test](https://hallamlab-aspire.readthedocs.io/en/latest/reviewer-test.html) · [Issues and feature requests](https://github.com/hallamlab/ASPIRE/issues)

## Quick start

Use a 64-bit Linux system with Git and **Mamba** on `PATH`. The launcher creates its own Nextflow/Java controller and process environments. First-time setup needs internet access for packages and configured references. See [installation](https://hallamlab-aspire.readthedocs.io/en/latest/installation.html) for requirements.

```bash
git clone https://github.com/hallamlab/ASPIRE.git
cd ASPIRE
```

### Run the reviewer test

Download and extract the [public mock dataset (Zenodo DOI: 10.5281/zenodo.22906294)](https://doi.org/10.5281/zenodo.22906294). The data is downloaded separately; configuration and validation scripts are included in the repository. Edit the two paths below, then run from the ASPIRE directory:

```bash
DATASET=/absolute/path/to/mock_dataset
RESULTS=/absolute/path/to/aspire_mock_output

./examples/configure_mock_run.sh \
  --dataset "$DATASET" \
  --output "$RESULTS" \
  --config-out mock_run.generated.yml
./run_asv_pipeline.sh mock_run.generated.yml --no-resume
./examples/validate_mock_run.sh --dataset "$DATASET" --results "$RESULTS"
```

Success ends with `All mock-run checks passed.` Open
`<RESULTS>/summary/report/ASPIRE_run_report.html` to review accounting, module outputs and execution logs. See [remote report viewing](https://hallamlab-aspire.readthedocs.io/en/latest/outputs.html#view-reports-on-a-remote-server) if running over SSH.

### Run your own study

```bash
cp asv_pipeline_nextflow.yml my_study.yml
```

Edit the input/output paths, references, sample metadata, assay settings and enabled modules in `my_study.yml`. Follow the [study walkthrough](https://hallamlab-aspire.readthedocs.io/en/latest/getting-started.html) and [input guide](https://hallamlab-aspire.readthedocs.io/en/latest/inputs.html); the complete template contains placeholders and study-specific examples.

```bash
./run_asv_pipeline.sh my_study.yml
```

Use the same command to resume. See the [user guide](https://hallamlab-aspire.readthedocs.io/) for resources, configuration, outputs and troubleshooting.

## Workflow

[![ASPIRE workflow from amplicon reads through ASVs, taxonomy, optional decontamination and analyses, genome links and integrated reports.](docs/assets/workflow-brief.svg?v=6cbea7b6e524)](docs/assets/workflow-brief.svg)

[Vector SVG](docs/assets/workflow-brief.svg) · [PDF](docs/assets/workflow-brief.pdf) · [Workflow and data-flow diagrams](https://hallamlab-aspire.readthedocs.io/en/latest/workflow.html)

For the complete workflow, configuration reference and interpretation, use the **[Read the Docs guide](https://hallamlab-aspire.readthedocs.io/)**. Report bugs and request features through [GitHub issues](https://github.com/hallamlab/ASPIRE/issues).

## Cite ASPIRE

If you use ASPIRE, please cite:

> McLaughlin, R. J., Chen, S., Nag, A., Noonan, A. J. C., Bartolomeu, C., Borden, S. A., Lam, S., Myers, R., & Hallam, S. J. (2026). *ASPIRE: the Amplicon Sequencing Profiler for Investigating Respiratory Ecosystems*. bioRxiv, version 2 (12 August 2026). [DOI: 10.64898/2026.08.05.743000](https://doi.org/10.64898/2026.08.05.743000) · [Read version 2](https://www.biorxiv.org/content/10.64898/2026.08.05.743000v2).
