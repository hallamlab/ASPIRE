# ASPIRE

ASPIRE is a Nextflow DSL2 workflow for ASV generation, taxonomy assignment, decontamination, metadata-linked ASV summaries, ecological analyses, VOC association analyses, network/module analyses, and optional ASV-to-MAG linkage.

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

The test needs at least 8 CPU cores, 32 GB RAM and 50 GB storage as a starting provision. See the [reviewer guide](https://hallamlab-aspire.readthedocs.io/en/latest/reviewer-test.html) for inputs, expected results and coverage.

### Run your own study

```bash
cp asv_pipeline_nextflow.yml my_study.yml
```

Edit the input/output paths, references, sample metadata, assay settings and enabled modules in `my_study.yml`. Follow the [study walkthrough](https://hallamlab-aspire.readthedocs.io/en/latest/getting-started.html) and [input guide](https://hallamlab-aspire.readthedocs.io/en/latest/inputs.html); the complete template contains placeholders and study-specific examples.

```bash
./run_asv_pipeline.sh my_study.yml
```

Use the same command to resume. List supported restart points with `./run_asv_pipeline.sh --list-stages`. [Expert operations](https://hallamlab-aspire.readthedocs.io/en/latest/EXPERT_GUIDE.html) explains targeted reruns and cache management.

## Workflow

[![ASPIRE workflow from amplicon reads through ASVs, taxonomy, optional decontamination and analyses, genome links and integrated reports.](docs/assets/workflow.svg)](docs/assets/workflow.svg)

[Vector SVG](docs/assets/workflow.svg) · [PDF](docs/assets/workflow.pdf) · [Workflow and data-flow diagrams](https://hallamlab-aspire.readthedocs.io/en/latest/workflow.html)

## Full documentation

The [user guide](https://hallamlab-aspire.readthedocs.io/en/latest/index.html) covers [study preparation](https://hallamlab-aspire.readthedocs.io/en/latest/getting-started.html), [decontamination and analyses](https://hallamlab-aspire.readthedocs.io/en/latest/analyses.html), [configuration](https://hallamlab-aspire.readthedocs.io/en/latest/CONFIGURATION.html), [reports](https://hallamlab-aspire.readthedocs.io/en/latest/outputs.html), and [resuming or troubleshooting runs](https://hallamlab-aspire.readthedocs.io/en/latest/EXPERT_GUIDE.html).

Documentation source lives in `docs/`. Please use [GitHub issues](https://github.com/hallamlab/ASPIRE/issues) for bug reports and feature requests; include your Git revision and relevant logs.
