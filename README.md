# ASPIRE

ASPIRE is a Nextflow DSL2 workflow for ASV generation, taxonomy assignment, decontamination, metadata-linked ASV summaries, ecological analyses, VOC association analyses, network/module analyses, and optional ASV-to-MAG linkage.

## Quick start

Use a 64-bit Linux system with Git and **Mamba** on `PATH`. The launcher creates its own Nextflow/Java controller and process environments. First-time setup needs internet access for packages and configured references. See [installation](docs/installation.md) for requirements.

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
`<RESULTS>/summary/report/ASPIRE_run_report.html` to review accounting, module outputs and execution logs. See [remote report viewing](docs/outputs.md#view-reports-on-a-remote-server) if running over SSH.

The fixture has 60 biological samples plus four extraction controls. It exercises three-tier decontamination and the benchmarked statistical/network branches; optional ASV-to-MAG linkage needs separate genome inputs. See the [full reviewer guide](examples/MOCK_DATASET_TESTING.md) for inputs, expected outputs and validation details. Allow at least 8 CPU cores, 32 GB RAM and 50 GB storage as a starting provision for this test; runtime and storage depend on input depth and package caches.

### Run your own study

```bash
cp asv_pipeline_nextflow.yml my_study.yml
```

Edit the input/output paths, references, sample metadata, assay settings and enabled modules in `my_study.yml`. Follow the [study walkthrough](docs/getting-started.md) and [input guide](docs/inputs.md); the complete template contains placeholders and study-specific examples.

```bash
./run_asv_pipeline.sh my_study.yml
```

Use the same command to resume. List supported restart points with `./run_asv_pipeline.sh --list-stages`. [Expert operations](docs/EXPERT_GUIDE.md) explains targeted reruns and cache management.

## Workflow

[![ASPIRE workflow from amplicon reads through ASVs, taxonomy, optional decontamination and analyses, genome links and integrated reports.](docs/assets/workflow.svg)](docs/assets/workflow.svg)

[Vector SVG](docs/assets/workflow.svg) · [PDF](docs/assets/workflow.pdf) · [Workflow and data-flow diagrams](docs/workflow.md)

Numbered modules are conceptual groups, not a serial execution schedule. Configuration separates `core` ASV construction/taxonomy, `standard` final-table preparation, and `optional` analyses. Three-tier decontamination uses control-bearing raw counts to score contaminants and applies the flags to downstream microbial tables. Module selection and the available metadata determine which branches run.

[Brief appnote figure: SVG](docs/assets/workflow-brief.svg) · [PDF](docs/assets/workflow-brief.pdf)

## Documentation

| Goal | Guide |
|---|---|
| Install and test | [Installation](docs/installation.md) · [Reviewer test](examples/MOCK_DATASET_TESTING.md) |
| Prepare a study | [Study walkthrough](docs/getting-started.md) · [Inputs](docs/inputs.md) |
| Understand filtering | [Decontamination](docs/decontamination.md) · [Analysis guide](docs/analyses.md) |
| Configure every module | [Configuration guide](docs/CONFIGURATION.md) · [Parameter catalogue](docs/CONFIG_PARAMETERS.md) |
| Trace the workflow | [Diagrams](docs/workflow.md) · [Process I/O reference](docs/PROCESS_REFERENCE.md) |
| Review results | [Outputs and report](docs/outputs.md) · [VOC statistics](docs/VOC_STATISTICS.md) |
| Resume or troubleshoot | [Expert guide](docs/EXPERT_GUIDE.md) · [Troubleshooting](docs/troubleshooting.md) |
| Check manuscript methods | [SPARK functionality](docs/SPARK_FUNCTIONALITY.md) |
| Build or publish the docs | [Documentation maintenance](docs/documentation.md) |

The same guides build as a Read the Docs site. [GitHub issues](https://github.com/hallamlab/ASPIRE/issues) welcomes bug reports and feature requests; include the Git revision, relevant configuration and error logs.
