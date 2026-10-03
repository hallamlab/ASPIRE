# Installation

ASPIRE runs from its GitHub checkout. Mamba builds the controller and per-process environments automatically; you do not need to install Nextflow, Java or every analysis tool separately.

```bash
git clone https://github.com/hallamlab/ASPIRE.git
cd ASPIRE
```

Continue with the [reviewer test](reviewer-test.md) or [your own study](getting-started.md).

## Requirements

ASPIRE is developed for a 64-bit Linux environment. Before starting, install:

- Bash and standard GNU command-line utilities.
- Conda or Mamba, with `mamba` available on `PATH`.
- Git for obtaining and identifying the workflow revision.
- Internet access on the first run to solve Conda environments and download
  configured SINA and QIIME2/SILVA references. Fully offline runs require
  pre-populated package caches and local reference paths.

The wrapper creates a repository-local controller environment containing
Nextflow and its Java runtime, then creates process-specific environments from
the committed YAML definitions. Users should not manually combine all process
dependencies into one environment. ASPIRE isolates the run's package cache and
serializes Conda environment creation to prevent concurrent repodata-lock
failures. Environment solves use strict channel priority to avoid pathological
cross-channel backtracking, and module-specific environments avoid the legacy
all-in-one dependency search space. Builds are terminated after 30 minutes by
default rather than hanging indefinitely; set `ASPIRE_MAMBA_BUILD_TIMEOUT` only
when a slower package source is expected. Analysis tasks remain parallel.

Resource needs depend on sample count and sequencing depth. For the complete
mock benchmark, provision at least 8 CPU cores, 32 GB RAM, and 50 GB of free
storage for input data, Conda environments, Nextflow work files, downloaded
references, and final outputs. `core.resources.threads` controls per-task CPU use; it
does not limit the total storage used by cached tasks.

Verify the entrypoint prerequisites:

```bash
command -v bash
command -v git
command -v mamba
```

If `mamba` is unavailable, install a current Miniforge distribution from
<https://github.com/conda-forge/miniforge> and open a new shell before running
ASPIRE.


## Key Files

- `run_asv_pipeline.sh`: main wrapper for routine runs.
- `asv_pipeline.nf`: current Nextflow workflow.
- `asv_pipeline_nextflow.yml`: authoritative complete config template with every supported parameter.
- `examples/mock.local.yml`: full-module template used by the portable mock config generator.
- `examples/configure_mock_run.sh`: validates a supplied mock fixture and writes a machine-local YAML.
- `examples/validate_mock_run.sh`: validates the completed mock run against its truth contract.
- [reviewer guide](reviewer-test.md): expanded mock test instructions.
- `docs/CONFIGURATION.md`: configuration and dependency reference.
- `docs/CONFIG_PARAMETERS.md`: exhaustive key-by-key definitions and template values.
- `docs/EXPERT_GUIDE.md`: restart, cache, resource, and diagnostic guidance.
- `docs/PROCESS_REFERENCE.md`: process purpose, inputs, and published outputs.
- `docs/SPARK_FUNCTIONALITY.md`: manuscript-to-workflow functionality crosswalk.
- `processes/`: scripts and conda environment YAMLs used by individual stages.
