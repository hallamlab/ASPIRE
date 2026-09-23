# ASPIRE Expert Guide

This guide covers execution and diagnosis on the `main` branch. Start with the
[README](../README.md) for installation and quick starts, the
[configuration reference](CONFIGURATION.md) for YAML semantics, and the
[process reference](PROCESS_REFERENCE.md) for process I/O.

## Preserve Run Identity

Record the repository revision and keep one configuration per analytical run:

```bash
git rev-parse HEAD
cp asv_pipeline_nextflow.yml my_study.yml
./run_asv_pipeline.sh my_study.yml
```

Use a new `core.paths.output_dir` for an independent parameter experiment. Reusing
an output directory means the wrapper may reuse its `.aspire` runtime cache and
will replace the public module tree after a successful run.

## Resume or Start Fresh

The normal command resumes the best retained run history automatically:

```bash
./run_asv_pipeline.sh my_study.yml
```

Run without cache reuse:

```bash
./run_asv_pipeline.sh my_study.yml --no-resume
```

List valid stages and force the selected stage plus later registered stages to
rerun while allowing unaffected upstream work to remain cached:

```bash
./run_asv_pipeline.sh --list-stages
./run_asv_pipeline.sh my_study.yml --rerun-from standard:FILTER_COUNTS
```

The live ANSI dashboard uses the real `core:PROCESS`, `standard:PROCESS`, and
`optional:PROCESS` Nextflow scopes. These are execution identities rather than
post-processed display labels. The wrapper defaults `TERMINAL_WIDTH` to 160 so
the qualified names are not shortened with an ellipsis; an explicitly supplied
environment value takes precedence.

Choose a specific retained Nextflow run only when automatic selection is not
appropriate:

```bash
cd /path/to/ASPIRE
./.controller_env/bin/nextflow log -q
./run_asv_pipeline.sh my_study.yml --resume-run RUN_NAME
```

`--resume-policy last-with-tasks` is the default. `latest` is available for a
deliberate preference for the newest history entry.

## Choosing the Earliest Rerun Stage

Start at the earliest process whose inputs, parameters, script, or environment
changed. Common examples:

| Change | Suggested stage |
|---|---|
| FASTQ trimming or merge parameters | `FASTP_QC` or `MERGE_READS` |
| technical sample/ASV threshold | `FILTER_TABLE` |
| taxonomy reference | `TAXONOMY` |
| host/non-target rules | `FILTER_COUNTS` |
| control subtraction or metadata cohort | `PLOT_METADATA` |
| three-tier thresholds | `THREE_TIER_DECONTAM` |
| ISA settings | `INDICSPECIES` |
| VOC correlation/plot thresholds | `VOC_CORRELATION` |
| network cohort or inference settings | `SPIECEASI` |
| topology null count/seed | `NETWORK_TOPOLOGY` |

The stage registry uses `core:PROCESS`, `standard:PROCESS`, and
`optional:PROCESS` identifiers. Bare process names remain accepted for existing
automation. The registry is a rerun boundary list, not a promise that every later
independent branch biologically depends on the selected process. Nextflow still
uses the declared DAG and task hashes.

## Runtime and Storage

By default, runtime state is retained under:

```text
<output_dir>/.aspire/
├── nf_work/
├── conda_cache/
└── publication_staging/
```

Keep `core.paths.keep_runtime_dir: true` until validation and any planned reruns are
complete. `core.paths.work_dir` can move task work to fast scratch or network
storage; `core.paths.conda_cache_dir` can move environments independently. Do not
point concurrent runs at the same Conda cache: the controller serializes access
and may reject unsafe overlap.

Disk use includes input staging, per-task work products, Conda environments,
references and published output. `core.resources.threads` controls per-task CPU use,
not aggregate disk or the number of concurrently scheduled independent tasks.

## Nextflow Options

Pass native options after `--` so the wrapper still performs environment setup,
logging and atomic output publication:

```bash
./run_asv_pipeline.sh my_study.yml -- \
  -with-report custom_report.html \
  -with-trace custom_trace.tsv
```

Direct `nextflow run` is useful for workflow debugging but bypasses part of the
supported wrapper lifecycle and should not be the normal user entrypoint.

## Reading a Failed Task

The terminal error reports a work directory. Inspect its generated command and
logs:

```bash
cd /reported/nextflow/work/directory
sed -n '1,240p' .command.sh
sed -n '1,240p' .command.err
sed -n '1,240p' .command.out
```

After a successful run, equivalent task records are copied under
`<output_dir>/logs/tasks/`, and `<output_dir>/logs/task_execution.tsv` maps
process names and hashes to statuses and durations. The controller stream is
stored in `<output_dir>/logs/controller.log`.

## Environment and Network Failures

The wrapper creates `.controller_env` and process-specific Conda environments.
Messages about a serialized mamba slot indicate controlled environment creation.
If a solve or download fails, fix connectivity or package availability and run
the same command again; completed work remains resumable.

MITOMASTER is an external service. For a deliberately offline analysis, set
`standard.mito.run_mitomaster: false` and document that choice. Do not silently treat an
API failure as a valid no-hit response. Local mitochondrial and contaminant
BLAST screens can still run when compatible database prefixes or FASTAs are
configured.

## Analytical Review

- Confirm `summary/tables/run_manifest.tsv` contains the intended samples.
- Review read/sample loss before interpreting group differences.
- Treat `.micro.tsv`, `.decon.tsv`, and raw matrices as audit intermediates;
  identify the configured final downstream table.
- Review control classifications and read loss when three-tier decontamination
  is enabled.
- Review effect sizes and adjusted p-values, not only nominal significance.
- Treat network edges and ASV-VOC correlations as associations.
- Confirm repeated observations use the configured participant/block column.
- Check `module_output_manifest.tsv` and reference checksums before archiving.

## Extending the Workflow

New processes should have a process-specific environment, explicit inputs and
outputs, configuration validation, stable publication paths, and tests. Add the
process to `run_asv_pipeline.sh --list-stages` when it is a supported rerun
boundary. Update `PROCESS_REFERENCE.md`, `CONFIGURATION.md`, the full YAML
template, mock configuration and mock validator together so the reviewer path
continues to exercise the advertised functionality.
