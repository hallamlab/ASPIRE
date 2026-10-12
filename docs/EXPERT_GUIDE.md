# ASPIRE Expert Guide

This guide covers execution, resource use, recovery and diagnosis. Start with the
[installation guide](installation.md) for installation and quick starts, the
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
./run_asv_pipeline.sh my_study.yml --rerun-from standard:FILTER_ASVS
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

V2 selects a successful run with submitted tasks for the requested execution phase.
Use `--resume-run` to select another retained run explicitly.

## Choosing the Earliest Rerun Stage

Start at the earliest process whose inputs, parameters, script, or environment
changed. Common examples:

| Change | Suggested stage |
|---|---|
| Primer definitions or screening | `core:PRIMER_TRIM` |
| FASTQ trimming or merge parameters | `core:FASTP_QC` or `core:MERGE_READS` |
| per-sample ASV relative-abundance threshold | `FILTER_ASVS` |
| taxonomy reference | `TAXONOMY` |
| host/non-target rules | `FILTER_ASVS` |
| plot annotations or displayed biological types | `PLOT_METADATA` |
| biological inclusion cutoff, sample classes, TECH/BIO thresholds | `CONTROL_DECONTAM` |
| Groups excluded from optional analyses | `optional:ANALYSIS_COHORT_DIAGNOSTICS` |
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
- Review control classifications and read loss when control decontamination
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
template, mock configuration and mock validator together so the test path
continues to exercise the advertised functionality.

## Selected analysis cohorts

`optional.analysis_cohort.exclude_groups` applies to grouping diagnostics and to the synchronized analysis tables selected after batch correction. Metadata plots, Sankey accounting and diversity retain those biological groups. Review the cohort sample-selection audits when changing this setting.

## V2 execution phases and cache management

## Restart, overwrite, or isolate an experiment

| Goal | Command |
|---|---|
| Resume ordinary work | `./run_asv_pipeline.sh my_run.yml` |
| Run sequence processing and canonical handoff | `./run_asv_pipeline.sh my_run.yml --phase preprocess` |
| Analyze an existing canonical dataset | `./run_asv_pipeline.sh my_run.yml --phase analysis` |
| Force mapping and dependent work | `./run_asv_pipeline.sh my_run.yml --rerun-from CREATE_COUNT_MATRIX` |
| Force one downstream analysis | `./run_asv_pipeline.sh my_run.yml --phase analysis --rerun-from DIVERSITY_ANALYSIS` |
| Execute without task-cache reuse | `./run_asv_pipeline.sh my_run.yml --no-resume` |
| Select a particular resume history | `./run_asv_pipeline.sh my_run.yml --resume-run RUN_NAME` |

`--rerun-from` forces the named task; true descendants rerun through dependency
invalidation. Independent branches can remain cached. It does not mean “only
run this task.” Required missing upstream work may also execute. Do not combine
`--no-resume` with `--rerun-from` or `--resume-run`.

The controller prefers the successful resume pointer for the selected phase.
An explicit `--resume-run` overrides it. In all-phase mode, absent a usable
pointer, the controller can select a managed run with the most completed/cached
tasks. Consult `cache list` rather than assuming the last terminal attempt is
always the baseline.

A rerun targeting the same output directory replaces published workflow outputs
on successful finalization. During execution, visible published results may
still describe the previous completed run. Current progress is in
`<runtime_dir>/publication_staging/logs/controller.log` and
`nextflow_trace.tsv`; task-specific `.command.sh`, `.command.out`, `.command.err`
and `.exitcode` live in the trace's hashed work directory.

For a side-by-side experiment, use a separate YAML, output directory and runtime
directory. Update every derived input path that points back into the old run;
changing only `paths.output_dir` is insufficient when a study recipe has explicit
cross-module paths. For analysis-only work, provide an intact canonical dataset
under the configured output directory rather than pointing at loose count files.

## Resources and Nextflow overrides

`sample_threads` is per task; `max_parallel_sample_tasks` is per sample-process
limit; `analysis_threads` controls dataset-level CPU requests. Concurrent
branches can have large combined memory demands. Work directories and published
intermediates can consume substantially more space than the compressed inputs.

Use an explicit Nextflow config for executor-specific settings, for example:

```groovy
// local_resources.config: adjust to your machine, not a required default.
executor {
    name = 'local'
    cpus = 8
    memory = '32 GB'
}
```

Set YAML analysis threads consistently with that budget, then launch:

```bash
./run_asv_pipeline.sh my_run.yml -- -c local_resources.config
```

Scheduler deployment additionally requires a site-appropriate executor, queue,
resource directives, shared input/work access, and Conda availability on compute
nodes. There is no universally configured cluster profile in this guide.

The wrapper pins its default Nextflow version through `NXF_VER` and uses the
v1 syntax parser. Override versions only as a separately validated change.
`ASPIRE_MAMBA_BUILD_TIMEOUT` controls the serialized package-build timeout
(default 900 seconds); it is not a timeout for biological analyses.

## Storage and cache lifecycle

```bash
./run_asv_pipeline.sh my_run.yml cache usage
./run_asv_pipeline.sh my_run.yml cache list --limit 20
./run_asv_pipeline.sh my_run.yml cache prune --keep 3
# Apply only after reviewing the dry run:
./run_asv_pipeline.sh my_run.yml cache prune --keep 3 --force
```

Keeping three runs when only three are eligible removes nothing. Protected phase
resume pointers may also prevent a requested retention count. Pruning operates
on managed run history; it does not empty all work directories, remove orphaned
work, delete publication staging, or clear Conda environments.

When you no longer need to resume this output's runs:

```bash
./run_asv_pipeline.sh my_run.yml cache clear --all --include-conda
# Apply the destructive cleanup after reviewing its targets:
./run_asv_pipeline.sh my_run.yml cache clear --all --include-conda --force
```

Clearing forfeits affected resume state and requires environments/tasks to be
rebuilt. Verify `paths.output_dir`, `runtime_dir`, and cache overrides first;
never clear a cache needed by an active run. Retain `.nextflow` history and
runtime state together when resume is required. `--no-resume` skips cache reuse
but is not a disk-cleaning command.

Use explicit cache commands in this revision. The launcher's
`keep_runtime_dir` YAML default expression can treat false as true; do not rely
on automatic deletion by setting false. Do not delete publication staging
without verifying which canonical/published inputs or custom links need it.
