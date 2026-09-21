# ASPIRE expert guide

This guide complements the [README](../README.md) and
[YAML reference](CONFIGURATION.md). There is no separate expert-mode flag.
Use the supported wrapper for every production run. Examples below assume
that you are in the ASPIRE repository and have configured `my_run.yml`.

## Plan a run and preserve its identity

Record `git rev-parse HEAD`, the launch command, supplied YAML, sample manifest,
reference versions/checksums, and final output manifests. Conda environment YAMLs
can contain unpinned dependencies: a later environment solve may differ from an
existing cached environment. Retain the actual run's software records.

The full configuration template enables many analyses and contains example
paths. It is a configuration reference, not a minimal runnable dataset. Replace
active placeholders and turn off branches you cannot supply inputs for. The
mock configuration generator is the supported portable first-run example.
Local `examples/si.local.yml` and `examples/set1-2.local.yml` are study recipes,
not recommended universal settings.

Before a production run, confirm:

- Manifest IDs are unique and match the configured metadata ID column exactly.
- R1/R2 files refer to the correct paired sample; single-end assays use the
  separate single-end/full-length instructions.
- Primer removal and fixed trimming do not remove the same bases twice.
- Length and merging settings fit the assay, not merely another study's YAML.
- Every enabled module has its metadata, reference files, and upstream inputs.
- CPU, RAM, disk space, network access, and the required font are available.

```bash
./run_asv_pipeline.sh --help
./run_asv_pipeline.sh --list-stages
./run_asv_pipeline.sh my_run.yml -- -preview
```

Preview assembles/validates workflow construction without executing analysis
tasks. It is not proof that every input schema, tool invocation, API or environment
will work. Controller initialization may still occur.

## Parameter semantics that matter

| Setting | Meaning and review point |
|---|---|
| `fastp.trim_front_r1` and related fields | Fixed numbers of bases removed, not primer identification. The standalone [primer tool](../scripts/PRIMER_TRIMMING.md) is optional and outside Nextflow. |
| `merge.min_overlap` | Minimum aligned overlap in bases. Raising it may remove valid shortened pairs. |
| `merge.max_diffs` | Maximum disagreements in the overlap, not the final merged-read error count. |
| `filter.max_ee` | Maximum expected errors calculated from read qualities; distinct from observed overlap disagreements. |
| `filter.min_len`, `max_len` | Read-length bounds; inspect length distributions and taxon-specific losses. |
| `unoise.min_size` | Pooled abundance of a unique input sequence required for denoising. It is not per-sample ASV prevalence. |
| `read_mapping.min_identity` | Read/reference alignment identity fraction. Current count construction uses first-qualifying-hit assignment. |
| `table_filter.min_sample_sum` | Minimum mapped counts in a sample. |
| `table_filter.min_asv_sum` | Percent relative abundance attained in at least one retained sample. Zero is honored. |
| `filter_counts.abundance_threshold` | Later percent-abundance threshold, after non-target processing. |
| `spieceasi.min_rel_abund` | Fractional maximum relative abundance across samples: 0.005 means 0.5% in at least one sample. |
| `spieceasi.min_prevalence` | Fraction of samples with positive counts: 0.005 means 0.5% of samples. |

Force-retained ASV lists can bypass network abundance/prevalence gates;
zero-variance removal still applies. Inspect the filtering audit rather than
assuming all retained ASVs passed every standard cutoff.

### Read mapping and ambiguity

```yaml
read_mapping:
  min_identity: 0.999
```

The workflow fallback is 0.999. The SI recipe currently uses 0.995. For a
253-bp ungapped A/C/G/T alignment, these allow zero and one substitutions.
An identity fraction is not a general one-error cap: alignment length, gaps,
ambiguous bases, and the identity definition matter.

VSEARCH orders candidates by shared k-mers, aligns them, and by default stops
at the first qualifying hit (`maxaccepts=1`, `maxrejects=32`). The current
workflow does not expose search-exhaustiveness or ambiguity policy in YAML.
Adding unknown YAML fields does not enable those behaviors.

VSEARCH itself supports `--maxaccepts 0 --maxrejects 0` to consider all candidates
surviving its prefilter, and `--top_hits_only` to report highest-identity hits.
These are not equivalent to rejecting ambiguous assignments. In VSEARCH 2.31.0,
`--otutabout` assigns a read to one top-ranked accepted reference; equal-identity
ties favor earlier database entries. An all-best-hits report and a count table
therefore have different semantics. Unique-best-only counting would require a
separate implementation and validation. See the [VSEARCH manual](https://torognes.github.io/vsearch/commands/vsearch-usearch_global.1.html)
and [versioned count-table implementation](https://github.com/torognes/vsearch/blob/v2.31.0/src/search.cc).

Report mapping recovery against the actual filtered-read denominator. Denoising,
chimera removal, and mapping are separate filters: fewer candidate ASVs can
coexist with stable mapped abundance. A higher mapping rate is not proof of
more accurate abundance estimates. Match ASVs between runs by sequence, not
numbered IDs.

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

## Statistical and ecological review

- Verify the actual sample cohort for each test. Missing measurements and rare
  groups can change it; overlays may intentionally include additional samples.
- Set blocking/pairing to the sampling design. Samples from one cruise or
  participant are not automatically independent replicates.
- Interpret corrected p-values together with effect sizes and test scope.
  An omnibus association does not identify a significant specific contrast.
- A fitted compartment association is not a predictive advantage. Use paired
  held-out comparisons and check dispersion diagnostics where provided.
- Batch correction and outlier removal require their own audits. Check for
  confounding between technical batches and biological strata.
- Inspect SPIEC-EASI stability and the selected lambda before extending its
  search range. Compare module membership, not just module numbers, across runs.
- Prefer prespecified sensitivity comparisons and report changes transparently;
  do not select parameter settings to obtain significance.

For generated BASIN/BASINS-linked diversity cohorts, use
`diversity.matched_cohort.source: community_predictor_comparison` and enable
that producer. This is an explicit task dependency. A pre-existing external
cohort instead uses `matched_cohort.table`; do not specify both. Statistical
cohort restriction does not globally discard ASV samples from unrelated work.

## Failure handling and useful diagnostics

| Symptom | Check / next action |
|---|---|
| `[MITOMASTER STOPPED]` | Read `mitomaster_output.tsv.status.json` in the task work directory. Check endpoint connectivity/status, HTTP error, timeout or response format. Partial results are quarantined, not accepted. |
| Missing required font | Install the requested font and confirm `fc-match` returns it rather than a substitute. |
| Missing generated cohort | Enable its producer and use `matched_cohort.source`; do not point at a not-yet-published file. |
| Very low mapping recovery | Check primers, sequence boundaries, alignment settings, chimera/denoising losses and count-table denominators. |
| Apparently unchanged published results | Inspect current runtime logs; final publication happens after successful completion. |
| Cache prune freed nothing | Compare retained versus eligible run records and protected pointers; inspect the separate staging/Conda categories. |
| No new process after editing a helper | Use the appropriate `--rerun-from` if that script is not represented in task hashing. |
| Environment build failed | Review network/repodata/solver errors. Waiting for the serialized build slot alone is not a failure. |

MITOMASTER HTTP/network failures retry according to the configured limits, then
stop. Valid header-only results mean no reported hits; an empty response cannot
safely be interpreted that way. If you intentionally omit the API, set
`mito.run_mitomaster: false` and document that remote evidence was not used.
Resuming after a failed task is preferable to removing all cached work.

## Extending and testing the workflow

Process definitions live in `workflow/modules/`; configuration, channels and
subworkflow calls live in `asv_pipeline.nf`. The wrapper assembles them into a
same-directory generated launch script. Add new processes to the relevant
module, wire explicit dependencies, document YAML keys, provide an environment,
and update the controller stage registry and output organizer as needed.

Keep sample identifiers, canonical-dataset schemas, output inventories and
hash-based provenance consistent. Use shared figure styling and avoid hidden
inputs that bypass task dependencies.

The tests combine pytest functions, unittest cases, external-tool fixtures and
Nextflow integration fixtures. Use an appropriate development environment with
pytest and the dependencies required by the selected tests:

```bash
python -m pytest tests/test_cache_management.py tests/test_phase_architecture.py
python -m unittest discover -s tests -p 'test_mitomaster_output_contract.py'
python -m unittest discover -s tests -p 'test_trim_amplicon_primers.py'
```

Some integration tests skip without Nextflow, Cutadapt or analysis dependencies.
A skipped test is not successful tool validation. The full biological smoke test
is the mock workflow plus its supplied validator; syntax tests do not replace it.
