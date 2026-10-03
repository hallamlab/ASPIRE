# Run your own study

Start from a working [installation](installation.md). A successful [reviewer test](reviewer-test.md) checks the installation before you introduce study-specific settings. See [input organization](inputs.md) before editing the YAML template.

## General Quick Start

Create a run config from the full template, then edit all paths for your environment:

```bash
cp asv_pipeline_nextflow.yml my_run.yml
```

At minimum, review:

- `core.paths.input_dir`
- `core.paths.output_dir`
- `core.paths.manifest`
- `core.paths.runtime_dir`
- `core.paths.keep_runtime_dir`
- `core.paths.work_dir`
- `core.paths.conda_cache_dir`
- any `/abs/path/...` placeholder
- enabled optional branches that require metadata, reference databases, VOC tables, or genome/MAG inputs

Run the pipeline:

```bash
./run_asv_pipeline.sh my_run.yml
```

List valid stage names for targeted reruns:

```bash
./run_asv_pipeline.sh --list-stages
```

Force a rerun from one stage onward using the retained default runtime cache:

```bash
./run_asv_pipeline.sh my_run.yml --rerun-from standard:PLOT_METADATA
```

Pass extra Nextflow options after `--`:

```bash
./run_asv_pipeline.sh my_run.yml -- -with-report report.html -with-trace trace.tsv
```

Direct Nextflow invocation is intended only for debugging because public-output
finalization is performed by the wrapper:

```bash
nextflow run asv_pipeline.nf --params-file my_run.yml --pipeline_config my_run.yml
```

**Review the outputs:** After a successful wrapper run, open
`<output_dir>/summary/report/ASPIRE_run_report.html` in a web browser. This is
the recommended starting point before biological interpretation. Review the
opening **Data Accounting Summary** for sample and group composition, sequence
and ASV retention, and the Sankey, swarmplot, UpSet, and collector's-curve
figures. Use the **Output Inventory** to open each enabled module's tables and
plots under `<output_dir>/modules/`. Use **Nextflow run details** to inspect task
status and resource use, execution timing, the workflow DAG, trace data,
software/runtime records, and controller or per-task logs under
`<output_dir>/logs/`. Machine-readable inventories, exact paths, and checksums
are stored under `<output_dir>/summary/tables/`.
