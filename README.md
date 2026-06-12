# ASPIRE

ASPIRE is a Nextflow DSL2 workflow for ASV generation, taxonomy assignment, decontamination, metadata-linked ASV summaries, ecological analyses, VOC association analyses, network/module analyses, and optional ASV-to-MAG linkage.

The supported entrypoint is `run_asv_pipeline.sh`. It bootstraps the controller environment, launches `asv_pipeline.nf`, manages resume behavior, and supports stage-aware reruns with `--rerun-from`.

The canonical user documentation is this README. `ASPIRE.ipynb` is kept as a short run notebook; the old duplicate technical notebook has been removed.

## Key Files

- `run_asv_pipeline.sh`: main wrapper for routine runs.
- `asv_pipeline.nf`: current Nextflow workflow.
- `asv_pipeline_nextflow.yml`: full config template.
- `examples/set1-2.local.yml`: local example config with absolute paths for the UBC/LMP test dataset.
- `processes/`: scripts and conda environment YAMLs used by individual stages.

## Quick Start

Install runtime prerequisites:

```bash
command -v mamba
```

Create a run config from the full template, then edit all paths for your environment:

```bash
cp asv_pipeline_nextflow.yml my_run.yml
```

At minimum, review:

- `paths.input_dir`
- `paths.output_dir`
- `paths.manifest`
- `paths.work_dir`
- `paths.conda_cache_dir`
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

Force a rerun from one stage onward while preserving cacheability for future resumes:

```bash
./run_asv_pipeline.sh my_run.yml --rerun-from PLOT_METADATA
```

Pass extra Nextflow options after `--`:

```bash
./run_asv_pipeline.sh my_run.yml -- -with-report report.html -with-trace trace.tsv
```

Direct Nextflow invocation is supported, but the wrapper is preferred:

```bash
nextflow run asv_pipeline.nf --params-file my_run.yml --pipeline_config my_run.yml
```

## Inputs

FASTQs can be discovered from `paths.input_dir`, but manifest mode is preferred for reproducible sample names.

Manifest format:

- Tab-separated, no header.
- Column 1: `sample_id`.
- Column 2: R1 FASTQ.
- Column 3: R2 FASTQ, optional for single-end data.
- Lines starting with `#` are ignored.
- Relative FASTQ paths are resolved relative to the manifest file.

Metadata is required by enabled metadata-aware branches such as metadata plots, Sankey, diversity, indicator species, VOC correlation, power analysis, taxonomy patient-aware analysis, lung-status analysis, and several network overlays. The configured sample column must match the manifest sample IDs.

Reference inputs depend on enabled branches:

- `sina.reference` and taxonomy references are used for SINA alignment and taxonomy assignment. The config can point at local files or URLs.
- `mito.mito_db` and `mito.biof_db` are required when mitochondrial/contaminant decontamination is enabled.
- `voc_correlation.voc_table` is required when VOC correlation is enabled.
- `asv_mag_link.*` inputs are required only when ASV-to-MAG linkage is enabled.

## Workflow Stages

The wrapper's current stage order is:

1. `FASTP_QC`
2. `MERGE_READS`
3. `FILTER_READS`
4. `RELABEL_FILTERED`
5. `CONCAT_FASTAS`
6. `DEREPLICATE`
7. `DENOISE`
8. `CHIMERA_CHECK`
9. `CREATE_COUNT_MATRIX`
10. `FILTER_TABLE`
11. `SINA_TRIM`
12. `TAXONOMY`
13. `MITOMASTER`
14. `MITO_DECONTAM`
15. `FILTER_COUNTS`
16. `GENERAL_STATS`
17. `PLOT_METADATA`
18. `PLOT_UPSET`
19. `ASV_BATCH_CORRECTION`
20. `ASV_META_FROM_CORRECTED`
21. `BUBBLEPLOTTER`
22. `UMAP_CLUSTERING`
23. `OUTLIER_CHECKER`
24. `COLLECTORS_CURVE`
25. `DIVERSITY_ANALYSIS`
26. `INDICSPECIES`
27. `INDICSPECIES_PLOTS`
28. `VOC_CORRELATION`
29. `CLUSTERMAPS`
30. `POWER_ANALYSIS_PIPELINE`
31. `SPIECEASI`
32. `NETWORK_MODULES`
33. `ASV_MAG_LINK`
34. `GRAPH_NETWORK`
35. `MODULE_MAG_ANCHORS`
36. `SANKEY`
37. `MASTER_SUMMARY`

Disabled optional branches are skipped based on the YAML config.

## Count Filtering

`FILTER_COUNTS` produces the ASV count tables used by downstream analyses. The important outputs are:

- `ASV_target.tsv`: final microbial ASV table after contaminant removal, mitochondrial removal, abundance/prevalence filtering, taxonomy-quality filtering, and explicit taxon exclusions.
- `ASV_target.micro.tsv`: intermediate microbial table before final abundance/taxonomy filtering; kept for audit and data-loss summaries.
- `ASV_target.decon.tsv`: intermediate decontaminated table.
- `ASV_target.mito.tsv`: mitochondrial table.

Downstream metadata and VOC analyses use `ASV_target.tsv`, not the intermediate `.micro.tsv`, so taxa excluded by final filtering should not re-enter later outputs.

Explicit taxon exclusions are configured with `filter_counts.exclude_taxa`. Entries are exact, case-insensitive matches against parsed taxonomy ranks, and underscores in configured values are normalized to spaces. Supported forms include strings, maps, and mixed lists:

```yaml
filter_counts:
  enabled: true
  exclude_taxa:
    - "Species:Homo sapiens"
    - "Class:Mammalia"
    - Order: Primates
```

Use this for host or other known non-target ranks that should be removed even if they pass sequence and abundance filters.

## Metadata And ASV Outputs

`PLOT_METADATA` builds the run's metadata-linked ASV products. Typical outputs include:

- `metadata/ASV_meta.tsv`
- `metadata/ASV_meta_micro.tsv`
- `metadata/ASV_final.tsv`
- `metadata/ASV_final.micro.tsv`
- run metadata summaries and plots

When batch correction is enabled, corrected count and metadata tables are produced and downstream branches that support corrected inputs use them.

## VOC Correlation

`VOC_CORRELATION` links VOC abundances to filtered ASV abundances. It requires `voc_correlation.enabled: true`, a metadata table, the final filtered ASV counts, and `voc_correlation.voc_table`.

Direction filtering is controlled by:

```yaml
voc_correlation:
  enabled: true
  correlation_direction: positive  # positive, negative, or both
```

This setting applies to ASV-VOC correlation tables and ASV-VOC correlation heatmaps. For example, `positive` keeps only positive ASV-VOC correlations in the reported long table and correlation clustermap, allowing statements such as "VOC abundance was positively correlated with ASV X." Multiple-testing q-values are computed across the tested ASV-VOC pairs before direction filtering.

VOC abundance plots are different from correlation plots:

- `asv_voc_clustermap*` shows ASV-VOC correlation values and respects `correlation_direction`.
- `sample_voc_brush_clustermap*` shows per-sample VOC abundance z-scores, not correlations. Blue indicates lower-than-average VOC abundance for that VOC, white is near the VOC mean, and orange indicates higher-than-average abundance.
- `patient_case_voc_barplots_brush*` displays per-VOC patient z-scores so VOCs are visually comparable on one axis; statistical tests are still run on original patient-level VOC values.

## Optional Analysis Branches

Major optional modules are controlled by YAML `enabled` flags:

- `mito`: BLAST-based mitochondrial and contaminant screening.
- `filter_counts`: count filtering, intermediate count audit tables, and explicit taxon exclusions.
- `general_stats`: run-level ASV and sample summaries.
- `metadata_plots`: metadata-linked ASV summary tables and plots.
- `plot_upset`, `bubbleplotter`, `umap_clustering`: metadata visualization branches.
- `batch_correction` and `outlier_detection`: corrected ASV tables and outlier checks.
- `collectors_curve`: rarefaction/collector curve summaries.
- `diversity`: Shannon, Bray-Curtis, Jaccard, and optional patient-aware diversity workflows.
- `indicspecies`: indicator species analysis and aligned indicator plots.
- `voc_correlation`: VOC-ASV association analysis and VOC abundance visualizations.
- `clustermaps`: ASV and metadata heatmaps.
- `spieceasi`, `network_modules`, `graph_network`: SPIEC-EASI network inference, module detection, and network visualization.
- `power_analysis`: patient-aware power analysis using metadata-linked ASV tables.
- `taxonomy_patient_aware` and `lung_status_analysis`: patient-aware taxonomic and lung-status comparisons.
- `asv_mag_link` and `module_mag_anchors`: ASV-to-MAG/barrnap linkage and module anchoring.
- `sankey`: data-loss and filtering Sankey summaries.
- `master_summary`: final combined ASV summary export.

## Output Structure

Typical top-level output directories include:

- `fastp/`, `merged/`, `filtered/`, `concat/`, `derep/`, `denoise/`, `nochimeras/`
- `ASVs/`, `sina/`, `taxonomy/`, `mito/`, `stats/`, `logs/`
- `metadata/`, `batch_correction/`, `outliers_corrected/`
- `diversity/`, `indicspecies/`, `voc_correlation/`, `clustermaps/`
- `spieceasi/`, `network/`, `asv_mag_link/`
- `power_analysis/`, `taxonomy_patient_aware/`, `lung_status_analysis/`
- `sankey/`, `master_summary/`

The actual directory set depends on which branches are enabled.

## Runtime And Cache Behavior

The wrapper reads `paths.work_dir` and `paths.conda_cache_dir` from the YAML and exports `NXF_WORK` and `NXF_CONDA_CACHEDIR` when configured. It also creates the controller environment at `.controller_env`.

Normal repeated wrapper runs resume from Nextflow cache. If a branch does not rerun because cached outputs are valid, use `--rerun-from STAGE_NAME`.

Several stages include checksums of external process scripts in their task commands, so edits to important Python/R helper scripts invalidate the relevant Nextflow task cache. This avoids stale outputs when a script changes but the input filenames stay the same.

## Troubleshooting

- `No usable entries detected in manifest`: check tab separation, sample IDs, and FASTQ paths.
- `metadata file not found`: set the branch-specific metadata path or disable that branch.
- Host taxa still appear downstream: confirm `filter_counts.exclude_taxa` is set and rerun from `FILTER_COUNTS` or at least from `PLOT_METADATA` if the final `ASV_target.tsv` is already corrected.
- VOC direction did not change outputs: rerun from `VOC_CORRELATION`.
- Sankey complains about intermediates: set `filter_counts.save_intermediates: true`.
- BLAST database errors: set `mito.mito_db` and `mito.biof_db` to valid database prefixes or compatible FASTA paths.
- Conda solve errors: confirm `mamba` is available and review the relevant `environments.*` config entry.
