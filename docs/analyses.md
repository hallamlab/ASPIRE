# Choose and interpret analyses

The optional branches use the selected downstream tables. Read [decontamination](decontamination.md) before enabling control-based filtering and [VOC statistics](VOC_STATISTICS.md) before interpreting repeated-participant associations. Use the [process reference](PROCESS_REFERENCE.md) to trace individual inputs and outputs.

## Indicator Species Analysis

`INDICSPECIES` runs the primary analyses listed in `optional.indicspecies.group_cols` when enabled, preserving the standard outputs such as `Type_Group_indicator_species_summary.tsv` and `Case_indicator_species_summary.tsv`.

Additional nested ISA runs can be requested with `optional.indicspecies.stratified`. Each analysis tests `group_col` separately within each selected `within_col` value. For example, the VOC-enabled example config tests cancer/control indicators within each respiratory sample type:

```yaml
optional:
  indicspecies:
    stratified:
      enabled: true
      analyses:
        - within_col: Type_Group
          group_col: Case
          levels:
            - Oral Rinse
            - BAL
            - Bronchial Brush
```

These stratified runs write per-level and pooled tables named like `stratified_Case_within_Type_Group_Bronchial_Brush_indicator_species_summary.tsv` and `stratified_Case_within_Type_Group_indicator_species_summary.tsv`.


## Count Filtering

`FILTER_COUNTS` produces the ASV count tables used by downstream analyses. The important outputs are:

- `ASV_target.tsv`: final microbial ASV table after contaminant removal, mitochondrial removal, abundance/prevalence filtering, taxonomy-quality filtering, and explicit taxon exclusions.
- `ASV_target.micro.tsv`: intermediate microbial table before final abundance/taxonomy filtering; kept for audit and data-loss summaries.
- `ASV_target.decon.tsv`: intermediate decontaminated table.
- `ASV_target.mito.tsv`: mitochondrial table.

Downstream metadata and VOC analyses use `ASV_target.tsv`, not the intermediate `.micro.tsv`, so taxa excluded by final filtering should not re-enter later outputs.

Explicit taxon exclusions are configured with `standard.filter_counts.exclude_taxa`. Entries are exact, case-insensitive matches against parsed taxonomy ranks, and underscores in configured values are normalized to spaces. Supported forms include strings, maps, and mixed lists:

```yaml
standard:
  filter_counts:
    enabled: true
    exclude_taxa:
      - "Species:Homo sapiens"
      - "Class:Mammalia"
      - Order: Primates
```

Use this for host or other known non-target ranks that should be removed even if they pass sequence and abundance filters.


## Optional Three-Tier Decontamination

Set `optional.three_tier_decontam.enabled: true` to run the SPARK manuscript contamination screen after `PLOT_METADATA` and before batch correction, diversity, indicator-species, SPIEC-EASI/network, VOC, clustermaps, and the other downstream modules. The statistical model uses the raw pre-filter ASV count matrix so negative controls remain available, while its decisions are applied to the final host-filtered long and wide microbial tables.

```yaml
optional:
  three_tier_decontam:
    enabled: true
    metadata: /path/to/sample_metadata.tsv
    output_dir: three_tier_decontam
    metadata_sample_col: Sample
    sample_types: [BAL, Bronchial Brush, Oral Rinse]
    pooled_threshold: 0.1
    within_type_threshold: 0.1
    aggressive_threshold: 0.5
    combine_mode: min
    biological_plausibility: true
```

The complete score tables, flags, taxonomy audit, read-loss summary, and filtered long/wide tables are published under the configured `output_dir`. The module is disabled by default. For `RUN_FROM_FINAL_CHECKPOINT`, set `checkpoint.raw_asv_counts` to the original control-bearing count matrix; falling back to `checkpoint.asv_counts` is only valid when that table still contains the negative controls.


## Metadata And ASV Outputs

`PLOT_METADATA` builds the run's metadata-linked ASV products. Typical outputs include:

- `modules/metadata_plots/tables/ASV_meta_micro.tsv`
- `modules/metadata_plots/tables/ASV_final.micro.tsv`
- `modules/metadata_plots/tables/metadata_updated_micro.tsv`
- `modules/metadata_plots/tables/master_table_micro.tsv`
- run metadata summaries and plots

When batch correction is enabled, corrected count and metadata tables are produced and downstream branches that support corrected inputs use them.


## VOC Correlation

`VOC_CORRELATION` links filtered ASV abundances to measurements of Volatile Organic Compounds (VOCs). It requires `optional.voc_correlation.enabled: true`, a metadata table, the final filtered ASV counts, and `optional.voc_correlation.voc_table`.

Direction filtering is controlled by:

```yaml
optional:
  voc_correlation:
    enabled: true
    correlation_direction: positive  # positive, negative, or both
```

This setting applies to ASV-VOC correlation tables and ASV-VOC correlation heatmaps. For example, `positive` keeps only positive ASV-VOC correlations in the reported long table and correlation clustermap, allowing statements such as "VOC abundance was positively correlated with ASV X." Multiple-testing q-values are computed across the tested ASV-VOC pairs before direction filtering.

VOC abundance plots are different from correlation plots:

- `asv_voc_clustermap*` shows ASV-VOC correlation values and respects `correlation_direction`.
- `sample_voc_brush_clustermap*` shows per-sample VOC abundance z-scores, not correlations. Blue indicates lower-than-average VOC abundance for that VOC, white is near the VOC mean, and orange indicates higher-than-average abundance.
- `patient_case_voc_barplots_brush*` displays per-VOC patient z-scores so VOCs are visually comparable on one axis; statistical tests are still run on original patient-level VOC values.


## Standard and Optional Branches

The standard final-table sections are `standard.mito`,
`standard.filter_counts`, `standard.general_stats`, and
`standard.metadata_plots`. They perform non-target screening, final count
filtering, run summaries, and metadata-linked final-table construction.

Selectable analyses live under `optional` and are controlled by their
`enabled` flags:

- `optional.three_tier_decontam`: control-based prevalence, frequency, and biological-plausibility decontamination.
- `optional.plot_upset`, `optional.bubbleplotter`, `optional.umap_clustering`: metadata visualization branches.
- `optional.batch_correction` and `optional.outlier_detection`: corrected ASV tables and outlier checks.
- `optional.collectors_curve`: rarefaction/collector curve summaries.
- `optional.diversity`: Shannon, Bray-Curtis, Jaccard, and patient-aware diversity workflows.
- `optional.indicspecies`: indicator species tables, standard plots, and aligned indicator plots.
- `optional.voc_correlation`: VOC-ASV association analysis and VOC abundance visualizations.
- `optional.clustermaps`: ASV and metadata heatmaps.
- `optional.spieceasi`: SPIEC-EASI inference, module detection, network visualization, and optional module/MAG anchor tables.
- `optional.network_topology`: SPARK-compatible global topology statistics and degree-preserving null-network comparison.
- `optional.power_analysis`: patient-aware power analysis using metadata-linked ASV tables.
- `optional.taxonomy_patient_aware`: patient-aware taxonomic comparisons.
- `optional.lung_status_analysis`: patient-aware lung-status comparisons.
- `optional.asv_mag_link`: ASV-to-MAG/barrnap linkage and module anchoring.
- `optional.sankey`: data-loss and filtering Sankey summaries.
- `optional.master_summary`: final combined ASV summary export.


## SPARK Supplement Functionality

The mapping from analyses reported in the SPARK and ASPIRE
manuscripts to workflow switches, implementations, and outputs is documented in
[`docs/SPARK_FUNCTIONALITY.md`](SPARK_FUNCTIONALITY.md). The private
manuscript-scale run uses 1,000 degree-preserving null networks. The public mock
configuration uses 100 draws to demonstrate the same method in less time.
