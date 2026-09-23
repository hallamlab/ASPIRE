# ASPIRE Configuration Reference

This guide describes the configuration accepted by the `main` branch. The
authoritative complete template is
[`asv_pipeline_nextflow.yml`](../asv_pipeline_nextflow.yml). Its exhaustive
[parameter catalogue](CONFIG_PARAMETERS.md) defines every key and template
value. Copy the template for a production run and retain the resolved copy with
the results. Private study YAMLs are intentionally excluded from the repository
because they contain machine-specific paths and study provenance.

## How Configuration Is Applied

Run ASPIRE through the supported wrapper:

```bash
./run_asv_pipeline.sh my_run.yml
```

The wrapper supplies the same resolved YAML to the controller and Nextflow.
The complete template and generated mock configuration organize process
sections into three top-level namespaces:

- `core`: always-scheduled ASV construction and taxonomy settings.
- `standard`: canonical final-table preparation required by the usual
  downstream analysis graph.
- `optional`: independently selectable analytical and reporting branches.

`environments` remains top-level because it is an operational software mapping
shared by all tiers. A branch with `enabled: false` is not scheduled. An
optional branch may intentionally be `enabled: true` in the mock benchmark;
enabled-by-default does not mean mandatory. Enabling a branch does not create
missing biological inputs: its metadata columns, reference files, control
labels, VOC table, or genome inputs must also be supplied.

Legacy flat YAML configurations remain supported. A section must appear either
flat or inside one tier, never in both places.

## Module Dependency Guide

| Section | Required upstream data | Additional required input | Main published module |
|---|---|---|---|
| `standard.mito` | filtered ASVs and taxonomy | BLAST databases/FASTAs; MITOMASTER network access when enabled | `non_target_filtering` |
| `standard.filter_counts` | taxonomy and non-target evidence | metadata | `non_target_filtering` |
| `optional.three_tier_decontam` | raw count matrix and final microbial tables | control-bearing metadata and DNA concentration | `contamination_filtering` |
| `standard.metadata_plots` | final microbial/mitochondrial tables | metadata | `metadata_plots` |
| `optional.batch_correction` | metadata-linked ASV tables | batch and biological-covariate columns | `batch_correction` |
| `optional.outlier_detection` | metadata-linked ASV table | configured grouping columns | `outlier_detection` |
| `optional.diversity` | downstream ASV counts | metadata/group columns | `diversity` |
| `optional.indicspecies` | downstream ASV counts | grouping columns; optional patient block | `indicator_analysis` |
| `optional.voc_correlation` | downstream ASV tables and optional ISA results | VOC table and sample/patient matching columns | `voc_correlation` |
| `optional.power_analysis` | downstream ASV tables and ISA completion | patient, case and type columns | `power_analysis` |
| `optional.taxonomy_patient_aware` | downstream ASV tables | patient/case/type metadata | `taxonomy` |
| `optional.lung_status_analysis` | downstream ASV tables | patient, cancer-site and lung-side metadata | `lung_status_analysis` |
| `optional.spieceasi` | downstream wide ASV table | none beyond configured filters | `network_analysis` |
| `optional.network_topology` | enabled SPIEC-EASI thresholded graph | none | `network_analysis` |
| `optional.asv_mag_link` | filtered ASV sequences | genome-QC, barrnap, and/or genome FASTA inputs | `asv_mag_link` |
| `optional.master_summary` | metadata-linked ASV tables | enabled branch outputs are incorporated when available | `summary` |

When three-tier decontamination is enabled, its filtered long and wide tables
replace the corresponding metadata-stage tables for downstream analyses. When
batch correction is enabled, its corrected count table replaces those counts
for the downstream branches that support correction. These are analytical
choices, not merely extra plots.

## Core Input and Runtime

### `core.paths`

- `input_dir`: directory searched for FASTQs when no manifest is supplied.
- `output_dir`: public result directory. Use a new directory for an independent
  analysis.
- `manifest`: optional TSV mapping `sample_id`, R1 and optional R2 FASTQs.
- `runtime_dir`: persistent Nextflow work, Conda cache and staging root;
  defaults to `<output_dir>/.aspire`.
- `keep_runtime_dir`: retain runtime state after success. Keep `true` for
  resume and targeted reruns.
- `work_dir`, `conda_cache_dir`: optional granular runtime overrides.

### `core.resources`

- `threads`: CPU request used by per-sample read-processing tasks.
- `single_end`: set `true` only for single-end data; paired-end is the default.

### `core.filename_patterns`

`r1_tokens`, `r2_tokens`, `ext_patterns`, and `sample_strip_regex` control FASTQ
discovery and sample-ID normalization. Prefer an explicit manifest when names
do not follow a consistent convention. Confirm the normalized
`summary/tables/run_manifest.tsv` before interpreting results.

## Read Processing and ASV Generation

| Section | Important fields | Effect |
|---|---|---|
| `core.fastp` | `trim_front_r1`, `trim_tail_r1`, `trim_front_r2`, `trim_tail_r2` | Fixed end trimming before merging. Values are read-specific base counts. |
| `core.merge` | `max_diffs`, `min_overlap`, `trunc_quality`, `allow_stagger` | Paired-read merge requirements. |
| `core.filter` | `max_ee`, `min_len`, `max_len` | Expected-error and merged-length filters. |
| `core.concat` | `relabel`, `label_sep` | Preserves sample identity when filtered FASTAs are concatenated. |
| `core.unoise` | `min_size` | Minimum dereplicated abundance supplied to denoising. |
| `core.swarm` | `distance` | ASV clustering/mapping distance setting. |
| `core.table_filter` | `min_sample_sum`, `min_asv_sum` | Early technical count filter. `min_sample_sum` is a read count, not a percentage. |

Changing any of these settings changes the inferred ASV cohort and normally
requires a fresh or appropriately targeted upstream run.

## Alignment and Taxonomy

### `core.sina`

Set either `reference` to a local SINA ARB reference or `reference_url` to a
downloadable archive. `regions` lists accepted variable-region annotations and
`trim_to` selects the desired region. `batch_size`, `threads`, `keep_gaps`, and
`verbose` control execution rather than biological filtering.

### `core.taxonomy`

`ref_taxonomy` and `ref_sequences` point to local QIIME2 artifacts; their URL
counterparts permit retrieval when local files are absent. Output names and
`download_subdir` control organization. Pin local reference artifacts when an
exact historical classification must be reproduced.

## Non-Target and Host Filtering

### `standard.mito`

- `enabled`: run mitochondrial/contaminant evidence generation.
- `run_mitomaster`: contact the external MITOMASTER service. Set `false` for an
  intentionally offline run; local taxonomy and BLAST evidence still run.
- `chunk_size`, `mitomaster_workers`, `mitomaster_retries`,
  `mitomaster_timeout`, `mitomaster_header_mode`: API chunking and retry policy.
- `min_pident`, `min_percov`: local BLAST identity and coverage cutoffs.
- `mito_db`, `biof_db`: existing nucleotide database prefixes.
- `mito_fasta`, `contaminant_fasta`: FASTA alternatives that take precedence
  and are rebuilt into run-specific databases.
- `steps`, `host_first_step`, `mitochondria_substring`: evidence labels and
  reporting order.

### `standard.filter_counts`

- `metadata`, `sample_id_col`, `group_col`: sample/group mapping.
- `min_group_size`: minimum group size considered by group-aware filtering.
- `abundance_threshold`: percentage-scale relative-abundance threshold used by
  this process; inspect the recorded configuration when comparing runs.
- `min_consensus`: minimum taxonomy consensus value.
- `exclude_taxa`: explicit exact rank/value exclusions such as
  `Species:Homo sapiens` or `Class:Mammalia`.
- `taxon_col`, `consensus_col`, `biofactorial_col`, `mito_cols`: input evidence
  column names.
- `save_intermediates`: retain `.decon`, `.micro`, and mitochondrial audit
  tables required by data-loss reporting.

The final downstream microbial table is `ASV_target.tsv`. The `.micro.tsv` and
`.decon.tsv` files are audit intermediates, not authoritative replacements.

### `optional.three_tier_decontam`

This optional branch requires extraction controls in the raw count matrix.
`metadata` must include the configured sample identifier, negative/positive
control flags, DNA concentration, and sample type.

- `negative_control_col`, `positive_control_col`: Boolean/control-label fields.
- `negative_control_labels`, `positive_control_labels`: accepted labels when
  controls are encoded categorically.
- `concentration_col`: DNA concentration for frequency modeling.
- `type_col`, `sample_types`: biological strata used by within-type models.
- `pooled_threshold`, `within_type_threshold`, `aggressive_threshold`:
  contaminant score thresholds.
- `combine_mode`: rule used to combine evidence tiers.
- `biological_plausibility`: apply the final plausibility screen.

The raw control-bearing table is used for scoring; decisions are applied to the
host-filtered microbial tables before all downstream analyses.

## Metadata, QC and Visualization

### `standard.metadata_plots`

`metadata` and `sample_col` establish sample matching. `type_col`, `color_col`,
`palette_file`, `group_order`, and `include_rank` control annotations.
`subtraction_group_col` and `subtraction_groups` identify groups subtracted or
removed before the biological table is finalized. `keep_types` restricts
retained biological types. `input_table` chooses the count-filter checkpoint;
`min_pre_correction_sample_sum` is an additional sample read-count threshold.
`drop_zero_asvs` removes features that become all zero. `run_micro` and
`run_mito` select table families.

### Other descriptive branches

- `standard.general_stats`: enables run-level FASTQ and sequence summaries.
- `optional.sankey`: maps sample/group columns and ordering into data-loss diagrams.
- `optional.plot_upset`: configures domain, groups, Venn/UpSet behavior and formats.
- `optional.bubbleplotter`: controls count/group columns, output formats and sizing.
- `optional.umap_clustering`: controls normalization, UMAP and HDBSCAN settings.
- `optional.collectors_curve`: controls grouping, permutations, seed and presence rule.
- `optional.clustermaps`: selects ranks, top-N limits, ISA threshold columns, palettes,
  sizes, formats and mitochondrial-table behavior.

## Batch Correction and Outliers

`optional.batch_correction.enabled` activates ConQuR-based correction using `batch_col`
and `biological_covariates`. Its `conqur_*` fields configure the correction;
the UMAP/HDBSCAN fields configure diagnostics. Enabling it changes the count
table used by downstream analyses and should be justified in the study design.

`optional.outlier_detection` configures Isolation Forest, one-class SVM and HDBSCAN.
`group_cols` defines the strata, `vote_threshold` controls consensus, and the
`iso_*`, `svm_*`, and `hdbscan_*` fields tune individual detectors. Outlier
results are diagnostic unless a separately documented process removes samples.

## Ecological and Association Analyses

### `optional.diversity`

`sample_col`, `group1_col`, `group2_col`, `color_col`, `group_order`, and
`exclude_groups` define cohorts. `block_col` supplies the participant/block for
repeated measures. `permanova_perms` and `random_state` control inference;
`umap_neighbors` and `umap_min_dist` control visualization. The nested
`patient_aware` block configures patient/case/type columns, contralateral rules,
permutations, transformation and complete-type requirements.

### `optional.indicspecies`

`group_cols` lists primary groupings and `block_col` optionally restricts
permutations. `perms`, `seed`, `q_threshold`, and `min_n` control tests.
`stratified.analyses` requests a `group_col` analysis separately within selected
levels of another column. Palette/order fields control display only.
`plot_enabled` and `aligned_plot_enabled` independently schedule the two plot
families; aligned plot thresholds are configured with `aligned_*` fields.

### `optional.voc_correlation`

- `voc_table`, `voc_sample_col`: VOC matrix and its identifier.
- `metadata_sample_col`, `patient_col`, `case_col`, `type_col`: matching and
  grouping fields.
- `sample_types`: eligible sample types.
- `sample_id_mode`, `use_legacy_voc_subset`: matching/subsetting behavior for
  established datasets.
- `correlation_direction`: reported direction for the baseline ASV-VOC output.
- `isa_correlation_direction`: direction for ISA-focused correlation outputs.
- `isa_brush_groups`, `isa_all_type_groups`,
  `isa_exclude_all_types_from_brush`: membership rules for ISA/VOC figures.
- `isa_min_abs_rho`: minimum absolute Spearman magnitude for focused rows and
  columns.
- `sample_min_abs_z`: minimum absolute sample VOC z-score for focused displays.
- `voc_columns`: explicit VOC list; empty permits numeric discovery.

Patient-level inference is enabled by default within VOC analysis. Configure
`patient_inference`, `patient_permutations`, `patient_seed`,
`patient_min_patients`, `patient_min_nonzero`, and `clr_pseudocount`.
See [VOC statistics](VOC_STATISTICS.md) for normalization, permutation tests,
FDR families, outputs, and limitations of the legacy sample-level plots.

### Study-design modules

- `optional.power_analysis`: cancer/type-group sample-size grids, simulation and
  permutation counts, alpha, seed, transformation and contralateral handling.
- `optional.taxonomy_patient_aware`: taxonomic levels, prevalence, patient/case/type
  columns, contralateral handling, transformation, alpha and top-N display.
- `optional.lung_status_analysis`: patient, case, cancer-site, lung-side and derived
  tumour-side/contralateral/healthy label columns.

#### Power-analysis execution and precision

The sample-type taxonomic power stage uses the CPUs allocated to its Nextflow
task for independent, seed-stable bootstrap simulations. It precomputes patient/
sample-type means and accelerates the small exact signed-rank tests on validated
SciPy 1.17.1; other versions fall back to their native SciPy implementation.
The statistical procedure, simulation counts, permutation counts, and existing
early-stopping rule are unchanged. Cancer PERMANOVA also reuses its centered
distance matrix within each permutation test. Other power stages, including
ISA, remain serial; these optimizations do not eliminate their runtime.

Sample-type taxonomic checkpoints are saved every 25 completed replicates under
`power_analysis/results/.taxonomic_sample_type_checkpoints/` in runtime staging.
They are keyed by inputs, settings, implementation, and library versions. A retry
with the same inputs reuses them regardless of worker count; changed inputs or
settings create new checkpoints. Earlier runs without checkpoints cannot recover
their in-memory progress. Python progress output is unbuffered on new runs.
Standalone execution supports `run_power_analysis_pipeline.sh --workers N`;
Nextflow supplies `--workers` from `task.cpus`, so no YAML parameter is needed.

There is no universal required number of simulations. Keep `n_simulations: 1000`
for the full analysis and the existing smaller mock setting for demonstrations.
For an estimated power of 0.80, 1,000 simulations give Monte Carlo standard error
about 0.013 (approximately ±2.5 percentage points at 95%); 100 give about 0.040
(±7.8 points). This measures simulation precision, not uncertainty from the
original study or adequacy of its sample size. The inner `n_perm` setting is
separate. See [Morris et al. (2019)](https://doi.org/10.1002/sim.8086).

## Networks and MAG Linkage

### `optional.spieceasi`

`min_rel_abund`, `min_prevalence`, and `remove_zero_var` define the network ASV
cohort. `method`, `lambda_min_ratio`, `nlambda`, `rep_num`, `thresh`, `ncores`,
and `seed` configure inference. `edge_threshold` and `keep_negative` control the
reported graph. `force_keep_indicator_asvs` should be enabled only when the
study explicitly requires indicator taxa to bypass the standard cohort filter.

`network_enabled` schedules graph rendering. `modules_enabled`,
`module_methods`, `module_resolutions`, `module_reps`,
`module_consensus_threshold`, and `module_seed` control consensus module
detection. The remaining palette, sizing, layout, and ISA-overlay fields affect
figures rather than inference.

### `optional.network_topology`

`n_null` is the number of seeded degree-preserving configuration-model draws;
`seed` makes them reproducible. `skip_null: true` reports observed topology
without the null comparison. The input is the thresholded SPIEC-EASI graph, so
`optional.spieceasi.enabled` is required.

### `optional.asv_mag_link`

Supply one or more of `genome_qc_dir`, `genome_qc_dirs`, or `barrnap_dir`, plus
genome FASTAs where required. `min_pident` and `min_qcov` control ASV-to-SSU
matching; `threads`, `top_n`, and `plot_top_n` control execution and reporting.
Leave the module disabled when no genome/MAG data exist.

## Summary and Environments

`optional.master_summary` builds integrated ASV tables. Its directory fields identify
the enabled analysis products and `whitelist` specifies files eligible for
direct integration. `max_direct_cols` protects the output from unbounded wide
joins.

`environments` maps logical process groups to committed Conda YAMLs. Most users
should not change these paths. A change alters software provenance and should
be recorded with the run configuration and Git commit.

## Production Preflight Checklist

1. Copy `asv_pipeline_nextflow.yml` and use a new `core.paths.output_dir`.
2. Validate every absolute path and remove unused placeholder paths.
3. Confirm manifest sample IDs against metadata and all external tables.
4. Confirm primer trimming, read length and expected-error settings.
5. Pin taxonomy/SINA references when exact reproducibility matters.
6. Decide explicitly whether MITOMASTER, host filtering, three-tier
   decontamination and batch correction apply.
7. Disable branches lacking required columns or inputs.
8. Record the Git commit and preserve the resolved YAML, manifest, report,
   checksums and logs.
