# ASPIRE Configuration Reference

This document covers every public section and key in
`asv_pipeline_nextflow.yml`. Copy the template into a run-specific file. Empty
values mean “not set.” Lists are YAML lists unless a field explicitly uses a
comma-separated string. Unknown keys may be ignored and are not a stable
extension interface.

Values in the template demonstrate the respiratory study and mock workflows.
They are not universal biological defaults. In particular, adapt trimming and
length limits, metadata columns, group labels/orders, palettes, patient-pairing
fields, prevalence cutoffs, and statistical thresholds to the study design.

## Module dependency guide

Core FASTQ processing, ASV generation, SINA trimming, and taxonomy form the
base path. Most later modules consume `FILTER_COUNTS` and/or `PLOT_METADATA`.
Metadata plots require metadata with manifest-matching sample IDs. Diversity,
indicator, clustermap, grouping, power, paired, and network overlays require
their configured grouping columns. `VOC_CORRELATION` additionally needs a VOC
table. `MEASUREMENT_ASSOCIATION` needs measurement columns in metadata or a
joinable external table. MAG branches need genome-QC/barrnap or genome FASTA
inputs. Disable a branch when its conditional inputs are absent.

## Core input and sequence processing

### `paths`, `resources`, and `filename_patterns`

- `paths.input_dir`: directory searched for FASTQs when no manifest is given.
- `paths.output_dir`: required public output directory.
- `paths.manifest`: optional TSV with `sample_id`, `fastq_r1`, and `fastq_r2`.
- `paths.runtime_dir`: persistent runtime root; defaults to
  `<output_dir>/.aspire`.
- `paths.keep_runtime_dir`: retain runtime/cache after success when true.
- `paths.work_dir`, `paths.conda_cache_dir`: optional runtime subdirectory
  overrides. Use storage with sufficient space and safe locking.
- The launcher accepts `--phase all`, `--phase preprocess`, or
  `--phase analysis`; `all` is the default. Preprocessing publishes the
  checksum-validated canonical handoff under
  `<output_dir>/preprocessing_dataset`. Analysis reads only this handoff and
  cannot schedule raw-read processing.
- The supported `--rerun-from STAGE` controller option advances persistent
  per-stage cache generations under `paths.runtime_dir`. It preserves upstream
  task identities and work directories; newly executed tasks remain normally
  cacheable and resumable.
- `cache usage`, `cache list`, `cache prune`, and `cache clear` report or manage
  Nextflow run metadata and work directories. Destructive actions require
  `--force`; `cache prune --keep 10` is dry-run by default and protects the
  recorded phase resume targets. `--prune-cache-after-run 10` applies that
  retention policy after a successful workflow. Conda environments are not
  tied to individual runs and require the separate `--include-conda` option
  with an explicit full-cache clear.
- `resources.sample_threads`: CPUs assigned to each per-sample FASTP, merge,
  and read-filter task. Use `1` to maximize sample concurrency.
- `resources.analysis_threads`: CPU budget requested by dataset-level modules.
- `resources.max_parallel_sample_tasks`: maximum concurrent instances of each
  per-sample read-processing process. Nextflow's CPU scheduler still enforces
  the total available CPU budget.
- `resources.threads`: deprecated compatibility alias for `sample_threads`.
- `resources.single_end`: allow samples without R2 when true.
- `filename_patterns.r1_tokens`, `r2_tokens`: tokens identifying read mates.
- `filename_patterns.ext_patterns`: regular expressions for accepted FASTQs.
- `filename_patterns.sample_strip_regex`: removes lane/read suffixes to form
  sample IDs. Verify generated IDs before matching metadata.

### Read processing and ASV generation

- `fastp.trim_front_r1`, `trim_tail_r1`, `trim_front_r2`, `trim_tail_r2`:
  non-negative base counts removed from read ends; set for the actual assay.
- `merge.max_diffs`: maximum overlap mismatches; `min_overlap`: minimum overlap
  bases; `trunc_quality`: quality truncation threshold; `allow_stagger`: permit
  staggered pairs.
- `table_filter.min_sample_sum`, `min_asv_sum`: non-negative total-count gates.
  `table_filter.script` is an optional developer script override.
- `filter.max_ee`: maximum expected errors; `min_len` and `max_len`: retained
  sequence length bounds in bases.
- `concat.relabel`: relabel sequences; `concat.label_sep`: label separator.
- `unoise.min_size`: minimum UNOISE abundance.
- `swarm.distance`: Swarm clustering distance.

## Alignment, taxonomy, and non-target filtering

### `sina`

`reference` is a local ARB reference and `reference_url` its download fallback;
`download_subdir` names the cache. `regions` lists tested variable regions and
`trim_to` selects the retained region. `batch_size`, `threads`, `keep_gaps`, and
`verbose` control execution and output.

### `taxonomy`

`ref_taxonomy` and `ref_sequences` are local QIIME artifacts; corresponding
`_url` keys are download fallbacks and `_filename` keys name cached files.
`download_subdir` locates the cache. `output_dir`, `output_tsv`, `stats_tsv`,
and `uppercase_fasta` name products; `threads` sets CPUs.

### `mito`

`enabled` enables screening. `run_mitomaster`, `chunk_dir`, `chunk_size`,
`mitomaster_workers`, `mitomaster_retries`, `mitomaster_timeout`, and
`mitomaster_header_mode` control MITOMASTER. `mito_db`/`biof_db` are BLAST
database prefixes; `mito_fasta`/`contaminant_fasta` are preferred FASTA
alternatives. `blast_threads`, `min_pident`, and `min_percov` control BLAST
acceptance. `mitochondria_substring`, `feature_col`, `taxon_col`,
`consensus_col`, `steps`, and `host_first_step` define evidence fields/order.
`output_dir`, `prefix`, `formats`, `figsize`, `style`, `dpi`, and `no_plots`
control products.

### `filter_counts` and `general_stats`

`filter_counts.enabled` enables final filtering. `metadata` and `sample_id_col`
identify samples; `group_col` and `min_group_size` support group prevalence;
`abundance_threshold` and `min_consensus` are filtering cutoffs.
`exclude_taxa` is the rank-aware exclusion list. `taxon_col`, `consensus_col`,
`biofactorial_col`, and `mito_cols` name evidence fields. `output`,
`mito_output_dir`, and `save_intermediates` control products.
`general_stats.enabled` toggles run statistics.

## Metadata and visualization

For all `sample_col` fields, values must match normalized manifest sample IDs.
`color_col` names a metadata color field; `palette_file` supplies an external
two-column palette; `*_palette` accepts `label=#RRGGBB` mappings.

- `sankey`: `enabled`, `metadata`, `sub_dir`, `sample_col`, `group1_col`,
  `color_col`, `palette_file`, `keep_types`, `vertical_order`, `arrangement`,
  `output_prefix`, `title`, `make_labeled`, and `make_unlabeled` configure the
  data-loss diagram.
- `metadata_plots`: `enabled`, `metadata`, `sub_dir`, `sample_col`, `type_col`,
  `color_col`, `palette_file`, `subtraction_group_col`, `subtraction_groups`,
  `keep_types`, `group_order`, `include_rank`, `run_micro`, and `run_mito`.
  BASIN integration fields include `biochem_assignments` and its join/include
  fields, `stratification_timeseries` and its join/include fields, plus
  `cruise_group_assignments`, `cruise_group_meta_join_col`,
  `cruise_group_join_col`, and `cruise_group_include_cols`. Cruise-group fields
  are renamed on import to `cruise_group`, `cruise_group_max_prob`,
  `cruise_group_entropy`, `cruise_group_uncertain`, and optional prefixed PCs.
  Nested `group_normalization.enabled`, `columns`, `pattern`, `replacement`,
  and `preserve_source` normalize group labels while retaining provenance.
- `asv_time_depth_curtain`: `enabled`, `output_dir`, sample/cruise/date/depth
  fields, `group_col`, `palette`, and `group_order` configure a
  hybrid-compartment time-depth curtain. With `source_mode: asv`, the surface
  is estimated from ASPIRE metadata. With `source_mode: basin`, `basin_grid`
  and `basin_cells` provide BASIN's canonical rendered surface and measurement
  positions; ASPIRE preserves both and enlarges positions containing ASV data.
  `base_point_size`, `asv_point_size`, and `asv_point_edge_width` control the
  background BASIN markers and the larger, compartment-colored ASV markers
  with black outlines.
  `renewal_events` and `renewal_date_col` optionally add vertical dashed lines
  for algorithmic renewal-onset dates. `exclude_label_pattern` removes
  reporting-only labels such as `other` or `outlier` from the display support
  surface without deleting sample audit rows. `maximum_depth`, `depth_step_m`,
  `time_subdivisions_per_month` and `time_sigma_months` control the calendar
  grid and longitudinal categorical-support smoothing. Every calendar month in
  each represented year is included. No smoothing is applied across depth;
  observed depth states define midpoint-bounded vertical layers.
  `contour_visual_depth_sigma_m` optionally rounds the displayed categorical
  boundaries only and does not modify assignments or the exported support
  grid. `formats` controls plot export. BASIN mode consumes published BASIN
  data products but does not import or execute BASIN plotting code.
- `collectors_curve`: `enabled`, sample/group/color fields, `group_order`,
  `permutations`, `seed`, `out_prefix`, `title`, `formats`, `xpad`, `max_cols`,
  `show_perms`, and `presence_threshold`.
- `plot_upset`: `enabled`, `sub_dir`, `domain`, `taxonomy_path`, sample/group/
  color fields, `group_order`, `subset_groups`, `skip_venn`, `raw_only`,
  `final_only`, `formats`, and `font_size`.
- `bubbleplotter`: `enabled`, `output_prefix`, count/sample/group/color fields,
  `group1_order`, `formats`, `figsize`, `bubble_scale`, and `no_auto_size`.
- `umap_clustering`: the common output/count/group fields plus `normalize`,
  `transform`, `no_scale`, `n_neighbors`, `min_dist`, `umap_metric`,
  `min_cluster_size`, `min_samples`, and `hdbscan_metric`.

## Batch correction and outliers

### `batch_correction`

`enabled`, `output_dir`, `sample_id_col`, `batch_col`, `asv_orientation`,
`biological_covariates`, `biological_color_col`, `color_palette_col`, and
`biological_palettes` define inputs and biological structure. ConQuR controls
are `conqur_mode`, `conqur_num_core`, `conqur_batch_ref`,
`conqur_logistic_lasso`, `conqur_quantile_type`, `conqur_simple_match`,
`conqur_lambda_quantile`, `conqur_interplt`, `conqur_delta`, and
`conqur_auto_install`. `correction_policy` selects raw, corrected, or `auto`;
auto gates are `auto_min_sample_rho`, `auto_min_bray_rho`,
`auto_max_batch_eta_ratio`, `auto_min_batch_eta_drop`, and
`auto_min_bio_eta_ratio`. Diagnostic clustering uses `umap_neighbors`,
`umap_min_dist`, `hdbscan_min_cluster_size`, `hdbscan_min_samples`,
`hdbscan_selection_method`, `optimize_clustering`, `target_clusters`,
`n_features_plot`, and `random_state`.

Every biological grouping intended for downstream interpretation should be
listed in `biological_covariates`; this prevents correction from treating its
signal as unwanted batch structure. In the Saanich Inlet configuration these
are depth, hybrid oxygen/GMM compartment, cruise group, and nitrate-qualified
renewal phase.

### `outlier_detection`

`enabled`, `output_dir`, `sample_col`, `group_cols`, `transform`,
`asv_orientation`, `pre_transformed`, and `scale` define data. `use_iso`,
`use_svm`, and `use_hdb` select detectors; `vote_threshold` is required
agreement. Isolation Forest uses `iso_contamination`, `iso_estimators`, and
`iso_random_state`; SVM uses `svm_kernel`, `svm_gamma`, and `svm_nu`; HDBSCAN
uses `hdbscan_min_cluster_size`, `hdbscan_min_samples`, and `hdbscan_metric`.

## Ecological and association analyses

### `diversity`

`enabled`, `output_dir`, `mito_output_dir`, `mito_input`, sample/group/color/
block fields, `exclude_groups`, `group_order`, `run_mito`, UMAP settings,
`permanova_perms`, `random_state`, and `verbose` configure the general branch.
Nested `patient_aware` fields enable the paired design and name sample, patient,
case, type, and lung-side fields; `sample_types`, contralateral controls,
`transform`, `permutations`, `seed`, and `require_complete_types` define its
cohort and test.

### `indicspecies`

`enabled`, sample/color fields, `group_cols`, palettes/orders, focus labels,
`block_col`, `blocked_cols`, `perms`, `seed`, `q_threshold`, and `min_n` control tests.
`blocked_cols` lists the grouping analyses that use the repeated-measures design.
If a listed group varies within a block, ASPIRE averages to one block-by-group
profile and restricts permutations within blocks. If it is constant within a
block (for example, a cruise-level environmental group), ASPIRE averages to one
profile per block and permutes the independent block profiles.
`stratified.enabled` and each `analyses` entry (`within_col`, `group_col`,
`levels`) request nested tests. Plot fields are `focus_group1_label`,
`label_focused_asvs`, `plot_enabled`, `plot_pairs_mode`, `plot_output_dir`,
`aligned_plot_enabled`, `aligned_plot_output_dir`, `aligned_alpha`,
`aligned_min_stat`, `aligned_top_n`, `venn`, and `taxonomy`. Legacy
`group1_*`/`group2_*` fields remain supported.

`INDICSPECIES` runs general multigroup `multipatt` (`duleg = FALSE`) only.
Candidate associations include every single group and multigroup proper subset,
but exclude the non-informative union of all represented groups during model
fitting (`max.order = G - 1`). The output reports `association_group_count`, `association_groups`,
`association_scope`, and `indicator_class`; DULEG-specific configuration and
outputs are not generated. Successful runs replace the complete ISA table
directory, removing obsolete files from earlier configurations.

`INDICSPECIES_PLOTS` is deliberately bounded to the primary `group_cols` and
their general multigroup summaries. With `plot_pairs_mode: all`, it
renders all pairs among those primary groupings only; stratified per-state and
pooled summaries are not crossed against unrelated ISA tables. The exact plan
is saved as `isa_plot_plan.tsv`. Stratified ISA remains available to the
network/guild/MAG/function evidence workflows with its state provenance intact.

### `voc_correlation` and `measurement_association`

- `voc_correlation`: `enabled`, metadata/patient/case/type fields,
  `sample_types`, required `voc_table`, `output_dir`, `voc_sample_col`,
  `sample_id_mode`, `use_legacy_voc_subset`, `correlation_direction`
  (`positive`, `negative`, or `both`), and optional `voc_columns`.
- `measurement_association`: `enabled`, `output_dir`, optional
  `measurement_table`, sample/ASV IDs, paired metadata/measurement join lists,
  `measurement_cols`, `exclude_cols`, group/palette, `asv_subset_source`
  (`internal_filters` or `spieceasi_standard_filter`), optional `max_asvs`,
  `min_total`, `min_prevalence`, `top_correlations`, direction,
  `ordination_methods` (`cca,rda,dbrda`), `ordination_measurement_cols`,
  `permutations`, `top_vectors`, and `formats`. `max_asvs: 0` disables a
  fixed-number cap. The
  `spieceasi_standard_filter` source requires `spieceasi.enabled: true` and
  selects only rows with `passes_standard_filter=TRUE`; forced-retention and
  final network-connectivity status do not affect this cohort.
- `titan`: `enabled`, `output_dir`, optional `environmental_variables`,
  `min_split`, `minimum_samples`, `minimum_occurrence`, `minimum_prevalence`,
  `minimum_mean_relative_abundance`, `permutations`,
  `bootstrap_count`, `seed`, `imax`, `iv_total`, `purity_cutoff`,
  `reliability_cutoff`, `ncpus`, `memory`, `plot_enabled`,
  `ranked_plot_top_n`, and `formats`. Enabling TITAN requires both
  `measurement_association.enabled` and `spieceasi.enabled`. An empty
  `environmental_variables` list inherits
  `measurement_association.measurement_cols`; a nonempty list overrides it.
  The ASV cohort is exactly SPIEC-EASI `retained_final`, not a second
  prevalence-ranked subset. `minimum_occurrence` is applied independently
  after restricting samples to observed values for each measurement and cannot
  be below TITAN2's required occurrence of three. `minimum_samples` cannot be
  below TITAN2's hard minimum of 10, while cohorts below 20 are retained with a
  warning. `min_split` must be at least three and also imposes a structural
  requirement of at least twice that many samples. `minimum_prevalence` is the
  fraction of variable-specific measured samples with nonzero abundance, and
  `minimum_mean_relative_abundance` is the taxon's mean relative abundance
  across those samples. Both default to zero and therefore add no filter to the
  inherited SPIEC-EASI cohort. `ncpus` cannot exceed
  `resources.analysis_threads`. `ranked_plot_top_n: 0` displays all tested
  taxa; this affects only the ranked figure, never computation or tables.
- `microbial_compartments`: `enabled`, `output_dir`, `sample_col`, optional
  `transpose`, `minimum_prevalence`, `minimum_max_relative_abundance`,
  `zero_replacement` (`multiplicative` or `pseudocount`),
  `multiplicative_delta`, `pseudocount`, `k_min`, `k_max`,
  `min_cluster_size`, `min_cluster_fraction`, `min_mean_silhouette`,
  `min_stability_ari`, `min_cluster_jaccard`, `stability_replicates`,
  `stability_sample_fraction`, `stability_block_col`,
  `stability_stratum_col`, `stability_primary_quantile`,
  `stability_near_tie_tolerance`, `blocked_robustness_enabled`,
  `prediction_strength_enabled`, `hierarchical_enabled`, `seed`, `ncpus`,
  `environmental_compartment_cols`, `depth_col`, `date_col`,
  `dominant_top_n`, `plot_enabled`, and `formats`. Enabling the branch requires
  SPIEC-EASI and metadata processing. It inherits the ASV-table orientation and
  uses exactly the SPIEC-EASI `retained_final` cohort before any optional
  module-specific filtering. With `multiplicative_delta: 0`, the replacement
  value is calculated as `1 / D^2`, where D is the retained ASV count; every
  replacement is recorded in the zero-replacement audit. PAM is primary and
  average-linkage hierarchical clustering is sensitivity-only. Primary stability
  retains a configured fraction of samples separately within every block, thereby
  preserving the cruise and seasonal composition of every replicate. A PAM result
  is supported only when it satisfies minimum size, the configured lower-tail ARI,
  and minimum per-cluster median-Jaccard criteria. Whole blocks are additionally
  withheld within configured seasonal strata as a robustness analysis; blocked ARI
  and prediction strength are reported but do not select K. The supported K with the
  largest lower-tail ARI is selected, using the configured tolerance to define
  near ties. Prediction strength is an independent held-out diagnostic; silhouette
  is descriptive and does not gate selection. The block identifier controls
  resampling only and is never supplied to PAM. Environmental labels, depth,
  and date are accepted only by the separate post hoc stage and cannot
  influence clustering.
- `microbial_state_interpretation`: `enabled`, `output_dir`, `sample_col`,
  `cruise_col`, `depth_col`, `season_col`, `date_col`,
  `environmental_compartment_cols`, `hybrid_col`, `permutations`,
  `bootstrap_replicates`, `seed`, `hybrid_palette`, `hybrid_order`,
  `mc_palette`, `mc_order`, `asv_top_n`, `selected_asvs`,
  `asv_selection_metric`, `maximum_depth`, `depth_step`,
  `time_subdivisions_per_month`, `time_sigma_months`,
  `depth_visual_sigma_m`, `contour_time_step_days`,
  `contour_max_time_support_days`, `contour_max_depth_support_m`,
  `contour_reference_grid`,
  `exclude_curtain_label_pattern`, and `formats`.
  This post hoc branch requires supported microbial compartments and network
  modules. BASIN labels, depth, season, and time never enter microbial-state
  inference. Correspondence permutations are restricted within cruise, and
  predictive comparisons leave one complete cruise out. The selected-ASV audit
  records whether curtain overlays came from the configured list or the ranked
  list. In addition to the legacy longitudinal curtain, a BASIN-comparable
  microbial-compartment contour is generated by converting each observed PAM
  label to one-hot membership, interpolating those memberships jointly over
  date and depth, renormalizing them, and displaying the maximum membership.
  The three `contour_*` settings define its regular time step and the maximum
  temporal and vertical distances over which interpolation is supported. When
  supplied, `contour_reference_grid` fixes the displayed dates and depths to a
  published BASIN continuous-compartment grid so the figures share identical
  axes; unsupported ASPIRE regions remain gray.
- `ecological_context_atlas`: `enabled`, `output_dir`, `sample_col`,
  `cruise_col`, `depth_col`, `date_col`, `context_cols`, `network_variables`,
  `linkage_q_threshold`, `selected_asvs`, `context_sheet_source`, and `formats`. This reporting-only
  branch runs after TITAN, microbial-compartment inference, measurement/module
  association, compartment/module association, SPIEC-EASI, and ASV–MAG network analysis. Every retained ASV is
  preserved in normalized sample-context, continuous-response, categorical-
  indicator, ASV–MAG/function, and multi-omics tables. `network_variables`
  controls both SPIEC-EASI-only and heterogeneous ASV–MAG TITAN network views
  but never filters tabular results. An empty or omitted list inherits every
  TITAN-configured measurement; measurements without completed TITAN results
  are omitted rather than rendered as empty networks.
  `linkage_q_threshold` (0.05 in the Saanich Inlet run) identifies the reporting
  subset in which the module–compartment, module–measurement, and direct
  ASV–measurement tests all pass their source-process Benjamini–Hochberg
  thresholds. It does not recompute or alter those source statistics.
  `context_sheet_source: mag_linked` renders all eligible MAG-linked ASVs;
  `explicit` renders only `selected_asvs`. No composite evidence score is used.
- `community_turnover`: `enabled`, `output_dir`, `sample_col`, `profile_col`,
  `date_col`, `depth_col`, optional `transpose`, `minimum_prevalence`,
  `minimum_max_relative_abundance`, `zero_replacement`,
  `multiplicative_delta`, `pseudocount`, `distance_metrics`, `primary_metric`,
  `fixed_depth_min_profile_fraction`, `minimum_time_difference_days`,
  `environmental_compartment_cols`, `primary_environmental_compartment_col`,
  `lcbd_permutations`, `boundary_permutations`,
  `boundary_bootstrap_replicates`, `within_profile_lcbd_enabled`,
  `within_profile_lcbd_min_samples`, `seed`, `plot_enabled`, and `formats`.
  The branch requires SPIEC-EASI and metadata processing and independently
  inherits the exact `retained_final` ASV cohort. Fixed-depth temporal turnover
  uses exact numeric equality and retains depths represented in at least
  `fixed_depth_min_profile_fraction` of eligible profiles; no depth mismatch or
  interpolation is permitted. Aitchison is primary, while Bray-Curtis is an
  optional sensitivity result. Global LCBD uses the configured number of
  `adespatial` permutations and BH correction; within-profile LCBD is
  descriptive and requires at least the configured sample count. Cross-boundary
  versus within-compartment turnover is paired at the cruise level; inference
  uses the mean paired difference, a cruise bootstrap confidence interval, and
  a two-sided paired sign-flip test with BH correction within distance metric.

The supplied measurement table is authoritative when a measurement name also
occurs in metadata. Empty, null, and nonnumeric cells remain unmeasured;
negative values are also treated as unmeasured, while zero is retained as a
measured value. ASV and ecological-module Spearman tests use pairwise-complete
samples for each measurement and do not impute chemistry. Constrained
ordinations use complete cases across `ordination_measurement_cols`, which
should contain sufficiently covered core variables rather than sparse
chemistry.

Set `module_analysis_enabled: true` to relate SPIEC-EASI ecological modules to
the same physicochemical measurements. This requires the network and network
module branches. `sparse_measurement_cols` labels passive sparse chemistry,
while `module_cruise_col`, `module_depth_col`, and `module_q_threshold` control
coverage auditing and member-ASV support summaries. Optional `pca_scores`,
`pca_loadings`, `pca_explained`, `hybrid_assignments`, and `hybrid_centroids`
must be supplied together to render the BASIN-aligned environmental biplot.
The module analysis writes pairwise-complete module correlations, per-variable
coverage, member-ASV direction and enrichment summaries, module scores, and
auditable module/vector coordinates. Core vectors use BASIN PCA loadings;
sparse vectors are passive Spearman projections and do not alter the PCA.
It also renders an axis-free environmental network in which each module's
exact within-module SPIEC-EASI topology is laid out locally around that
module's abundance-weighted BASIN PCA centroid. The configured
`asv_mag_network.anchor_top_n` is shared with this renderer: an ASV is an
anchor when it falls within that dense rank for degree, eigenvector centrality,
or betweenness within its module. Between-module edges are retained in the
audit edge table but omitted from the figure. Small display-only shifts prevent
module islands and hybrid centroids from obscuring one another; tables retain
both unmodified environmental centroids and displayed coordinates.
The same fixed layout is also rendered with node area scaled independently by
raw-count prevalence, degree, weighted strength, betweenness, eigenvector
centrality, and among-module participation coefficient. Prevalence is computed
from the post-filter observed count table, not a batch-correction pseudocount
table. `module_network_metric_top_n` controls the number of globally leading
ASVs labeled in each metric plot and reported per metric in the top-player
table. The complete node table retains all metric values, global ranks,
taxonomic assignments, module-anchor flags, and display coordinates.
The sample-level correlations and their BH-adjusted rank-correlation p-values
are descriptive associations; they are not adjusted for depth or repeated
observations within cruise.

### `grouping_diagnostics`

General fields are `enabled`, output/sample/group fields, `baseline_group`,
`primary_group`, palettes/orders, `distance_metrics`, `transform`,
`permutations`, `random_state`, and `formats`. `soft_labeling` uses `enabled`,
`k`, `target_cols`, `exclude_labels`, `min_class_samples`, `distance_quantile`,
`apply_downstream`, `target_col`, `min_confidence`,
`min_neighbor_agreement`, and `min_cv_balanced_accuracy`. Nested `power` uses
`enabled`, `sample_sizes`, `simulations`, `permutations`, `alpha`, and
`min_groups`. Observed labels are not overwritten; excluded labels are neither
training classes nor predictions. For publication-primary analyses, keep
`apply_downstream: false`. Inferred labels and cross-validated performance then
remain an explicit grouping-diagnostics sensitivity analysis, while ISA,
networks, and guild inference use observed BASIN assignments only.

`cruise_level_group_cols` adds categorical cruise labels without treating every
bottle as an independent replicate. The module constructs one equal-depth-
weighted Hellinger water-column profile per cruise using `cruise_col`,
`depth_col`, and `cruise_depth_min_prevalence`. Output rows explicitly report
`analysis_unit=cruise`; bottle groupings report `analysis_unit=bottle`.
Set `require_cruise_level_groups: true` when cruise-level results are required.
Validation then stops immediately if `cruise_level_group_cols` is empty instead
of silently producing bottle-only results.

### `community_predictor_comparison`

`enabled`, `output_dir`, `sample_col`, `cruise_col`, `year_col`, `date_col`, `season_col`,
and `depth_col` define the matched cohorts and blocking structure.
`cruise_depth_min_prevalence` controls which anchored depths enter cruise-level
microbial water-column profiles; missing retained depths are interpolated and
reported in `cruise_depth_profile_coverage.tsv`. `pea_col`
and `centroid_col` select the cruise-level continuous predictors;
`cruise_group_col`, `cruise_group_probability_col`, and
`cruise_group_uncertain_col` select the neutral BASIN cruise assignment and its
confidence fields. `renewal_group_col` selects BASIN's cruise-level,
nitrate-qualified renewal phase (normally `renewal_phase`). The expected
levels are `baseline`, `renewal`, `post-renewal`, and `unknown`; the latter
retains unresolved nitrate coverage or separately classified O2-only anomalies
without promoting them to renewal. The phase is evaluated from one
depth-standardized microbial profile per cruise, with
season adjustment, within-year permutations, PERMDISP, and leave-year-out
cross-validation;
the SI configuration also imports the renewal event identifier, onset and
termination flags, nitrate-support measurements, inference provenance, and raw
O2-anomaly fields for audit, but these fields are not substituted for
`renewal_phase` in the fitted comparisons.
`o2_group_col`, `gmm_group_col`, and `hybrid_group_col` select the bottle-level organizations.
`assignment_source_col` plus `observed_source_label` restrict the categorical
comparison to observed BASIN assignments. `min_group_n` excludes unstable rare
hybrid categories from both matched models. `permutations`, `seed`, and
`formats` configure inference and output. Cruise models use leave-year-out
validation; bottle models use cruise-restricted permutations and
leave-one-cruise-out validation. When `date_col` is available, consecutive-cruise
turnover tables and leave-year-out prediction plots are also generated.
When cruise groups are present, the module compares season, PEA, centroid,
combined PEA+centroid, and cruise-group models on the same cruise cohort, reports
season-adjusted PERMANOVA and PERMDISP for the groups, plots their Bray-Curtis
ordination, and repeats cross-validation after excluding uncertain assignments.
It also tests hybrid-compartment community structure separately within each
environmental cruise group and within each nitrate-qualified renewal phase. These
stratified tests adjust for depth and season,
restrict permutations within cruise, report PERMDISP, and use leave-one-cruise-out
validation. Pairwise compartment contrasts are reported with within-cruise-group
Benjamini-Hochberg correction and pair-specific PERMDISP. The corresponding stratified ISA identifies taxa associated with
compartments inside each cruise state without using ASVs to define the groups.
ISA is run in both directions for both state systems: compartments within each
cruise state, and cruise states within each compartment.

## Study-design analyses

- `power_analysis`: `enabled`, output and sample/patient/case/type fields,
  `sample_sizes_cancer`, `sample_sizes_stype`, `n_simulations`, `n_perm`,
  `alpha`, `seed`, `skip_estimate`, `skip_plot`, `transform`, and contralateral
  controls.
- `taxonomy_patient_aware`: `enabled`, output and cohort fields, `count_col`,
  `tax_levels`, `sample_types`, `min_prevalence`, contralateral/lung-side
  controls, `skip_omnibus`, `transform`, `alpha`, and `top_n`.
- `lung_status_analysis`: `enabled`, output and sample/type/case/patient fields,
  site/side/status fields, `status_a_value`, `status_b_value`,
  `reference_status_value`, `permutations`, and `seed`.

Generalized aliases are `group_power_analysis` for `power_analysis`,
`taxonomy_group_association` for `taxonomy_patient_aware`, and
`paired_group_contrast` for `lung_status_analysis`. Use only one name for each
module in one configuration.

## Clustermaps and networks

### `clustermaps`

Fields cover enable/output/mitochondrial/ISA paths; sample, sample-code, ASV,
count, and up to four grouping columns; orders, exclusions and palettes;
`ranks`, per-rank `topN`, ISA statistic/significance fields and threshold;
`formats`, `figwidth`, `row_height`, `min_height`, `max_height`, and
`mito_sample_mode`.

### `spieceasi`

`enabled`, output/prefix, and `network_enabled` control the branch. Inference
uses `transpose`, abundance/prevalence/zero-variance filters, `method`,
`lambda_min_ratio`, `nlambda`, `rep_num`, `thresh`, `pulsar_criterion`,
`ncores`, and `seed`. Module detection uses `modules_enabled`,
`module_methods`, `module_primary_method`, `module_resolutions`, `module_reps`,
consensus/stability/min-size gates, best-only options, and ISA filters. Graph
inputs/overlays include module/position/node paths, `network_modes`, metadata,
group palettes/orders/focus, and ISA groups. Rendering uses edge, layout,
degree, ISA, and abundance size/scale fields plus `keep_negative` and
`force_filter`, `force_spieceasi`, and `force_graphs`.

`module_subnetworks_enabled` adds one phylum-colored induced subnetwork per
ecological module. These figures retain the exact complete-network SPIEC-EASI
node coordinates in both the figure and audit table; network nodes are never
repositioned. Integers remain centered on isolated nodes, whereas integers for
substantially overlapping nodes are moved to external circular callouts joined
to their source nodes by short leader lines. Taxonomy controls node fill, an analysis-eligible
ASV-genome link controls border width, and node area represents a consortium
prominence index. The index is the equal-weight mean of within-module empirical
percentiles for `module_subnetwork_prominence_metrics`; the SI configuration
uses maximum relative abundance, eigenvector centrality, and participation
coefficient. ASVs at or above `module_subnetwork_prominence_threshold` are
marked as high-prominence ASVs with diamond symbols.
`module_subnetwork_label_top_n` labels the leading index values in addition to
all high-prominence ASVs and every genome-linked ASV. Integer node labels map to an
ASV key containing only phylum, family, and genus. Marker areas are bounded by
`module_subnetwork_prominence_min_area` and
`module_subnetwork_prominence_max_area`. The node-audit and plot-manifest
tables preserve every component percentile, the composite score and rank,
legacy anchor ranks, label counts, encodings, and coordinates.
`module_high_prominence_asvs.tsv` is the thresholded reporting subset;
`module_subnetwork_node_audit.tsv` remains the complete unfiltered inventory.

`pulsar_criterion: stars` requests ordinary StARS. The ASPIRE compatibility
value `bstars` requests bounded StARS by passing `criterion=stars` with lower
and upper StARS bounds enabled; `bstars` itself is never passed to SPIEC-EASI
as an unsupported criterion.

### MAG and summary modules

- `genome_cooccurrence`: independently infers metagenomic and metatranscriptomic
  species-level proportionality co-occurrence networks from the non-normalized
  `read_count` field in the configured genome-level long-format recruitment
  tables. With `input_feature_level: genome`, `genome_col` identifies genome
  records and `genome_qc` supplies their species mapping. A complete
  genome-by-sample grid is constructed, absent combinations are assigned zero,
  positive counts below `minimum_read_count` are set to zero, and the remaining
  values are averaged across every eligible genome within each species. Thus,
  species values represent mean per-genome recruitment rather than summed or
  deduplicated species read counts. `genome_qc` also supplies the shared quality/taxonomy table;
  `expected_species` is a fail-fast assertion on the deterministic one-genome-
  per-species selection; and the sample/genome/value fields identify the exact
  abundance columns. Species-level counts drive inference, while representative
  selection supplies one stable node identifier and its metadata. Selection
  prioritizes MIMAG high over
  medium, then greater assembly length, greater completeness, lower
  contamination, greater N50, and genome identifier. With
  `input_scale: raw_counts`, the assay-specific matrices retain raw counts until
  inference. Each sample is then closed to unit sum, zeros are replaced
  multiplicatively using `zero_replacement_fraction` times that sample's
  smallest positive part, and the composition is CLR transformed. Pairwise
  association is Lovell/Erb proportionality, rho-p. `min_rel_abund` and
  `min_prevalence` control species filtering. Whole cruises are resampled with
  replacement within season for `bootstrap_iterations`; an edge must reach
  `min_abs_rho`, `min_bootstrap_recovery`, and `min_sign_consistency`. Empirical
  p-values use `permutations` independent species-wise shuffles within the
  configured `permutation_strata` and are Benjamini-Hochberg adjusted; `max_q`
  is the adjusted-significance cutoff. `seed` makes resampling, permutation,
  community detection, and plotting deterministic. `closure_total` applies
  only to the optional `input_scale: closed_abundance` preparation mode and is
  not used by the configured SI analysis. Degree-zero species are retained in
  the 24-species GraphML files and plots so a lack of supported co-occurrence is
  not hidden. Positive edges indicate concordant relative-recruitment profiles;
  negative edges indicate contrasting profiles. Neither direction establishes
  direct ecological interaction. Descriptive native proportionality clusters use
  absolute selected-edge weights and are labeled `MG-PC#` or `MT-PC#`. Separate
  fixed-layout overlays report phylum, family, mean abundance, proportionality
  cluster, degree, betweenness, eigenvector centrality,
  closeness, and descriptive annotations propagated through eligible ASV-
  genome links: ASV ecological module, hybrid ISA, microbial-compartment ISA,
  and TITAN response. Propagated overlays are annotations of linked ASVs, not
  direct genome-level tests.

- `asv_mag_link`: `enabled`, optional `master_tsv`, `genome_qc_dir` or
  `genome_qc_dirs`, `id_token_indexes`, `barrnap_dir` or `genome_fasta_dir`,
  `output_dir`, `threads`, `min_pident`, `min_qcov`, optional inclusive
  `asv_taxonomy_min_confidence` (0-1), optional hard
  `min_completeness` (inclusive) and `max_contamination` (exclusive) MAG
  quality gates, optional `gunc_assessment_value` requiring both
  `gunc_assessment` and `gunc_strict_assessment` to equal that value,
  optional `require_species_assignment` to exclude MAGs without a classified
  species, optional `min_rrna_marker_count` requiring each of the 16S, 23S,
  and 5S rRNA marker counts to meet the configured minimum,
  `top_n`, and `plot_top_n`. When configured, missing assessment columns are
  fatal rather than silently retaining MAGs. A usable retained 16S/barrnap
  reference remains required independently of the selected quality gate.
- `asv_mag_network`: `enabled`, output/prefix, `graph_variant`, `anchor_top_n`, link identity/
  coverage gates, taxonomy sources, `mag_id_mode`, optional `mag_abundance`,
  its format and genome/sample/value fields, optional `cruise_metadata`
  containing cruise dates for all-metagenome monthly recruitment profiles,
  `min_shared_samples`,
  `mag_knn` (up to this many strongest positive and this many strongest
  negative abundance-correlation neighbors selected per MAG within an
  ASV-linked MAG subnetwork; the undirected union is plotted),
  `abundance_transform`, `mag_abundance_normalization`,
  `functional_module_min_fraction`, and
  `functional_annotations`.

`asv_mag_network.anchor_top_n` is shared by the network renderer and the
module-anchor summary. ASVs are ranked independently within each ecological
module by degree, eigenvector centrality, and betweenness centrality. An ASV is
an anchor when its dense rank is at most `anchor_top_n` for at least one of the
three metrics; it is not required to qualify for all three. The default is 1,
and metric ties are retained.

`asv_mag_network.biochemical_groupings` selects biochemical ISA overlays for
the locked ecological-module network (for example,
`[o2_subcompartment_final, cruise_group]`), and `isa_q_threshold` controls
which indicator associations receive a colored node rim. Node fill remains the
ecological-module palette, so fill/rim overlap can be compared without changing
ASV positions or edges. Single-group indicators use that biochemical level's
configured color, multigroup indicators use a black rim, and non-significant
ASVs retain a thin gray rim. These shared, modality-independent panels are
written as
`asv_mag_network_ecological_modules_by_<grouping>.{pdf,png,svg}` and indexed in
the network panel manifest.

Direct grouping diagnostics are configured with
`grouping_diagnostic_sample_col`, `grouping_diagnostic_cruise_col`,
`grouping_diagnostic_groups`, `grouping_diagnostic_cruise_level_groups`, and
`grouping_diagnostic_permutations`. The default comparison set is
`Depth`, `Season`, legacy `o2_compartment`, `gmm_component`, hybrid
`o2_subcompartment_final`, and `cruise_group`. Put labels that describe an
entire cruise, normally `Season` and `cruise_group`, in
`grouping_diagnostic_cruise_level_groups`; those tests aggregate to one profile
per cruise rather than treating bottles from the same cruise as independent. Feature
selection is independent of ASV–MAG linkage: the ASV model uses every processed
ASV, the MetaG model uses every MAG recruitment profile, and the MetaT model
uses every MAG expression profile. Joint ASV+MetaG, ASV+MetaT, and
ASV+MetaG+MetaT models use all features from their component blocks on the
samples shared by those blocks. Each standardized omics block receives equal
total weight in joint distances so feature count alone cannot make one modality
dominate.

Bottle-level grouping tests retain samples as the analysis units and restrict
label permutations and cross-validation folds by cruise. Cruise-level tests first
average all sampled depths into one MAG/community profile per cruise, then
permute and cross-validate those independent cruise profiles. Outputs in
`asv_mag_network/group_diagnostics/` include multivariate effect size,
permutation significance, blocked nearest-centroid balanced accuracy,
chance-adjusted balanced accuracy, classification permutation significance,
per-feature effect/FDR tables, group-mean feature heatmaps, and a common
grouping-comparison figure. Chance adjustment permits predictive comparisons
among groupings with different numbers of levels. These results test whether each grouping captures
repeatable community, population-abundance, or expression structure; they do
not require the grouping to agree with every individual feature.

For nonexclusive MAG recruitment counts, supply the SeqKit paired-end summary
used for recruitment as `mag_abundance_seqkit` and, independently for RNA,
`mag_transcript_abundance_seqkit`. ASPIRE parses common `_pe.1/_pe.2`,
`_R1/_R2`, and `_1/_2` filename suffixes, requires equal R1/R2 `num_seqs`, and
uses the number of read pairs—not their sum—as the input-fragment denominator.
With `mag_abundance_normalization: auto`, a supplied SeqKit table selects
input-fragment FPM; without one, ASPIRE stops rather than silently changing the
normalization. Set `provided_fpkm` or `provided_tpm` to validate and retain
an upstream-normalized value column without further library-size correction.
Median-ratio normalization remains available only through an explicit
`median_ratio` setting. Do not directly divide values having different units;
the ecological-context atlas accepts paired TPM/TPM log2-ratio and standard-error
matrices through `rna_dna_log2_tpm_ratio` and `rna_dna_log2_tpm_se`.
After the configured transformation, each MAG profile is standardized across
samples. MAG-MAG subnetwork edges are Spearman correlations of these normalized
profiles. FPM is noncompositional: sums across genomes may exceed one million
when a fragment maps to multiple genomes. `relative` is therefore retained only
for mutually exclusive assignments, while `raw` skips sequencing-depth
normalization.
- `master_summary`: `enabled`, `output_dir`, source-module directories,
  `max_direct_cols`, and the filename `whitelist` eligible for integration.

## `environments`

Every value is a process-specific Conda YAML. Keys are `main`, `sina`,
`taxonomy`, `mitomaster`, `mito_checker`, `filter_counts`, `general_stats`,
`sankey`, `plot_metadata`, `asv_time_depth_curtain`, `batch_correction`, `outlier_checker`,
`collectors_curve`, `plot_upset`, `bubbleplotter`, `umap_clustering`,
`diversity`, `indicspecies`, `voc_correlation`, `measurement_association`, `titan`,
`grouping_diagnostics`, `community_predictor_comparison`,
`group_label_augmentation`, `clustermaps`,
`power_analysis`, `taxonomy_patient_aware`, `lung_status_analysis`,
`spieceasi`, `network_modules`, `network`, `master_summary`, `asv_mag_link`, and
`asv_mag_network`. Keep committed values for ordinary runs. Overrides are for
dependency development and change the reproducibility environment.

## Remaining field glossary

The following less-common fields are listed explicitly so configuration review
does not require searching the workflow source:

- Taxonomy cache names: `ref_taxonomy_url`, `ref_taxonomy_filename`,
  `ref_sequences_url`, and `ref_sequences_filename`.
- General join and identifier fields: `metadata_sample_col`, `asv_id_col`,
  `metadata_join_cols`, `measurement_join_cols`, and `sample_code_col`.
- Patient/lung design fields: `cancer_site_col`, `contralateral_col`,
  `contralateral_sample_types`, `contralateral_value`,
  `exclude_contralateral_in_cancer`, `keep_contralateral_in_cancer`,
  `lung_code_col`, `lung_side_col`, `lung_status_col`, `tumor_side_col`, and
  `healthy_col`.
- Clustermap fields: `isa_file`, `exclude_group1`, `group2_col`,
  `group3_col`, `group4_col`, `group1_palette`, `group2_palette`,
  `group3_palette`, `group4_palette`, `group2_order`, `isa_min_stat`,
  `isa_significance_cols`, `isa_stat_cols`, and `plot_asv_level`.
  Set `plot_asv_level: false` to retain `clustermap_ASV_ID_plot.tsv` while
  skipping the expensive ASV-resolution figures; higher taxonomic-rank
  clustermaps are unaffected.
- Shared grouping maps: `group_palettes`, `group_orders`, and `focus_labels`.
- Master-summary source locations: `clustermaps_dir`, `indicspecies_dir`,
  `spieceasi_dir`, and `asv_mag_dir`.
- Network module fields: `module_seed`, `module_consensus_threshold`,
  `module_best_only`, `module_best_min_size`, `module_best_min_stability`,
  `module_best_top_n`,
  `module_isa_only`, `module_color_by_isa`, `module_isa_source`,
  `module_isa_min_stat`, and `module_isa_max_q`.
- Optional network input/overlay fields: `modules_sub`, `modules_all`,
  `graph_pos_sub`, `graph_pos_all`, `node_features`, and `isa_overlay_groups`.
- Network filtering/rendering fields: `min_rel_abund`, `min_prevalence`,
  `force_keep_isa_asvs`, `remove_zero_var`,
  `edge_threshold`, `edge_width_scale`, `layout_iters`, `layout_seed`,
  `layout_scale`, `degree_scale`, `degree_size_mode`, `degree_min_area`,
  `isa_scale`, `abundance_size_mode`, `abundance_reference`,
  `abundance_reference_area`, `abundance_min_area`, `abundance_max_area`, and
  `abundance_scale_power`. `max_labels` caps collision-adjusted labels per
  network plot (default `100`; use `0` only when an unlimited labeled plot is
  known to fit in memory).

`min_rel_abund` is the minimum maximum within-sample relative abundance and
`min_prevalence` is the minimum fraction of samples containing an ASV; an ASV
must satisfy both. `force_keep_isa_asvs` defaults to `false`. Enabling it restores
significant primary-ISA ASVs that fail network filtering and can substantially
increase model size, so it should be treated as an explicit sensitivity choice
rather than a default.

For example, `min_rel_abund: 0.005` and `min_prevalence: 0.01` retain an ASV only
when it is observed in at least 1% of samples and reaches at least 0.5% relative
abundance in one or more samples. The abundance parameter is therefore a
minimum threshold applied to each ASV's maximum observed relative abundance; it
is not an upper abundance limit.

Degree-zero ASVs in the primary ASV network remain documented in the filtering
and adjacency tables but are removed from GraphML files, module detection, and
rendered networks because they have no inferred association. The explicitly
bounded genome networks retain their degree-zero species for audit. Multigroup
ISA nodes are rendered as
segmented circles using the exact component colors shown in the legend; colors
are never averaged into an unlisted blended color.

When `module_best_only` is enabled, `module_best_min_size` and
`module_best_min_stability` are hard eligibility filters. Eligible modules are
then ranked using internal edge weight, edge density, conductance, weighted
modularity contribution, and mean node stability. `module_best_top_n` limits
highlighted module figures to the strongest N modules by this composite quality
ranking. The complete module assignment and diagnostic tables are still
written, so this setting changes visualization focus rather than discarding
module membership evidence. When `module_best_only` is disabled (the default),
every inferred module receives its own label and color in both the SPIEC-EASI
and downstream ASV–MAG ecological-module figures. The quality-ranking table is
still written for audit and review but does not filter the displayed modules.

`ncores` controls the number of independent pulsar/StARS R workers. ASPIRE
forces BLAS, OpenMP, BLIS, MKL, NumExpr, and RcppParallel to one thread inside
each worker, preventing nested numerical-library threading from oversubscribing
the host. Thus `ncores: 25` requests approximately 25 computational workers,
not 25 workers each attempting to use every CPU.

For the SI production network, `thresh: 0.05` is the StARS instability target,
`lambda_min_ratio: 0.01` expands the fitted sparsity path, and
`edge_threshold: 0.0` retains the nonzero edges selected by SPIEC-EASI rather
than applying a second arbitrary partial-correlation cutoff. Review the selected
lambda index and stability path after fitting; a boundary solution indicates
that the lambda path should be widened or sampled more finely.
- MAG taxonomy and abundance schema fields: `asv_taxonomy_source`,
  `mag_taxonomy_source`, `mag_abundance_format`, `mag_abundance_genome_col`,
  `mag_abundance_sample_col`, and `mag_abundance_value_col`.

Path/column fields name inputs or exact case-sensitive columns. Boolean fields
toggle the behavior named by the key. Size, scale, threshold, and iteration
fields are numeric and should retain template values unless the corresponding
analysis design is being deliberately changed.

## Production preflight checklist

- All active paths exist and reference compatible files.
- Manifest IDs are unique, pairing is correct, and metadata IDs match exactly.
- Every active group, color, patient, lung-side, and measurement field exists.
- Palette labels and orders cover observed metadata values.
- Trimming, merge, length, and error parameters match the assay.
- Optional modules without their conditional inputs are disabled.
- A representative run has been reviewed for read retention, filtering,
  taxonomy, sample accounting, and metadata joins before the full cohort run.
### `group_guild_function`

Integrates ISA group evidence with ecological network modules, ASV–MAG mappings,
and MAG functional modules. It requires `indicspecies`, `network_modules`, and
`asv_mag_network`. `groupings` controls all tabled comparisons, while
`primary_groupings` limits only the focused plot and does not remove results.
Set `target_modules: [auto_n_s]` to audit all KEGG modules whose definitions
describe environmental nitrogen or sulfur transformations and retain every
module represented in at least one linked MAG at
`asv_mag_network.functional_module_min_fraction`. The selection and excluded
modules are recorded in `tables/nitrogen_sulfur_module_selection.tsv`.
Explicit KEGG module IDs remain supported when a fixed hypothesis set is
required.

Supplying `pca_scores`, `pca_explained`, `hybrid_assignments`, and
`hybrid_centroids` adds a BASIN PC1–PC2 overlay for the ecological-network
ASVs; the workflow also consumes the ASPIRE taxonomy-stage output. All four
external BASIN paths must be supplied together. Each network ASV is positioned
at the relative-abundance-weighted barycenter of the BASIN sample coordinates
where it was observed. Multiple ASV libraries mapped to the same cruise-depth
coordinate are consolidated by their median relative abundance before
positioning. For display, these ASV positions are pooled into one
abundance-weighted point per ecological-module/phylum combination; no
abundance or prevalence cutoff is applied. The plotted coordinates are
therefore observed distribution summaries, not predicted niches.
`ecological_module_asv_pca_positions.tsv` retains every ASV-level coordinate,
while `ecological_module_phylum_pca_positions.tsv` records the plotted
module–phylum coordinates, constituent ASVs, and abundance weights.

Both standalone and stratified ISA tables are propagated. Every stratified row
retains `stratified_within_col` and `stratified_within_value` through the
indicator-ASV, ecological-network-module, MAG, and target-function evidence
tables. Network-module abundance is also tested for hybrid-compartment effects
inside each cruise group and nitrate-qualified renewal phase using cruise-restricted
permutations.

Three complementary outputs quantify ecological-module agreement with the
biochemical groupings. `group_ecological_module_association.tsv` reports the
permutation-tested omnibus effect size (`effect_eta_squared`) and FDR-adjusted
significance for each module/grouping pair.
`group_level_ecological_module_abundance.tsv` reports effect direction using
group-level means, medians, differences, and ratios.
`indicator_asv_module_enrichment.tsv` tests whether significant indicator ASVs
are overrepresented within particular ecological modules. Finally,
`ecological_module_biochemical_group_concordance.tsv` reports Cramér's V and a
permutation p-value across significant single-group indicators. Cramér's V is
an overlap measure, not evidence that the classifications are interchangeable;
multigroup indicators remain visible in the network but are excluded from this
one-label contingency statistic.

The mapping input contains every sequence- and taxonomy-eligible MAG candidate.
Ambiguous ASVs therefore remain in the evidence chain; `eligible_mag_count` and
`asv_link_weight` preserve mapping multiplicity for auditing and weighted
summaries. Unique-only mappings remain available as a sensitivity subset.

The ASV–MAG heterogeneous network retains the complete SPIEC-EASI ASV network
as pale context and adds every eligible MAG link. MAG-paired ASVs are emphasized
as larger gray circles and MAGs are smaller black circles. ASV–MAG identity is
defined by the sequence evidence; Spearman correlation of relative-abundance
profiles is attached to each sequence edge as supporting evidence, not used to
invent or discard identity links. Terminal ASV sample suffixes such as
`_plate8` are removed automatically only when doing so increases overlap with
the MAG recruitment sample IDs. The selected harmonization mode and overlap
counts are recorded in
`qc/asv_mag_network_metagenome_abundance_sample_harmonization.tsv`.
Metagenome recruitment produces
`network/asv_mag_network_metagenome_heterogeneous_network.{pdf,png,svg}`.
When `mag_transcript_abundance` is configured, an otherwise identical
metatranscriptome network is written as
`network/asv_mag_network_metatranscriptome_heterogeneous_network.{pdf,png,svg}`.
The former uses MAG–MAG abundance-profile correlations; the latter uses
MAG–MAG expression-profile correlations. Both retain the same SPIEC-EASI ASV
backbone and sequence-derived ASV–MAG links. Modality is explicit in every
associated network, correlation-table, and lineage-profile filename.

`max_modules` limits guild abundance-heatmap figures to the highest-ranked
ecological modules; it does not remove modules from the evidence tables.
`heatmap_top_asvs` caps the all-ASV heatmap at the most variable ASVs within
those modules so labels remain reviewable. Guild analysis does not redraw a
separate network. Its canonical network panels are emitted by
`ASV_MAG_NETWORK`, preserving the exact SPIEC-EASI ASV backbone, eligible
ASV–MAG links, MAG–MAG associations, central ASV coordinates, external MAG
subnetwork coordinates, ellipses, and anchor lines used by the base
heterogeneous-network panel.

`asv_mag_network.ambiguity_target_modules` selects functions displayed in the
one-ASV-to-many-MAG figures and functional network overlays. With
`[auto_n_s]`, the same data-driven N/S coverage audit described above is used.
The network directory contains one ecological-module/anchor overlay and one
KEGG contributor overlay per selected module and abundance modality. Every
overlay is rendered by the same specialized heterogeneous-network renderer,
so only its evidence styling changes. Contributor tables and
`network/asv_mag_network_network_panel_manifest.tsv` provide the evidence and
file index for assembling a multi-panel figure. PDF, PNG, and SVG ambiguity
output is generated for every ASV having at least two eligible MAG candidates.
Dense mappings are paginated at 10 MAGs per page without filtering candidates.
No example is automatically selected for the main text.

In each KEGG contributor panel, linked contributing ASVs retain the same
selected ecological-module colors used in panels 1–2. A dark black outline
marks functional contribution on those ASVs. The external MAG subnetworks keep
the exact MAG-node and MAG-association styling of the corresponding base
heterogeneous network. A contributing ASV's external ellipse and its straight
ASV-to-ellipse anchor are bold black, but nothing inside that subnetwork is
restyled. All ASVs in selected ecological modules retain pale versions of their
module colors as structural context, while contributing ASVs use the saturated
version; unselected modules remain pale gray. Thus hue identifies
ecological-module membership and the black outline independently identifies
KEGG-module evidence. Orange and blue remain reserved for positive and negative
MAG abundance or expression associations.

MAG marker area is constant across every external subnetwork; marker size does
not vary with subnetwork density or any unreported quantity. Ecological-module
legend entries use natural numeric order (`M1`, `M2`, …, `M10`), and each
linked-ASV callout includes its module label, for example `ASV1 [M1]`.

Monthly lineage-recruitment profiles are generated for every
analysis-eligible ASV–MAG-linked ASV, including ASVs with only one eligible MAG
link. The two-or-more-MAG requirement applies only to the ambiguity diagnostic
figures and tables. Each monthly profile renders the observed cruise-level
measurements and a cyclic smoothed trace for every linked lineage without
confidence bands, keeping single- and multi-lineage comparisons uncluttered.

When ASV counts are available, these profiles also underlay the linked ASV's
relative abundance on a labeled secondary percent axis. The ASV calculation
uses every ASV sample having a resolvable cruise identifier and date, not only
the samples overlapping the MetaG or MetaT recruitment table. Counts are
converted to within-sample relative abundance, averaged across sampled depths
within each cruise, and then displayed as cruise observations plus a cyclic
smoothed gray curve on a secondary percent axis. The cruise-level values are
written beside each figure as
`*_asv_relative_abundance_by_cruise.tsv`.

The canonical network figure progression is explicitly indexed in the
manifest: (1) all ASVs colored only by ecological module, with no MAG cue or
subnetwork; (2) the same ASV view with centrality-based ecological anchors
enlarged; (3) the full ASV–MAG heterogeneous network; and (4) one otherwise
identical heterogeneous-network overlay for every selected N/S KEGG module.
Panels 1–2 use the same seeded ASV coordinates as panels 3–4 but crop to the
ASV network so the MAG-free views remain readable.

Because panels 1–2 do not use MAG abundance or expression, they are emitted
once with modality-neutral filenames
`asv_mag_network_ecological_modules.*` and
`asv_mag_network_ecological_modules_and_anchors.*`. MetaG/MetaT filename
variants are reserved for panels whose MAG–MAG associations actually use the
corresponding abundance or expression profiles.

For panels 1–2, only modules marked `is_best=true` by the existing
`module_best_min_size`, `module_best_min_stability`, and
`module_best_top_n` selection are assigned distinct colors; other modules
remain as gray network context. An ecological anchor is an ASV ranked first
within its selected module for degree, eigenvector centrality, or betweenness
centrality. In panel 2, ordinary ASVs are small and all anchors share one
larger size. A dark black outline on that large marker indicates that the
anchor has at least one analysis-eligible MAG link; anchors without such a
link retain the same size with only the standard light outline.

### ASV-level ecological-module PCA overlays

With the four BASIN PCA inputs configured, `group_guild_function` creates the original module–phylum plot plus `ecological_module_asv_pca_overlay_asv_level` and `ecological_module_asv_pca_overlay_asv_level_abundant_asvs_labeled` in each configured format. Both new plots show individual observed ASV positions. The labeled variant places ASV IDs locally within the point axes, searches nearby collision-free positions with short noncrossing leader lines, and adds a diamond-keyed phylum/family/genus legend. Layout constraints reject intersecting leaders, overlapping labels, labels over scatter markers, and leaders through other labels. Point coordinates never move. A `.layout.json` audit accompanies the labeled figure.

```yaml
group_guild_function:
  pca_label_min_mean_abundance_pct: 1.0
```

The threshold is a percentage: 1.0 selects ASVs with mean relative abundance **strictly greater than 1%**, calculated over the full input ASV count matrix before restricting to plotted network members. It does not select the top 1% of ranked ASVs. Valid values are 0 (all positive-mean ASVs) through values below 100. The default is 1.0. `tables/ecological_module_asv_pca_overlay_labels.tsv` records selected ASVs, abundances, coordinates, and whether a position was available. Individual ASV coordinates are abundance-weighted environmental PCA positions, not a new community ordination. Module–phylum colored circles in the original plot are pooled positions.

The ASV legend also reports module ID and “Module peak”: the hybrid compartment with the highest aggregate module mean in `compartment_matched_ecological_module_association.tsv`. This is a module-level maximum, not necessarily the individual ASV maximum. The selection audit records the module peak label and mean.

The ASV taxonomy legend is grouped under bold module/peak headers in numeric module order (M1, M2, …, M10). Member ASVs are ordered by their numeric ASV identifier within each module, retaining the colored diamonds and phylum/family/genus assignments.

### Optional matched cohort for diversity

`diversity.matched_cohort.table` optionally supplies an exact sample cohort TSV
using `diversity.sample_col`. The existing microbial diversity tests
use that cohort and write to the normal diversity outputs. Descriptive plots,
UMAP, and abundance inputs retain their normal sample coverage. Missing requested
samples cause an error. Leave it null for general ASPIRE sample eligibility.
The SI configuration supplies the saved 310-sample BASIN community cohort;
this setting does not restrict other ASPIRE processes or mitochondrial analyses.
