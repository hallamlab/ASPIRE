# V2 ecological analysis configuration

V2 extends the released ASPIRE preprocessing and shared analyses. The
[shared configuration reference](CONFIGURATION.md) and
[complete parameter catalogue](CONFIG_PARAMETERS.md) cover the merged template.
The sections below describe the additional ecological analyses and their dependencies.

## Environmental thresholds and community states

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


## Comparing community predictors

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


## Functional interpretation of ecological modules

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
