# Complete ASPIRE Configuration Parameter Catalogue

This file documents every explicit parameter in the canonical
[`asv_pipeline_nextflow.yml`](../asv_pipeline_nextflow.yml) template. It is generated
from that template; lists are documented as one parameter and their complete template
value is shown. Paths are resolved relative to the YAML file unless absolute.

The template value is a starting point, not a universal scientific recommendation.
Study-specific thresholds, metadata columns, references, and optional modules must be
chosen and reported for the study design. See [Configuration Reference](CONFIGURATION.md)
for dependencies and interpretation and [Process Reference](PROCESS_REFERENCE.md) for I/O.

## `core`

### `core.paths`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `core.paths.input_dir` | str | `/abs/path/input_fastqs` | Directory searched for FASTQ files when no manifest is supplied. |
| `core.paths.output_dir` | str | `/abs/path/output_dir` | Public output directory published atomically after a successful run. |
| `core.paths.manifest` | str | `/abs/path/sample_manifest.tsv` | Optional TSV declaring sample IDs and R1/R2 FASTQ paths. |
| `core.paths.runtime_dir` | str | `/abs/path/output_dir/.aspire` | Persistent private runtime root for work, staging, and environment caches. |
| `core.paths.keep_runtime_dir` | bool | `true` | Retain the private runtime after success so resume and targeted reruns remain possible. |
| `core.paths.work_dir` | null | `null` | Optional override for the Nextflow work directory; null derives it from runtime_dir. |
| `core.paths.conda_cache_dir` | null | `null` | Optional override for the per-process Conda cache; null derives it from runtime_dir. |

### `core.resources`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `core.resources.threads` | int | `6` | CPU threads requested by each parallel per-sample read-processing task. |
| `core.resources.single_end` | bool | `false` | Treat inputs as single-end reads instead of the default paired-end design. |

### `core.filename_patterns`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `core.filename_patterns.r1_tokens` | list | `["R1",1]` | Tokens used to recognize read-1 FASTQ filenames during discovery. |
| `core.filename_patterns.r2_tokens` | list | `["R2",2]` | Read-2 tokens paired positionally with r1_tokens. |
| `core.filename_patterns.ext_patterns` | list | `["\\.fastq\\.gz$","\\.fq\\.gz$","\\.fastq$","\\.fq$"]` | Regular expressions accepted as FASTQ filename suffixes. |
| `core.filename_patterns.sample_strip_regex` | str | `(_S[0-9]+)?(_L[0-9]{3})?(_R[12])?(_[12])?(_001)?$` | Regular expression removed from discovered filenames to derive sample IDs. |

### `core.fastp`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `core.fastp.trim_front_r1` | int | `19` | Number of bases removed from the 5′ end of read 1. |
| `core.fastp.trim_tail_r1` | int | `20` | Number of bases removed from the 3′ end of read 1. |
| `core.fastp.trim_front_r2` | int | `20` | Number of bases removed from the 5′ end of read 2. |
| `core.fastp.trim_tail_r2` | int | `20` | Number of bases removed from the 3′ end of read 2. |

### `core.merge`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `core.merge.max_diffs` | int | `20` | Maximum accepted diffs for the merge module. |
| `core.merge.min_overlap` | int | `5` | Minimum accepted overlap for the merge module. |
| `core.merge.trunc_quality` | int | `5` | Quality score used by the read-merging truncation rule. |
| `core.merge.allow_stagger` | bool | `true` | Permit staggered paired-read alignments during merging. |

### `core.table_filter`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `core.table_filter.min_sample_sum` | int | `5000` | Minimum total reads required for a sample to survive the early technical count filter. |
| `core.table_filter.min_asv_sum` | int | `0` | Minimum study-wide count required for an ASV to survive the early technical count filter. |
| `core.table_filter.script` | str | `processes/filter_table/filter_ASV_table.py` | Repository-relative implementation used for early count-table filtering. |

### `core.filter`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `core.filter.max_ee` | float | `1.0` | Maximum expected errors accepted for a merged read. |
| `core.filter.min_len` | int | `245` | Minimum accepted merged-read length in bases. |
| `core.filter.max_len` | int | `1500` | Maximum accepted merged-read length in bases. |

### `core.concat`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `core.concat.relabel` | bool | `true` | Rewrite concatenated FASTA headers so each sequence retains its sample identity. |
| `core.concat.label_sep` | str | `:` | Separator placed between sample labels and original sequence identifiers. |

### `core.unoise`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `core.unoise.min_size` | int | `3` | Minimum dereplicated abundance supplied to UNOISE denoising. |

### `core.swarm`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `core.swarm.distance` | int | `1` | Sequence-distance setting used for SWARM-style clustering/mapping behavior. |

### `core.sina`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `core.sina.reference` | null | `null` | Local reference artifact; null permits configured download behavior. |
| `core.sina.reference_url` | str | `https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/ARB_files/SILVA_138.2_SSURef_NR99_03_07_24_opt.arb.gz` | Download URL used when the local reference is unavailable. |
| `core.sina.download_subdir` | str | `sina_reference` | Directory or published subdirectory used for download subdir. |
| `core.sina.regions` | list | `["V4-V5","V5","V4"]` | Ordered values used for regions by the sina module. |
| `core.sina.trim_to` | str | `V4-V5` | Value selecting or naming trim to for the sina module. |
| `core.sina.batch_size` | int | `1000000` | Configured number or size for batch size in the sina module. |
| `core.sina.keep_gaps` | bool | `false` | Retain gaps when true. |
| `core.sina.threads` | int | `6` | CPU threads requested by this module. |
| `core.sina.verbose` | bool | `false` | Emit additional module diagnostics. |

### `core.taxonomy`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `core.taxonomy.ref_taxonomy` | null | `null` | Value selecting or naming ref taxonomy for the taxonomy module. |
| `core.taxonomy.ref_taxonomy_url` | str | `https://data.qiime2.org/2024.10/common/silva-138-99-tax.qza` | Download URL used when the local ref taxonomy is unavailable. |
| `core.taxonomy.ref_taxonomy_filename` | str | `silva-138_2-ssu-nr99-tax.qza` | Local filename assigned to the retrieved ref taxonomy artifact. |
| `core.taxonomy.ref_sequences` | null | `null` | Value selecting or naming ref sequences for the taxonomy module. |
| `core.taxonomy.ref_sequences_url` | str | `https://data.qiime2.org/2024.10/common/silva-138-99-seqs.qza` | Download URL used when the local ref sequences is unavailable. |
| `core.taxonomy.ref_sequences_filename` | str | `silva-138_2-ssu-nr99-seqs-DNA.qza` | Local filename assigned to the retrieved ref sequences artifact. |
| `core.taxonomy.download_subdir` | str | `taxonomy_reference` | Directory or published subdirectory used for download subdir. |
| `core.taxonomy.output_dir` | str | `taxonomy` | Directory or published subdirectory used for output dir. |
| `core.taxonomy.output_tsv` | str | `ASV_SILVA_tax.full-length.vsearch.tsv` | Value selecting or naming output tsv for the taxonomy module. |
| `core.taxonomy.stats_tsv` | str | `ASV_SILVA_stats.full-length.vsearch.tsv` | Value selecting or naming stats tsv for the taxonomy module. |
| `core.taxonomy.uppercase_fasta` | str | `ASVs.upper.fasta` | Path to the uppercase FASTA reference; null means it is not supplied. |
| `core.taxonomy.threads` | int | `6` | CPU threads requested by this module. |

## `standard`

### `standard.mito`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `standard.mito.enabled` | bool | `true` | Whether this module or nested analysis is scheduled. |
| `standard.mito.run_mitomaster` | bool | `true` | Contact MITOMASTER for mitochondrial evidence; disable for an intentionally offline run. |
| `standard.mito.chunk_dir` | null | `null` | Directory or published subdirectory used for chunk dir. |
| `standard.mito.chunk_size` | int | `10` | Configured number or size for chunk size in the mito module. |
| `standard.mito.mitomaster_workers` | int | `8` | Numeric setting for mitomaster workers in the mito module. |
| `standard.mito.mitomaster_retries` | int | `4` | Numeric setting for mitomaster retries in the mito module. |
| `standard.mito.mitomaster_timeout` | int | `90` | Numeric setting for mitomaster timeout in the mito module. |
| `standard.mito.mitomaster_header_mode` | str | `first` | Value selecting or naming mitomaster header mode for the mito module. |
| `standard.mito.output_dir` | str | `mito/mitomap` | Directory or published subdirectory used for output dir. |
| `standard.mito.prefix` | str | `nontarget` | Value selecting or naming prefix for the mito module. |
| `standard.mito.formats` | str | `svg,pdf` | Comma-separated output figure formats. |
| `standard.mito.min_pident` | float | `97.0` | Minimum BLAST percent identity for local non-target evidence. |
| `standard.mito.min_percov` | float | `51.0` | Minimum BLAST query coverage percentage for local non-target evidence. |
| `standard.mito.mitochondria_substring` | str | `mitochondria` | Value selecting or naming mitochondria substring for the mito module. |
| `standard.mito.feature_col` | str | `Feature ID` | Column containing ASV/feature identifiers. |
| `standard.mito.taxon_col` | str | `Taxon` | Column containing taxonomy strings. |
| `standard.mito.consensus_col` | str | `Consensus` | Column containing taxonomy consensus scores. |
| `standard.mito.steps` | str | `BioFactorial,Qiime_NB_FULL,MITOMASTER,BLAST_mito` | Value selecting or naming steps for the mito module. |
| `standard.mito.host_first_step` | str | `BioFactorial` | Value selecting or naming host first step for the mito module. |
| `standard.mito.figsize` | str | `10x6` | Figure width and height specification. |
| `standard.mito.style` | str | `whitegrid` | Plotting style applied to generated figures. |
| `standard.mito.dpi` | int | `300` | Raster figure resolution in dots per inch. |
| `standard.mito.no_plots` | bool | `false` | Enable or disable no plots behavior in the mito module. |
| `standard.mito.blast_threads` | int | `6` | Numeric setting for blast threads in the mito module. |
| `standard.mito.mito_db` | str | `/abs/path/blastdb/mito_ncbi` | BLAST database prefix used for mito screening. |
| `standard.mito.biof_db` | str | `/abs/path/blastdb/ssu_pipeline_contaminants` | BLAST database prefix used for biof screening. |
| `standard.mito.mito_fasta` | null | `null` | Path to the mito FASTA reference; null means it is not supplied. |
| `standard.mito.contaminant_fasta` | null | `null` | Path to the contaminant FASTA reference; null means it is not supplied. |

### `standard.filter_counts`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `standard.filter_counts.enabled` | bool | `true` | Whether this module or nested analysis is scheduled. |
| `standard.filter_counts.metadata` | str | `/abs/path/metadata.tsv` | Path to the metadata TSV consumed by this module. |
| `standard.filter_counts.output` | str | `ASV_target.tsv` | Output filename written by this module. |
| `standard.filter_counts.group_col` | str | `Type_Group` | Primary grouping column used by this module. |
| `standard.filter_counts.min_group_size` | int | `3` | Minimum number of samples required for a metadata group to participate in group-aware filtering. |
| `standard.filter_counts.abundance_threshold` | float | `0.5` | Minimum relative-abundance percentage used by the final abundance filter; 0.5 means 0.5%. |
| `standard.filter_counts.sample_id_col` | str | `Sample` | Column containing sample identifiers. |
| `standard.filter_counts.min_consensus` | float | `0.0` | Minimum taxonomy consensus score accepted by final count filtering. |
| `standard.filter_counts.exclude_taxa` | list | `["Species:Homo sapiens","Class:Mammalia"]` | Exact rank-qualified taxa removed even when they pass other filters. |
| `standard.filter_counts.taxon_col` | str | `Taxon` | Column containing taxonomy strings. |
| `standard.filter_counts.consensus_col` | str | `Consensus` | Column containing taxonomy consensus scores. |
| `standard.filter_counts.biofactorial_col` | str | `BioFactorial` | Input-table column containing biofactorial. |
| `standard.filter_counts.mito_cols` | list | `["MITOMASTER","BLAST_mito"]` | Input-table columns used for mito. |
| `standard.filter_counts.save_intermediates` | bool | `true` | Enable or disable save intermediates behavior in the filter counts module. |
| `standard.filter_counts.mito_output_dir` | str | `mito/ASVs` | Directory or published subdirectory used for mito output dir. |

### `standard.general_stats`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `standard.general_stats.enabled` | bool | `true` | Whether this module or nested analysis is scheduled. |

### `standard.metadata_plots`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `standard.metadata_plots.enabled` | bool | `true` | Whether this module or nested analysis is scheduled. |
| `standard.metadata_plots.metadata` | str | `/abs/path/metadata.tsv` | Path to the metadata TSV consumed by this module. |
| `standard.metadata_plots.sub_dir` | str | `.` | Directory or published subdirectory used for sub dir. |
| `standard.metadata_plots.sample_col` | str | `Sample` | Column containing sample identifiers. |
| `standard.metadata_plots.type_col` | str | `Type_Group` | Column containing sample-type labels. |
| `standard.metadata_plots.color_col` | str | `Color` | Metadata column containing or selecting display colors. |
| `standard.metadata_plots.palette_file` | null | `null` | Path to, or configured name of, the palette file input. |
| `standard.metadata_plots.subtraction_group_col` | str | `Type_Group` | Input-table column containing subtraction group. |
| `standard.metadata_plots.subtraction_groups` | list | `["Scope Flush","Skin Brush","Control"]` | Ordered values used for subtraction groups by the metadata plots module. |
| `standard.metadata_plots.keep_types` | list | `["Oral Rinse","BAL","Bronchial Brush"]` | Retain types when true. |
| `standard.metadata_plots.group_order` | list | `["Oral Rinse","BAL","Bronchial Brush"]` | Explicit category order used for group analysis and display. |
| `standard.metadata_plots.include_rank` | list | `[]` | Ordered values used for include rank by the metadata plots module. |
| `standard.metadata_plots.input_table` | str | `filtered` | Path to, or configured name of, the input table input. |
| `standard.metadata_plots.min_pre_correction_sample_sum` | int | `0` | Minimum accepted pre correction sample sum for the metadata plots module. |
| `standard.metadata_plots.drop_zero_asvs` | bool | `false` | Enable or disable drop zero ASVs behavior in the metadata plots module. |
| `standard.metadata_plots.run_micro` | bool | `true` | Run micro when true. |
| `standard.metadata_plots.run_mito` | bool | `true` | Run mito when true. |
| `standard.metadata_plots.group_normalization.enabled` | bool | `false` | Whether this module or nested analysis is scheduled. |
| `standard.metadata_plots.group_normalization.columns` | list | `[]` | Ordered values used for columns by the metadata plots module. |
| `standard.metadata_plots.group_normalization.pattern` | null | `null` | Value selecting or naming pattern for the metadata plots module. |
| `standard.metadata_plots.group_normalization.replacement` | str | `outlier` | Value selecting or naming replacement for the metadata plots module. |
| `standard.metadata_plots.group_normalization.preserve_source` | bool | `true` | Enable or disable preserve source behavior in the metadata plots module. |

## `optional`

### `optional.three_tier_decontam`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `optional.three_tier_decontam.enabled` | bool | `false` | Whether this module or nested analysis is scheduled. |
| `optional.three_tier_decontam.output_dir` | str | `three_tier_decontam` | Directory or published subdirectory used for output dir. |
| `optional.three_tier_decontam.metadata` | null | `null` | Path to the metadata TSV consumed by this module. |
| `optional.three_tier_decontam.metadata_sample_col` | str | `Sample` | Sample-identifier column in the metadata table. |
| `optional.three_tier_decontam.sample_col` | str | `Sample` | Column containing sample identifiers. |
| `optional.three_tier_decontam.negative_control_col` | str | `is_negative_control` | Input-table column containing negative control. |
| `optional.three_tier_decontam.positive_control_col` | str | `is_positive_control` | Input-table column containing positive control. |
| `optional.three_tier_decontam.negative_control_labels` | list | `["PBS","PBS_twz","Negative_96","Negative_man"]` | Ordered values used for negative control labels by the three tier decontam module. |
| `optional.three_tier_decontam.positive_control_labels` | list | `["Positive_96","Positive_man"]` | Ordered values used for positive control labels by the three tier decontam module. |
| `optional.three_tier_decontam.concentration_col` | str | `DNA_conc` | Input-table column containing concentration. |
| `optional.three_tier_decontam.type_col` | str | `Type_Group` | Column containing sample-type labels. |
| `optional.three_tier_decontam.sample_types` | list | `["BAL","Bronchial Brush","Oral Rinse"]` | Ordered values used for sample types by the three tier decontam module. |
| `optional.three_tier_decontam.pooled_threshold` | float | `0.1` | Decontam score cutoff for the pooled control-based tier. |
| `optional.three_tier_decontam.within_type_threshold` | float | `0.1` | Decontam score cutoff applied within each configured sample type. |
| `optional.three_tier_decontam.aggressive_threshold` | float | `0.5` | Secondary, more permissive contaminant score cutoff used by the configured combination rule. |
| `optional.three_tier_decontam.combine_mode` | str | `min` | Rule used to combine pooled and within-type contaminant evidence. |
| `optional.three_tier_decontam.biological_plausibility` | bool | `true` | Apply the taxonomy/biological-plausibility tier before producing downstream tables. |

### `optional.sankey`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `optional.sankey.enabled` | bool | `true` | Whether this module or nested analysis is scheduled. |
| `optional.sankey.metadata` | str | `/abs/path/metadata.tsv` | Path to the metadata TSV consumed by this module. |
| `optional.sankey.sub_dir` | str | `.` | Directory or published subdirectory used for sub dir. |
| `optional.sankey.sample_col` | str | `Sample` | Column containing sample identifiers. |
| `optional.sankey.group1_col` | str | `Type_Group` | Primary display or analysis grouping column. |
| `optional.sankey.color_col` | str | `Color` | Metadata column containing or selecting display colors. |
| `optional.sankey.palette_file` | null | `null` | Path to, or configured name of, the palette file input. |
| `optional.sankey.keep_types` | null | `null` | Retain types when true. |
| `optional.sankey.vertical_order` | list | `["Bronchial Brush","BAL","Oral Rinse","Skin Brush","Scope Flush","Control"]` | Explicit category order used for vertical analysis and display. |
| `optional.sankey.arrangement` | str | `freeform` | Value selecting or naming arrangement for the sankey module. |
| `optional.sankey.output_prefix` | str | `metadata/data_loss_sankey` | Output name, prefix, or destination used for output prefix. |
| `optional.sankey.title` | str | `Data Loss Flow` | Title printed on generated figures. |
| `optional.sankey.make_labeled` | bool | `true` | Enable or disable make labeled behavior in the sankey module. |
| `optional.sankey.make_unlabeled` | bool | `true` | Enable or disable make unlabeled behavior in the sankey module. |

### `optional.batch_correction`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `optional.batch_correction.enabled` | bool | `false` | Whether this module or nested analysis is scheduled. |
| `optional.batch_correction.output_dir` | str | `batch_correction` | Directory or published subdirectory used for output dir. |
| `optional.batch_correction.sample_id_col` | str | `Sample` | Column containing sample identifiers. |
| `optional.batch_correction.batch_col` | str | `DNA_plate` | Input-table column containing batch. |
| `optional.batch_correction.asv_orientation` | str | `features_rows` | Value selecting or naming ASV orientation for the batch correction module. |
| `optional.batch_correction.biological_covariates` | str | `Type_Group` | Value selecting or naming biological covariates for the batch correction module. |
| `optional.batch_correction.biological_color_col` | str | `Type_Group` | Input-table column containing biological color. |
| `optional.batch_correction.color_palette_col` | str | `Color` | Input-table column containing color palette. |
| `optional.batch_correction.conqur_mode` | str | `standard` | Value selecting or naming conqur mode for the batch correction module. |
| `optional.batch_correction.conqur_num_core` | int | `6` | Numeric setting for conqur num core in the batch correction module. |
| `optional.batch_correction.conqur_batch_ref` | null | `null` | Value selecting or naming conqur batch ref for the batch correction module. |
| `optional.batch_correction.conqur_logistic_lasso` | bool | `true` | Enable or disable conqur logistic lasso behavior in the batch correction module. |
| `optional.batch_correction.conqur_quantile_type` | str | `lasso` | Value selecting or naming conqur quantile type for the batch correction module. |
| `optional.batch_correction.conqur_simple_match` | bool | `true` | Enable or disable conqur simple match behavior in the batch correction module. |
| `optional.batch_correction.conqur_lambda_quantile` | str | `2p/logn` | Value selecting or naming conqur lambda quantile for the batch correction module. |
| `optional.batch_correction.conqur_interplt` | bool | `false` | Enable or disable conqur interplt behavior in the batch correction module. |
| `optional.batch_correction.conqur_delta` | float | `0.4999` | Numeric setting for conqur delta in the batch correction module. |
| `optional.batch_correction.conqur_auto_install` | bool | `true` | Enable or disable conqur auto install behavior in the batch correction module. |
| `optional.batch_correction.umap_neighbors` | int | `15` | Numeric setting for umap neighbors in the batch correction module. |
| `optional.batch_correction.umap_min_dist` | float | `0.1` | Numeric setting for umap min dist in the batch correction module. |
| `optional.batch_correction.hdbscan_min_cluster_size` | int | `5` | Configured number or size for hdbscan min cluster size in the batch correction module. |
| `optional.batch_correction.hdbscan_min_samples` | null | `null` | Value selecting or naming hdbscan min samples for the batch correction module. |
| `optional.batch_correction.hdbscan_selection_method` | str | `eom` | Value selecting or naming hdbscan selection method for the batch correction module. |
| `optional.batch_correction.optimize_clustering` | bool | `false` | Enable or disable optimize clustering behavior in the batch correction module. |
| `optional.batch_correction.target_clusters` | str | `3-8` | Value selecting or naming target clusters for the batch correction module. |
| `optional.batch_correction.n_features_plot` | int | `5` | Configured number or size for n features plot in the batch correction module. |
| `optional.batch_correction.random_state` | int | `42` | Random seed used to make stochastic behavior reproducible. |
| `optional.batch_correction.biological_palettes` | mapping | `{}` | Explicit label-to-color mapping used for biological displays. |
| `optional.batch_correction.correction_policy` | str | `always` | Select downstream counts: always uses corrected counts, never uses raw counts, and auto applies the count-space/batch/biological preservation gates. The template preserves the original always-correct behavior. |
| `optional.batch_correction.auto_min_sample_rho` | float | `0.85` | Minimum median per-sample Spearman count-preservation correlation required by automatic correction selection. |
| `optional.batch_correction.auto_min_bray_rho` | float | `0.75` | Minimum Bray-Curtis distance correlation required by automatic correction selection. |
| `optional.batch_correction.auto_max_batch_eta_ratio` | float | `0.95` | Numeric setting for auto max batch eta ratio in the batch correction module. |
| `optional.batch_correction.auto_min_batch_eta_drop` | float | `0.01` | Numeric setting for auto min batch eta drop in the batch correction module. |
| `optional.batch_correction.auto_min_bio_eta_ratio` | float | `0.7` | Numeric setting for auto min bio eta ratio in the batch correction module. |

### `optional.outlier_detection`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `optional.outlier_detection.enabled` | bool | `true` | Whether this module or nested analysis is scheduled. |
| `optional.outlier_detection.output_dir` | str | `outliers_corrected` | Directory or published subdirectory used for output dir. |
| `optional.outlier_detection.sample_col` | str | `Sample` | Column containing sample identifiers. |
| `optional.outlier_detection.group_cols` | list | `["type_group,Case"]` | Input-table columns used for group. |
| `optional.outlier_detection.transform` | str | `none` | Transformation applied to input abundance values. |
| `optional.outlier_detection.asv_orientation` | str | `Samples_rows` | Value selecting or naming ASV orientation for the outlier detection module. |
| `optional.outlier_detection.pre_transformed` | bool | `true` | Enable or disable pre transformed behavior in the outlier detection module. |
| `optional.outlier_detection.scale` | bool | `false` | Enable or disable scale behavior in the outlier detection module. |
| `optional.outlier_detection.use_iso` | bool | `true` | Enable or disable use iso behavior in the outlier detection module. |
| `optional.outlier_detection.use_svm` | bool | `true` | Enable or disable use svm behavior in the outlier detection module. |
| `optional.outlier_detection.use_hdb` | bool | `true` | Enable or disable use hdb behavior in the outlier detection module. |
| `optional.outlier_detection.vote_threshold` | int | `3` | Statistical or filtering cutoff for vote threshold in the outlier detection module. |
| `optional.outlier_detection.iso_contamination` | str | `auto` | Value selecting or naming iso contamination for the outlier detection module. |
| `optional.outlier_detection.iso_estimators` | int | `200` | Numeric setting for iso estimators in the outlier detection module. |
| `optional.outlier_detection.iso_random_state` | int | `42` | Numeric setting for iso random state in the outlier detection module. |
| `optional.outlier_detection.svm_kernel` | str | `rbf` | Value selecting or naming svm kernel for the outlier detection module. |
| `optional.outlier_detection.svm_gamma` | str | `scale` | Value selecting or naming svm gamma for the outlier detection module. |
| `optional.outlier_detection.svm_nu` | float | `0.1` | Numeric setting for svm nu in the outlier detection module. |
| `optional.outlier_detection.hdbscan_min_cluster_size` | int | `5` | Configured number or size for hdbscan min cluster size in the outlier detection module. |
| `optional.outlier_detection.hdbscan_min_samples` | int | `5` | Numeric setting for hdbscan min samples in the outlier detection module. |
| `optional.outlier_detection.hdbscan_metric` | str | `euclidean` | Value selecting or naming hdbscan metric for the outlier detection module. |

### `optional.collectors_curve`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `optional.collectors_curve.enabled` | bool | `true` | Whether this module or nested analysis is scheduled. |
| `optional.collectors_curve.sample_col` | str | `Sample` | Column containing sample identifiers. |
| `optional.collectors_curve.group_col` | str | `Type_Group` | Primary grouping column used by this module. |
| `optional.collectors_curve.color_col` | str | `Color` | Metadata column containing or selecting display colors. |
| `optional.collectors_curve.group_order` | list | `["Oral Rinse","BAL","Bronchial Brush"]` | Explicit category order used for group analysis and display. |
| `optional.collectors_curve.permutations` | int | `999` | Number of permutations used by the statistical test. |
| `optional.collectors_curve.seed` | int | `42` | Random seed used to make stochastic behavior reproducible. |
| `optional.collectors_curve.out_prefix` | str | `metadata/collectors_curve` | Output name, prefix, or destination used for out prefix. |
| `optional.collectors_curve.title` | str | `Collector's Curves` | Title printed on generated figures. |
| `optional.collectors_curve.formats` | str | `pdf,svg` | Comma-separated output figure formats. |
| `optional.collectors_curve.xpad` | float | `0.5` | Numeric setting for xpad in the collectors curve module. |
| `optional.collectors_curve.max_cols` | int | `3` | Input-table columns used for max. |
| `optional.collectors_curve.show_perms` | int | `10` | Numeric setting for show perms in the collectors curve module. |
| `optional.collectors_curve.presence_threshold` | float | `0.0` | Statistical or filtering cutoff for presence threshold in the collectors curve module. |

### `optional.plot_upset`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `optional.plot_upset.enabled` | bool | `true` | Whether this module or nested analysis is scheduled. |
| `optional.plot_upset.sub_dir` | str | `.` | Directory or published subdirectory used for sub dir. |
| `optional.plot_upset.domain` | str | `micro` | Value selecting or naming domain for the plot upset module. |
| `optional.plot_upset.taxonomy_path` | null | `null` | Value selecting or naming taxonomy path for the plot upset module. |
| `optional.plot_upset.sample_id_col` | str | `Sample` | Column containing sample identifiers. |
| `optional.plot_upset.group_col` | str | `Type_Group` | Primary grouping column used by this module. |
| `optional.plot_upset.color_col` | str | `Color` | Metadata column containing or selecting display colors. |
| `optional.plot_upset.group_order` | list | `["Oral Rinse","BAL","Bronchial Brush"]` | Explicit category order used for group analysis and display. |
| `optional.plot_upset.subset_groups` | null | `null` | Value selecting or naming subset groups for the plot upset module. |
| `optional.plot_upset.skip_venn` | bool | `false` | Skip venn when true. |
| `optional.plot_upset.raw_only` | bool | `false` | Enable or disable raw only behavior in the plot upset module. |
| `optional.plot_upset.final_only` | bool | `false` | Enable or disable final only behavior in the plot upset module. |
| `optional.plot_upset.formats` | str | `pdf,svg,png` | Comma-separated output figure formats. |
| `optional.plot_upset.font_size` | int | `10` | Configured number or size for font size in the plot upset module. |

### `optional.bubbleplotter`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `optional.bubbleplotter.enabled` | bool | `true` | Whether this module or nested analysis is scheduled. |
| `optional.bubbleplotter.output_prefix` | str | `metadata/bubble_plot_asv` | Output name, prefix, or destination used for output prefix. |
| `optional.bubbleplotter.count_col` | str | `corr_count` | Column containing ASV counts or transformed count values. |
| `optional.bubbleplotter.sample_col` | str | `Sample` | Column containing sample identifiers. |
| `optional.bubbleplotter.group1_col` | str | `Type_Group` | Primary display or analysis grouping column. |
| `optional.bubbleplotter.group2_col` | str | `Case` | Secondary display or analysis grouping column. |
| `optional.bubbleplotter.color_col` | str | `Color` | Metadata column containing or selecting display colors. |
| `optional.bubbleplotter.group1_order` | list | `["Oral Rinse","BAL","Bronchial Brush"]` | Explicit category order used for group1 analysis and display. |
| `optional.bubbleplotter.formats` | str | `pdf,png,svg` | Comma-separated output figure formats. |
| `optional.bubbleplotter.figsize` | str | `32,60` | Figure width and height specification. |
| `optional.bubbleplotter.bubble_scale` | int | `10` | Numeric setting for bubble scale in the bubbleplotter module. |
| `optional.bubbleplotter.no_auto_size` | bool | `false` | Configured number or size for no auto size in the bubbleplotter module. |

### `optional.umap_clustering`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `optional.umap_clustering.enabled` | bool | `true` | Whether this module or nested analysis is scheduled. |
| `optional.umap_clustering.output_prefix` | str | `metadata/umap_clustering` | Output name, prefix, or destination used for output prefix. |
| `optional.umap_clustering.sample_col` | str | `Sample` | Column containing sample identifiers. |
| `optional.umap_clustering.count_col` | str | `corr_count` | Column containing ASV counts or transformed count values. |
| `optional.umap_clustering.group1_col` | str | `Type_Group` | Primary display or analysis grouping column. |
| `optional.umap_clustering.group2_col` | str | `Case` | Secondary display or analysis grouping column. |
| `optional.umap_clustering.color_col` | str | `Color` | Metadata column containing or selecting display colors. |
| `optional.umap_clustering.group1_order` | list | `["Oral Rinse","BAL","Bronchial Brush"]` | Explicit category order used for group1 analysis and display. |
| `optional.umap_clustering.group1_palette` | null | `null` | Explicit label-to-color mapping used for group1 displays. |
| `optional.umap_clustering.group2_palette` | null | `null` | Explicit label-to-color mapping used for group2 displays. |
| `optional.umap_clustering.formats` | str | `pdf,png,svg` | Comma-separated output figure formats. |
| `optional.umap_clustering.normalize` | str | `clr` | Normalization method applied before this analysis or visualization. |
| `optional.umap_clustering.transform` | str | `none` | Transformation applied to input abundance values. |
| `optional.umap_clustering.no_scale` | bool | `true` | Enable or disable no scale behavior in the umap clustering module. |
| `optional.umap_clustering.n_neighbors` | int | `12` | Configured number or size for n neighbors in the umap clustering module. |
| `optional.umap_clustering.min_dist` | float | `0.05` | Minimum accepted dist for the umap clustering module. |
| `optional.umap_clustering.umap_metric` | str | `cosine` | Value selecting or naming umap metric for the umap clustering module. |
| `optional.umap_clustering.min_cluster_size` | int | `8` | Minimum accepted cluster size for the umap clustering module. |
| `optional.umap_clustering.min_samples` | int | `4` | Minimum accepted samples for the umap clustering module. |
| `optional.umap_clustering.hdbscan_metric` | str | `euclidean` | Value selecting or naming hdbscan metric for the umap clustering module. |

### `optional.diversity`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `optional.diversity.enabled` | bool | `true` | Whether this module or nested analysis is scheduled. |
| `optional.diversity.output_dir` | str | `diversity` | Directory or published subdirectory used for output dir. |
| `optional.diversity.mito_output_dir` | str | `mito/diversity` | Directory or published subdirectory used for mito output dir. |
| `optional.diversity.mito_input` | null | `null` | Path to, or configured name of, the mito input input. |
| `optional.diversity.sample_col` | str | `Sample` | Column containing sample identifiers. |
| `optional.diversity.group1_col` | str | `Type_Group` | Primary display or analysis grouping column. |
| `optional.diversity.group2_col` | str | `Case` | Secondary display or analysis grouping column. |
| `optional.diversity.color_col` | str | `Color` | Metadata column containing or selecting display colors. |
| `optional.diversity.block_col` | str | `Participant_ID` | Column defining repeated-measure or permutation blocks. |
| `optional.diversity.exclude_groups` | list | `[]` | Values or groups excluded according to groups. |
| `optional.diversity.group_order` | list | `["Oral Rinse","BAL","Bronchial Brush"]` | Explicit category order used for group analysis and display. |
| `optional.diversity.group1_palette` | null | `null` | Explicit label-to-color mapping used for group1 displays. |
| `optional.diversity.group2_palette` | null | `null` | Explicit label-to-color mapping used for group2 displays. |
| `optional.diversity.run_mito` | bool | `true` | Run mito when true. |
| `optional.diversity.umap_neighbors` | int | `30` | Numeric setting for umap neighbors in the diversity module. |
| `optional.diversity.umap_min_dist` | float | `0.01` | Numeric setting for umap min dist in the diversity module. |
| `optional.diversity.permanova_perms` | int | `999` | Numeric setting for permanova perms in the diversity module. |
| `optional.diversity.random_state` | int | `42` | Random seed used to make stochastic behavior reproducible. |
| `optional.diversity.verbose` | bool | `true` | Emit additional module diagnostics. |
| `optional.diversity.patient_aware.enabled` | bool | `false` | Whether this module or nested analysis is scheduled. |
| `optional.diversity.patient_aware.output_dir` | str | `patient_aware` | Directory or published subdirectory used for output dir. |
| `optional.diversity.patient_aware.sample_col` | str | `Sample` | Column containing sample identifiers. |
| `optional.diversity.patient_aware.patient_col` | str | `Participant_ID` | Column containing participant/patient identifiers for blocking or pairing. |
| `optional.diversity.patient_aware.case_col` | str | `Case` | Column containing case/control status. |
| `optional.diversity.patient_aware.type_col` | str | `Type_Group` | Column containing sample-type labels. |
| `optional.diversity.patient_aware.sample_types` | str | `Oral Rinse,BAL,Bronchial Brush` | Value selecting or naming sample types for the diversity module. |
| `optional.diversity.patient_aware.exclude_contralateral_in_cancer` | bool | `true` | Values or groups excluded according to contralateral in cancer. |
| `optional.diversity.patient_aware.contralateral_col` | str | `lung_status` | Input-table column containing contralateral. |
| `optional.diversity.patient_aware.cancer_site_col` | str | `Cancer_Site` | Input-table column containing cancer site. |
| `optional.diversity.patient_aware.lung_side_col` | str | `lung_code` | Input-table column containing lung side. |
| `optional.diversity.patient_aware.contralateral_value` | str | `Contralateral` | Value selecting or naming contralateral value for the diversity module. |
| `optional.diversity.patient_aware.contralateral_sample_types` | str | `Bronchial Brush,BAL` | Value selecting or naming contralateral sample types for the diversity module. |
| `optional.diversity.patient_aware.transform` | str | `none` | Transformation applied to input abundance values. |
| `optional.diversity.patient_aware.permutations` | int | `9999` | Number of permutations used by the statistical test. |
| `optional.diversity.patient_aware.seed` | int | `42` | Random seed used to make stochastic behavior reproducible. |
| `optional.diversity.patient_aware.require_complete_types` | bool | `false` | Enable or disable require complete types behavior in the diversity module. |

### `optional.indicspecies`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `optional.indicspecies.enabled` | bool | `true` | Whether this module or nested analysis is scheduled. |
| `optional.indicspecies.sample_col` | str | `Sample` | Column containing sample identifiers. |
| `optional.indicspecies.color_col` | str | `Color` | Metadata column containing or selecting display colors. |
| `optional.indicspecies.group_cols` | list | `["Type_Group","Case"]` | Input-table columns used for group. |
| `optional.indicspecies.stratified.enabled` | bool | `true` | Whether this module or nested analysis is scheduled. |
| `optional.indicspecies.stratified.analyses` | list | `[{"within_col":"Type_Group","group_col":"Case","levels":["Oral Rinse","BAL","Bronchial Brush"]}]` | Ordered values used for analyses by the indicspecies module. |
| `optional.indicspecies.group_palettes.Type_Group` | str | `Oral Rinse=#6A3D9A,BAL+Oral Rinse=#E78AC3,BAL=#0072B2,BAL+Bronchial Brush=#5CC8C8,Bronchial Brush=#009E73,Bronchial Brush+Oral Rinse=#B8E186,BAL+Bronchial Brush+Oral Rinse=#CBB6E9,not_indicator=#D3D3D3` | Value selecting or naming Type Group for the indicspecies module. |
| `optional.indicspecies.group_palettes.Case` | str | `Non-Cancer=#FFFFFF,Cancer=#A50026,Cancer+Non-Cancer=#000000,not_indicator=#D3D3D3` | Value selecting or naming Case for the indicspecies module. |
| `optional.indicspecies.group_orders.Type_Group` | str | `Oral Rinse,BAL+Oral Rinse,BAL,BAL+Bronchial Brush,Bronchial Brush,Bronchial Brush+Oral Rinse,BAL+Bronchial Brush+Oral Rinse,not_indicator` | Value selecting or naming Type Group for the indicspecies module. |
| `optional.indicspecies.group_orders.Case` | str | `Non-Cancer,Cancer,Cancer+Non-Cancer` | Value selecting or naming Case for the indicspecies module. |
| `optional.indicspecies.focus_labels` | mapping | `{}` | Value selecting or naming focus labels for the indicspecies module. |
| `optional.indicspecies.group1_order` | list | `["Oral Rinse","BAL+Oral Rinse","BAL","BAL+Bronchial Brush","Bronchial Brush","Bronchial Brush+Oral Rinse","BAL+Bronchial Brush+Oral Rinse","not_indicator"]` | Explicit category order used for group1 analysis and display. |
| `optional.indicspecies.group2_order` | list | `["Non-Cancer","Cancer","Cancer+Non-Cancer"]` | Explicit category order used for group2 analysis and display. |
| `optional.indicspecies.group1_palette` | str | `Oral Rinse=#6A3D9A,BAL+Oral Rinse=#E78AC3,BAL=#0072B2,BAL+Bronchial Brush=#5CC8C8,Bronchial Brush=#009E73,Bronchial Brush+Oral Rinse=#B8E186,BAL+Bronchial Brush+Oral Rinse=#CBB6E9,not_indicator=#D3D3D3` | Explicit label-to-color mapping used for group1 displays. |
| `optional.indicspecies.group2_palette` | str | `Non-Cancer=#FFFFFF,Cancer=#A50026,Cancer+Non-Cancer=#000000,not_indicator=#D3D3D3` | Explicit label-to-color mapping used for group2 displays. |
| `optional.indicspecies.block_col` | str | `Participant_ID` | Column defining repeated-measure or permutation blocks. |
| `optional.indicspecies.perms` | int | `9999` | Number of permutations used by the statistical test. |
| `optional.indicspecies.seed` | int | `42` | Random seed used to make stochastic behavior reproducible. |
| `optional.indicspecies.q_threshold` | float | `0.05` | Statistical or filtering cutoff for q threshold in the indicspecies module. |
| `optional.indicspecies.min_n` | int | `2` | Minimum accepted n for the indicspecies module. |
| `optional.indicspecies.focus_group1_label` | str | `Bronchial Brush` | Value selecting or naming focus group1 label for the indicspecies module. |
| `optional.indicspecies.label_focused_asvs` | bool | `true` | Enable or disable label focused ASVs behavior in the indicspecies module. |
| `optional.indicspecies.plot_enabled` | bool | `true` | Enable or disable plot enabled behavior in the indicspecies module. |
| `optional.indicspecies.plot_pairs_mode` | str | `all` | Value selecting or naming plot pairs mode for the indicspecies module. |
| `optional.indicspecies.plot_output_dir` | str | `indicspecies/plots` | Directory or published subdirectory used for plot output dir. |
| `optional.indicspecies.aligned_plot_enabled` | bool | `false` | Enable or disable aligned plot enabled behavior in the indicspecies module. |
| `optional.indicspecies.aligned_plot_output_dir` | str | `indicspecies/aligned` | Directory or published subdirectory used for aligned plot output dir. |
| `optional.indicspecies.aligned_alpha` | float | `0.05` | Statistical or filtering cutoff for aligned alpha in the indicspecies module. |
| `optional.indicspecies.aligned_min_stat` | float | `0.0` | Numeric setting for aligned min stat in the indicspecies module. |
| `optional.indicspecies.aligned_top_n` | int | `25` | Configured number or size for aligned top n in the indicspecies module. |
| `optional.indicspecies.venn` | null | `null` | Value selecting or naming venn for the indicspecies module. |
| `optional.indicspecies.taxonomy` | str | `/abs/path/output_dir/taxonomy/ASV_SILVA_tax.full-length.vsearch.tsv` | Value selecting or naming taxonomy for the indicspecies module. |

### `optional.voc_correlation`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `optional.voc_correlation.enabled` | bool | `false` | Whether this module or nested analysis is scheduled. |
| `optional.voc_correlation.metadata_sample_col` | str | `Sample` | Sample-identifier column in the metadata table. |
| `optional.voc_correlation.type_col` | str | `Type_Group` | Column containing sample-type labels. |
| `optional.voc_correlation.patient_col` | str | `Participant_ID` | Column containing participant/patient identifiers for blocking or pairing. |
| `optional.voc_correlation.case_col` | str | `Case` | Column containing case/control status. |
| `optional.voc_correlation.sample_types` | str | `Bronchial Brush,Lung Brush` | Value selecting or naming sample types for the VOC correlation module. |
| `optional.voc_correlation.voc_table` | str | `/abs/path/ref_db/VOC_table.tsv` | Path to, or configured name of, the VOC table input. |
| `optional.voc_correlation.output_dir` | str | `voc_correlation` | Directory or published subdirectory used for output dir. |
| `optional.voc_correlation.voc_sample_col` | str | `sample` | Input-table column containing VOC sample. |
| `optional.voc_correlation.sample_id_mode` | str | `legacy_patient_pair` | Value selecting or naming sample id mode for the VOC correlation module. |
| `optional.voc_correlation.use_legacy_voc_subset` | bool | `true` | Enable or disable use legacy VOC subset behavior in the VOC correlation module. |
| `optional.voc_correlation.correlation_direction` | str | `positive` | Direction retained in the baseline ASV–VOC correlation outputs: positive, negative, or both. |
| `optional.voc_correlation.isa_correlation_direction` | str | `both` | Direction retained in ISA-focused ASV–VOC outputs. |
| `optional.voc_correlation.isa_brush_groups` | list | `["Bronchial Brush","Lung Brush"]` | Ordered values used for isa brush groups by the VOC correlation module. |
| `optional.voc_correlation.isa_all_type_groups` | list | `["Oral Rinse","BAL","Bronchial Brush"]` | Ordered values used for isa all type groups by the VOC correlation module. |
| `optional.voc_correlation.isa_exclude_all_types_from_brush` | bool | `true` | Enable or disable isa exclude all types from brush behavior in the VOC correlation module. |
| `optional.voc_correlation.isa_min_abs_rho` | float | `0.35` | Minimum absolute Spearman rho required in focused ISA/VOC rows and columns. |
| `optional.voc_correlation.sample_min_abs_z` | float | `2.5` | Minimum absolute sample VOC z-score used for focused sample heatmaps. |
| `optional.voc_correlation.patient_inference` | bool | `true` | Add patient-level relative-abundance permutation correlations and CLR sensitivity analysis; legacy sample correlations remain exploratory. |
| `optional.voc_correlation.patient_permutations` | int | `9999` | Seeded permutations for patient correlations and case-status tests; case-status enumeration is exact when all allocations fit this budget. Use 999 for demonstrations, 9999 or more for analysis. |
| `optional.voc_correlation.patient_seed` | int | `42` | Random seed for patient-level VOC permutation tests. |
| `optional.voc_correlation.patient_min_patients` | int | `6` | Minimum patients with cognate measurements to test an ASV–VOC association (at least 3; default 6 is a feasibility gate, not a power guarantee). |
| `optional.voc_correlation.patient_min_nonzero` | int | `3` | Minimum patients with any detected counts for a candidate ASV; insufficient pairs are reported without p-values. |
| `optional.voc_correlation.clr_pseudocount` | float | `0.5` | Positive count pseudocount added to every supplied ASV before sample-wise centered log-ratios; sensitivity analysis only, default 0.5. |
| `optional.voc_correlation.voc_columns` | list | `[]` | Ordered values used for VOC columns by the VOC correlation module. |

### `optional.power_analysis`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `optional.power_analysis.enabled` | bool | `false` | Whether this module or nested analysis is scheduled. |
| `optional.power_analysis.output_dir` | str | `power_analysis` | Directory or published subdirectory used for output dir. |
| `optional.power_analysis.sample_col` | str | `Sample` | Column containing sample identifiers. |
| `optional.power_analysis.patient_col` | str | `Participant_ID` | Column containing participant/patient identifiers for blocking or pairing. |
| `optional.power_analysis.case_col` | str | `Case` | Column containing case/control status. |
| `optional.power_analysis.type_col` | str | `Type_Group` | Column containing sample-type labels. |
| `optional.power_analysis.sample_sizes_cancer` | str | `6,8,10,15,20,25,30` | Value selecting or naming sample sizes cancer for the power analysis module. |
| `optional.power_analysis.sample_sizes_stype` | str | `10,15,20,25,30,40,50` | Value selecting or naming sample sizes stype for the power analysis module. |
| `optional.power_analysis.n_simulations` | int | `1000` | Configured number or size for n simulations in the power analysis module. |
| `optional.power_analysis.n_perm` | int | `199` | Configured number or size for n perm in the power analysis module. |
| `optional.power_analysis.alpha` | float | `0.05` | Numeric setting for alpha in the power analysis module. |
| `optional.power_analysis.seed` | int | `42` | Random seed used to make stochastic behavior reproducible. |
| `optional.power_analysis.skip_estimate` | bool | `false` | Skip estimate when true. |
| `optional.power_analysis.skip_plot` | bool | `false` | Skip plot when true. |
| `optional.power_analysis.transform` | str | `none` | Transformation applied to input abundance values. |
| `optional.power_analysis.keep_contralateral_in_cancer` | bool | `false` | Retain contralateral in cancer when true. |
| `optional.power_analysis.contralateral_sample_types` | str | `Bronchial Brush,BAL` | Value selecting or naming contralateral sample types for the power analysis module. |

### `optional.taxonomy_patient_aware`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `optional.taxonomy_patient_aware.enabled` | bool | `false` | Whether this module or nested analysis is scheduled. |
| `optional.taxonomy_patient_aware.output_dir` | str | `taxonomy_patient_aware` | Directory or published subdirectory used for output dir. |
| `optional.taxonomy_patient_aware.sample_col` | str | `Sample` | Column containing sample identifiers. |
| `optional.taxonomy_patient_aware.patient_col` | str | `Participant_ID` | Column containing participant/patient identifiers for blocking or pairing. |
| `optional.taxonomy_patient_aware.case_col` | str | `Case` | Column containing case/control status. |
| `optional.taxonomy_patient_aware.type_col` | str | `Type_Group` | Column containing sample-type labels. |
| `optional.taxonomy_patient_aware.count_col` | str | `count` | Column containing ASV counts or transformed count values. |
| `optional.taxonomy_patient_aware.tax_levels` | str | `Phylum,Family` | Value selecting or naming tax levels for the taxonomy patient aware module. |
| `optional.taxonomy_patient_aware.sample_types` | str | `Oral Rinse,BAL,Bronchial Brush` | Value selecting or naming sample types for the taxonomy patient aware module. |
| `optional.taxonomy_patient_aware.min_prevalence` | float | `0.1` | Minimum accepted prevalence for the taxonomy patient aware module. |
| `optional.taxonomy_patient_aware.exclude_contralateral_in_cancer` | bool | `true` | Values or groups excluded according to contralateral in cancer. |
| `optional.taxonomy_patient_aware.contralateral_col` | str | `lung_status` | Input-table column containing contralateral. |
| `optional.taxonomy_patient_aware.cancer_site_col` | str | `Cancer_Site` | Input-table column containing cancer site. |
| `optional.taxonomy_patient_aware.lung_side_col` | str | `lung_code` | Input-table column containing lung side. |
| `optional.taxonomy_patient_aware.contralateral_value` | str | `Contralateral` | Value selecting or naming contralateral value for the taxonomy patient aware module. |
| `optional.taxonomy_patient_aware.contralateral_sample_types` | str | `Bronchial Brush,BAL` | Value selecting or naming contralateral sample types for the taxonomy patient aware module. |
| `optional.taxonomy_patient_aware.skip_omnibus` | bool | `false` | Skip omnibus when true. |
| `optional.taxonomy_patient_aware.transform` | str | `none` | Transformation applied to input abundance values. |
| `optional.taxonomy_patient_aware.alpha` | float | `0.05` | Numeric setting for alpha in the taxonomy patient aware module. |
| `optional.taxonomy_patient_aware.top_n` | int | `12` | Configured number or size for top n in the taxonomy patient aware module. |

### `optional.lung_status_analysis`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `optional.lung_status_analysis.enabled` | bool | `false` | Whether this module or nested analysis is scheduled. |
| `optional.lung_status_analysis.output_dir` | str | `lung_status_analysis` | Directory or published subdirectory used for output dir. |
| `optional.lung_status_analysis.sample_col` | str | `Sample` | Column containing sample identifiers. |
| `optional.lung_status_analysis.type_col` | str | `Type_Group` | Column containing sample-type labels. |
| `optional.lung_status_analysis.sample_types` | str | `Bronchial Brush,BAL` | Value selecting or naming sample types for the lung status analysis module. |
| `optional.lung_status_analysis.case_col` | str | `Case` | Column containing case/control status. |
| `optional.lung_status_analysis.patient_col` | str | `Participant_ID` | Column containing participant/patient identifiers for blocking or pairing. |
| `optional.lung_status_analysis.cancer_site_col` | str | `Cancer_Site` | Input-table column containing cancer site. |
| `optional.lung_status_analysis.lung_code_col` | str | `lung_code` | Input-table column containing lung code. |
| `optional.lung_status_analysis.tumor_side_col` | str | `TumorSide` | Input-table column containing tumor side. |
| `optional.lung_status_analysis.contralateral_col` | str | `Contralateral` | Input-table column containing contralateral. |
| `optional.lung_status_analysis.healthy_col` | str | `Healthy` | Input-table column containing healthy. |
| `optional.lung_status_analysis.lung_status_col` | str | `lung_status` | Input-table column containing lung status. |
| `optional.lung_status_analysis.status_a_value` | str | `TumorSide` | Value selecting or naming status a value for the lung status analysis module. |
| `optional.lung_status_analysis.status_b_value` | str | `Contralateral` | Value selecting or naming status b value for the lung status analysis module. |
| `optional.lung_status_analysis.reference_status_value` | str | `Healthy` | Value selecting or naming reference status value for the lung status analysis module. |
| `optional.lung_status_analysis.permutations` | int | `9999` | Number of permutations used by the statistical test. |
| `optional.lung_status_analysis.seed` | int | `1` | Random seed used to make stochastic behavior reproducible. |

### `optional.clustermaps`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `optional.clustermaps.enabled` | bool | `true` | Whether this module or nested analysis is scheduled. |
| `optional.clustermaps.output_dir` | str | `clustermaps` | Directory or published subdirectory used for output dir. |
| `optional.clustermaps.mito_output_dir` | str | `mito/clustermaps` | Directory or published subdirectory used for mito output dir. |
| `optional.clustermaps.mito_input` | null | `null` | Path to, or configured name of, the mito input input. |
| `optional.clustermaps.isa_file` | str | `/abs/path/output_dir/indicspecies/` | Path to, or configured name of, the isa file input. |
| `optional.clustermaps.run_mito` | bool | `true` | Run mito when true. |
| `optional.clustermaps.sample_col` | str | `Sample` | Column containing sample identifiers. |
| `optional.clustermaps.sample_code_col` | str | `sample_code` | Input-table column containing sample code. |
| `optional.clustermaps.asv_id_col` | str | `ASV_ID` | Input-table column containing ASV id. |
| `optional.clustermaps.group1_col` | str | `Type_Group` | Primary display or analysis grouping column. |
| `optional.clustermaps.group2_col` | str | `Case` | Secondary display or analysis grouping column. |
| `optional.clustermaps.group3_col` | null | `null` | Input-table column containing group3. |
| `optional.clustermaps.group4_col` | null | `null` | Input-table column containing group4. |
| `optional.clustermaps.group1_order` | list | `["Oral Rinse","BAL+Oral Rinse","BAL","BAL+Bronchial Brush","Bronchial Brush","Bronchial Brush+Oral Rinse","BAL+Bronchial Brush+Oral Rinse","not_indicator"]` | Explicit category order used for group1 analysis and display. |
| `optional.clustermaps.exclude_group1` | null | `null` | Values or groups excluded according to group1. |
| `optional.clustermaps.group1_palette` | str | `Oral Rinse=#6A3D9A,BAL+Oral Rinse=#E78AC3,BAL=#0072B2,BAL+Bronchial Brush=#5CC8C8,Bronchial Brush=#009E73,Bronchial Brush+Oral Rinse=#B8E186,BAL+Bronchial Brush+Oral Rinse=#CBB6E9,not_indicator=#D3D3D3` | Explicit label-to-color mapping used for group1 displays. |
| `optional.clustermaps.group2_palette` | str | `Non-Cancer=#FFFFFF,Cancer=#A50026,Cancer+Non-Cancer=#000000,not_indicator=#D3D3D3` | Explicit label-to-color mapping used for group2 displays. |
| `optional.clustermaps.group3_palette` | null | `null` | Explicit label-to-color mapping used for group3 displays. |
| `optional.clustermaps.group4_palette` | null | `null` | Explicit label-to-color mapping used for group4 displays. |
| `optional.clustermaps.ranks` | str | `Phylum,Class,Order,Family,Genus,Species,ASV_ID` | Value selecting or naming ranks for the clustermaps module. |
| `optional.clustermaps.topN` | str | `Phylum=30,Class=30,Order=30,Family=30,Genus=30,Species=30,ASV_ID=200` | Value selecting or naming topN for the clustermaps module. |
| `optional.clustermaps.count_col` | str | `corr_count` | Column containing ASV counts or transformed count values. |
| `optional.clustermaps.isa_min_stat` | float | `0.6` | Numeric setting for isa min stat in the clustermaps module. |
| `optional.clustermaps.isa_significance_cols` | null | `null` | Input-table columns used for isa significance. |
| `optional.clustermaps.isa_stat_cols` | null | `null` | Input-table columns used for isa stat. |
| `optional.clustermaps.formats` | str | `pdf,png,svg` | Comma-separated output figure formats. |
| `optional.clustermaps.figwidth` | int | `32` | Numeric setting for figwidth in the clustermaps module. |
| `optional.clustermaps.row_height` | float | `0.08` | Numeric setting for row height in the clustermaps module. |
| `optional.clustermaps.min_height` | int | `8` | Minimum accepted height for the clustermaps module. |
| `optional.clustermaps.max_height` | int | `180` | Maximum accepted height for the clustermaps module. |
| `optional.clustermaps.mito_sample_mode` | str | `auto` | Value selecting or naming mito sample mode for the clustermaps module. |

### `optional.spieceasi`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `optional.spieceasi.enabled` | bool | `true` | Whether this module or nested analysis is scheduled. |
| `optional.spieceasi.output_dir` | str | `spieceasi` | Directory or published subdirectory used for output dir. |
| `optional.spieceasi.prefix` | str | `spieceasi` | Value selecting or naming prefix for the spieceasi module. |
| `optional.spieceasi.network_enabled` | bool | `true` | Enable or disable network enabled behavior in the spieceasi module. |
| `optional.spieceasi.modules_enabled` | bool | `true` | Enable or disable modules enabled behavior in the spieceasi module. |
| `optional.spieceasi.module_methods` | str | `leiden` | Value selecting or naming module methods for the spieceasi module. |
| `optional.spieceasi.module_primary_method` | str | `leiden` | Value selecting or naming module primary method for the spieceasi module. |
| `optional.spieceasi.module_resolutions` | str | `0.8,1.0,1.2` | Value selecting or naming module resolutions for the spieceasi module. |
| `optional.spieceasi.module_reps` | int | `50` | Numeric setting for module reps in the spieceasi module. |
| `optional.spieceasi.module_consensus_threshold` | float | `0.95` | Statistical or filtering cutoff for module consensus threshold in the spieceasi module. |
| `optional.spieceasi.module_seed` | int | `42` | Numeric setting for module seed in the spieceasi module. |
| `optional.spieceasi.module_best_only` | bool | `true` | Enable or disable module best only behavior in the spieceasi module. |
| `optional.spieceasi.module_best_min_size` | int | `5` | Configured number or size for module best min size in the spieceasi module. |
| `optional.spieceasi.module_best_min_stability` | float | `0.7` | Numeric setting for module best min stability in the spieceasi module. |
| `optional.spieceasi.module_isa_only` | bool | `false` | Enable or disable module isa only behavior in the spieceasi module. |
| `optional.spieceasi.module_color_by_isa` | bool | `false` | Enable or disable module color by isa behavior in the spieceasi module. |
| `optional.spieceasi.module_isa_source` | str | `group1` | Value selecting or naming module isa source for the spieceasi module. |
| `optional.spieceasi.module_isa_min_stat` | float | `0.25` | Numeric setting for module isa min stat in the spieceasi module. |
| `optional.spieceasi.module_isa_max_q` | float | `0.05` | Numeric setting for module isa max q in the spieceasi module. |
| `optional.spieceasi.modules_sub` | null | `null` | Value selecting or naming modules sub for the spieceasi module. |
| `optional.spieceasi.modules_all` | null | `null` | Value selecting or naming modules all for the spieceasi module. |
| `optional.spieceasi.graph_pos_all` | null | `null` | Value selecting or naming graph pos all for the spieceasi module. |
| `optional.spieceasi.graph_pos_sub` | null | `null` | Value selecting or naming graph pos sub for the spieceasi module. |
| `optional.spieceasi.node_features` | null | `null` | Value selecting or naming node features for the spieceasi module. |
| `optional.spieceasi.isa_overlay_groups` | null | `null` | Value selecting or naming isa overlay groups for the spieceasi module. |
| `optional.spieceasi.group_palettes.Type_Group` | str | `Oral Rinse=#6A3D9A,BAL+Oral Rinse=#E78AC3,BAL=#0072B2,BAL+Bronchial Brush=#5CC8C8,Bronchial Brush=#009E73,Bronchial Brush+Oral Rinse=#B8E186,BAL+Bronchial Brush+Oral Rinse=#CBB6E9,not_indicator=#D3D3D3` | Value selecting or naming Type Group for the spieceasi module. |
| `optional.spieceasi.group_palettes.Case` | str | `Non-Cancer=#FFFFFF,Cancer=#A50026,Cancer+Non-Cancer=#000000,not_indicator=#D3D3D3` | Value selecting or naming Case for the spieceasi module. |
| `optional.spieceasi.group_orders.Type_Group` | str | `Oral Rinse,BAL+Oral Rinse,BAL,BAL+Bronchial Brush,Bronchial Brush,Bronchial Brush+Oral Rinse,BAL+Bronchial Brush+Oral Rinse,not_indicator` | Value selecting or naming Type Group for the spieceasi module. |
| `optional.spieceasi.group_orders.Case` | str | `Non-Cancer,Cancer,Cancer+Non-Cancer` | Value selecting or naming Case for the spieceasi module. |
| `optional.spieceasi.focus_labels` | mapping | `{}` | Value selecting or naming focus labels for the spieceasi module. |
| `optional.spieceasi.network_modes` | list | `["all"]` | Ordered values used for network modes by the spieceasi module. |
| `optional.spieceasi.metadata` | str | `/abs/path/metadata.tsv` | Path to the metadata TSV consumed by this module. |
| `optional.spieceasi.color_col` | str | `Color` | Metadata column containing or selecting display colors. |
| `optional.spieceasi.group1_order` | list | `["Oral Rinse","BAL+Oral Rinse","BAL","BAL+Bronchial Brush","Bronchial Brush","Bronchial Brush+Oral Rinse","BAL+Bronchial Brush+Oral Rinse","not_indicator"]` | Explicit category order used for group1 analysis and display. |
| `optional.spieceasi.group2_order` | list | `["Non-Cancer","Cancer","Cancer+Non-Cancer"]` | Explicit category order used for group2 analysis and display. |
| `optional.spieceasi.group1_palette` | str | `Oral Rinse=#6A3D9A,BAL+Oral Rinse=#E78AC3,BAL=#0072B2,BAL+Bronchial Brush=#5CC8C8,Bronchial Brush=#009E73,Bronchial Brush+Oral Rinse=#B8E186,BAL+Bronchial Brush+Oral Rinse=#CBB6E9,not_indicator=#D3D3D3` | Explicit label-to-color mapping used for group1 displays. |
| `optional.spieceasi.group2_palette` | str | `Non-Cancer=#FFFFFF,Cancer=#A50026,Cancer+Non-Cancer=#000000,not_indicator=#D3D3D3` | Explicit label-to-color mapping used for group2 displays. |
| `optional.spieceasi.focus_group1_label` | str | `Bronchial Brush` | Value selecting or naming focus group1 label for the spieceasi module. |
| `optional.spieceasi.transpose` | bool | `true` | Transpose the configured matrix orientation before analysis. |
| `optional.spieceasi.min_rel_abund` | float | `0.005` | Minimum relative abundance required for an ASV to enter network inference. |
| `optional.spieceasi.min_prevalence` | float | `0.05` | Minimum prevalence required for an ASV to enter network inference. |
| `optional.spieceasi.remove_zero_var` | bool | `true` | Enable or disable remove zero var behavior in the spieceasi module. |
| `optional.spieceasi.force_keep_indicator_asvs` | bool | `true` | Force keep indicator ASVs instead of accepting a reusable cached product. |
| `optional.spieceasi.method` | str | `glasso` | Value selecting or naming method for the spieceasi module. |
| `optional.spieceasi.lambda_min_ratio` | float | `0.01` | Numeric setting for lambda min ratio in the spieceasi module. |
| `optional.spieceasi.nlambda` | int | `20` | Numeric setting for nlambda in the spieceasi module. |
| `optional.spieceasi.rep_num` | int | `50` | Numeric setting for rep num in the spieceasi module. |
| `optional.spieceasi.thresh` | float | `0.1` | Numeric setting for thresh in the spieceasi module. |
| `optional.spieceasi.ncores` | int | `8` | Parallel worker count requested by this module. |
| `optional.spieceasi.seed` | int | `10010` | Random seed used to make stochastic behavior reproducible. |
| `optional.spieceasi.edge_threshold` | float | `0.05` | Absolute association-weight threshold used to retain reported network edges. |
| `optional.spieceasi.keep_negative` | bool | `true` | Retain negative as well as positive inferred network edges. |
| `optional.spieceasi.layout_iters` | int | `1000` | Numeric setting for layout iters in the spieceasi module. |
| `optional.spieceasi.force_filter` | bool | `true` | Force filter instead of accepting a reusable cached product. |
| `optional.spieceasi.force_spieceasi` | bool | `true` | Force spieceasi instead of accepting a reusable cached product. |
| `optional.spieceasi.force_graphs` | bool | `true` | Force graphs instead of accepting a reusable cached product. |
| `optional.spieceasi.layout_seed` | int | `42` | Numeric setting for layout seed in the spieceasi module. |
| `optional.spieceasi.layout_scale` | float | `3.0` | Numeric setting for layout scale in the spieceasi module. |
| `optional.spieceasi.degree_scale` | float | `80.0` | Numeric setting for degree scale in the spieceasi module. |
| `optional.spieceasi.degree_size_mode` | str | `legacy` | Value selecting or naming degree size mode for the spieceasi module. |
| `optional.spieceasi.degree_min_area` | float | `0.0` | Numeric setting for degree min area in the spieceasi module. |
| `optional.spieceasi.edge_width_scale` | float | `10.0` | Numeric setting for edge width scale in the spieceasi module. |
| `optional.spieceasi.isa_scale` | float | `700.0` | Numeric setting for isa scale in the spieceasi module. |
| `optional.spieceasi.abundance_size_mode` | str | `legacy` | Value selecting or naming abundance size mode for the spieceasi module. |
| `optional.spieceasi.abundance_reference` | float | `5000.0` | Numeric setting for abundance reference in the spieceasi module. |
| `optional.spieceasi.abundance_reference_area` | float | `80.0` | Numeric setting for abundance reference area in the spieceasi module. |
| `optional.spieceasi.abundance_min_area` | float | `8.0` | Numeric setting for abundance min area in the spieceasi module. |
| `optional.spieceasi.abundance_max_area` | float | `420.0` | Numeric setting for abundance max area in the spieceasi module. |
| `optional.spieceasi.abundance_scale_power` | float | `1.6` | Numeric setting for abundance scale power in the spieceasi module. |
| `optional.spieceasi.pulsar_criterion` | str | `stars` | Model selection: stars uses ordinary StARS; bstars requests bounded StARS via lower/upper bounds. The template preserves ordinary StARS. |

### `optional.network_topology`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `optional.network_topology.enabled` | bool | `false` | Whether this module or nested analysis is scheduled. |
| `optional.network_topology.output_dir` | str | `spieceasi` | Directory or published subdirectory used for output dir. |
| `optional.network_topology.n_null` | int | `1000` | Number of seeded degree-preserving null-network draws. |
| `optional.network_topology.seed` | int | `42` | Random seed used to make stochastic behavior reproducible. |
| `optional.network_topology.skip_null` | bool | `false` | Report observed topology without generating the null ensemble. |

### `optional.asv_mag_link`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `optional.asv_mag_link.enabled` | bool | `false` | Whether this module or nested analysis is scheduled. |
| `optional.asv_mag_link.genome_qc_dir` | null | `null` | Directory or published subdirectory used for genome qc dir. |
| `optional.asv_mag_link.genome_qc_dirs` | null | `null` | Value selecting or naming genome qc dirs for the ASV mag link module. |
| `optional.asv_mag_link.id_token_indexes` | null | `null` | Value selecting or naming id token indexes for the ASV mag link module. |
| `optional.asv_mag_link.barrnap_dir` | null | `null` | Directory or published subdirectory used for barrnap dir. |
| `optional.asv_mag_link.genome_fasta_dir` | null | `null` | Directory or published subdirectory used for genome fasta dir. |
| `optional.asv_mag_link.output_dir` | str | `asv_mag_link` | Directory or published subdirectory used for output dir. |
| `optional.asv_mag_link.threads` | int | `8` | CPU threads requested by this module. |
| `optional.asv_mag_link.min_pident` | float | `97.0` | Minimum ASV-to-SSU percent identity accepted as a MAG link. |
| `optional.asv_mag_link.min_qcov` | float | `90.0` | Minimum ASV query coverage accepted as a MAG link. |
| `optional.asv_mag_link.top_n` | int | `5` | Configured number or size for top n in the ASV mag link module. |
| `optional.asv_mag_link.plot_top_n` | int | `20` | Configured number or size for plot top n in the ASV mag link module. |
| `optional.asv_mag_link.master_tsv` | null | `null` | Optional precomputed ASV-MAG master-link TSV, used instead of recomputing links from genome/Barrnap inputs. |

### `optional.master_summary`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `optional.master_summary.enabled` | bool | `true` | Whether this module or nested analysis is scheduled. |
| `optional.master_summary.output_dir` | str | `summary/tables` | Directory or published subdirectory used for output dir. |
| `optional.master_summary.clustermaps_dir` | str | `clustermaps` | Directory or published subdirectory used for clustermaps dir. |
| `optional.master_summary.indicspecies_dir` | str | `indicspecies` | Directory or published subdirectory used for indicspecies dir. |
| `optional.master_summary.spieceasi_dir` | str | `spieceasi` | Directory or published subdirectory used for spieceasi dir. |
| `optional.master_summary.asv_mag_dir` | str | `asv_mag_link` | Directory or published subdirectory used for ASV mag dir. |
| `optional.master_summary.max_direct_cols` | int | `300` | Input-table columns used for max direct. |
| `optional.master_summary.whitelist` | list | `["type_group_ISA_enriched.tsv","case_ISA_enriched.tsv","type_group_case_ISA_results.tsv","type_group_indicator_species_summary.tsv","Case_indicator_species_summary.tsv","type_group_indicator_species_results.tsv","Case_indicator_species_results.tsv","spieceasi_node_features.csv","spieceasi_modules_sub.tsv","spieceasi_modules_all.tsv","module_asv_anchor_table.tsv","clustermap_ASV_ID_plot.tsv","asv2mag_pairing.tsv","asv2mag_summary.tsv","asv2mag_genome_summary.tsv"]` | Ordered values used for whitelist by the master summary module. |

### `optional.measurement_association`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `optional.measurement_association.enabled` | bool | `false` | Whether this module or nested analysis is scheduled. |
| `optional.measurement_association.output_dir` | str | `measurement_association` | Directory or published subdirectory used for output dir. |
| `optional.measurement_association.measurement_table` | null | `null` | Optional sample-level measurement table; when absent, eligible numeric metadata columns supply measurements. |
| `optional.measurement_association.sample_col` | str | `Sample` | Column containing sample identifiers. |
| `optional.measurement_association.asv_id_col` | str | `ASV_ID` | Input-table column containing ASV id. |
| `optional.measurement_association.measurement_sample_col` | str | `Sample` | Input-table column containing measurement sample. |
| `optional.measurement_association.metadata_join_cols` | list | `[]` | Metadata join keys paired positionally with measurement_join_cols; empty lists use the configured sample identifiers. |
| `optional.measurement_association.measurement_join_cols` | list | `[]` | Measurement-table keys paired positionally with metadata_join_cols. |
| `optional.measurement_association.measurement_cols` | list | `[]` | Input-table columns used for measurement. |
| `optional.measurement_association.exclude_cols` | list | `[]` | Input-table columns used for exclude. |
| `optional.measurement_association.group_col` | str | `Type_Group` | Primary grouping column used by this module. |
| `optional.measurement_association.group_palette` | null | `null` | Explicit label-to-color mapping used for group displays. |
| `optional.measurement_association.max_asvs` | int | `300` | Maximum accepted ASVs for the measurement association module. |
| `optional.measurement_association.min_total` | float | `0.0` | Minimum accepted total for the measurement association module. |
| `optional.measurement_association.min_prevalence` | float | `0.0` | Minimum accepted prevalence for the measurement association module. |
| `optional.measurement_association.top_correlations` | int | `100` | Numeric setting for top correlations in the measurement association module. |
| `optional.measurement_association.correlation_direction` | str | `both` | Value selecting or naming correlation direction for the measurement association module. |
| `optional.measurement_association.ordination_methods` | str | `cca,rda,dbrda` | Comma-separated constrained ordination methods to run: cca, rda, and/or dbrda. |
| `optional.measurement_association.permutations` | int | `999` | Number of permutations used by the statistical test. |
| `optional.measurement_association.top_vectors` | int | `12` | Numeric setting for top vectors in the measurement association module. |
| `optional.measurement_association.formats` | str | `pdf,png,svg` | Comma-separated output figure formats. |

### `optional.grouping_diagnostics`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `optional.grouping_diagnostics.enabled` | bool | `false` | Whether this module or nested analysis is scheduled. |
| `optional.grouping_diagnostics.output_dir` | str | `grouping_diagnostics` | Directory or published subdirectory used for output dir. |
| `optional.grouping_diagnostics.sample_col` | str | `Sample` | Column containing sample identifiers. |
| `optional.grouping_diagnostics.group_cols` | list | `["Type_Group"]` | Input-table columns used for group. |
| `optional.grouping_diagnostics.baseline_group` | null | `null` | Value selecting or naming baseline group for the grouping diagnostics module. |
| `optional.grouping_diagnostics.primary_group` | str | `Type_Group` | Value selecting or naming primary group for the grouping diagnostics module. |
| `optional.grouping_diagnostics.group_palettes` | mapping | `{}` | Explicit label-to-color mapping used for group displays. |
| `optional.grouping_diagnostics.group_orders` | mapping | `{}` | Explicit category order used for group analysis and display. |
| `optional.grouping_diagnostics.distance_metrics` | list | `["bray"]` | Ordered values used for distance metrics by the grouping diagnostics module. |
| `optional.grouping_diagnostics.transform` | str | `relative` | Transformation applied to input abundance values. |
| `optional.grouping_diagnostics.permutations` | int | `999` | Number of permutations used by the statistical test. |
| `optional.grouping_diagnostics.random_state` | int | `42` | Random seed used to make stochastic behavior reproducible. |
| `optional.grouping_diagnostics.formats` | str | `pdf,png,svg` | Comma-separated output figure formats. |
| `optional.grouping_diagnostics.soft_labeling.enabled` | bool | `false` | Whether this module or nested analysis is scheduled. |
| `optional.grouping_diagnostics.soft_labeling.k` | int | `7` | Numeric setting for k in the grouping diagnostics module. |
| `optional.grouping_diagnostics.soft_labeling.target_cols` | list | `["Type_Group"]` | Input-table columns used for target. |
| `optional.grouping_diagnostics.soft_labeling.exclude_labels` | list | `["outlier"]` | Labels excluded from soft-label training and predictions; observed excluded labels are retained rather than overwritten. |
| `optional.grouping_diagnostics.soft_labeling.min_class_samples` | int | `3` | Minimum accepted class samples for the grouping diagnostics module. |
| `optional.grouping_diagnostics.soft_labeling.distance_quantile` | float | `0.95` | Numeric setting for distance quantile in the grouping diagnostics module. |
| `optional.grouping_diagnostics.soft_labeling.apply_downstream` | bool | `false` | Apply validated proposed labels to downstream metadata and long ASV annotations; false keeps diagnostic proposals separate from analysis labels. |
| `optional.grouping_diagnostics.soft_labeling.target_col` | str | `Type_Group` | Input-table column containing target. |
| `optional.grouping_diagnostics.soft_labeling.min_confidence` | float | `0.7` | Minimum accepted confidence for the grouping diagnostics module. |
| `optional.grouping_diagnostics.soft_labeling.min_neighbor_agreement` | float | `0.6` | Minimum accepted neighbor agreement for the grouping diagnostics module. |
| `optional.grouping_diagnostics.soft_labeling.min_cv_balanced_accuracy` | float | `0.6` | Minimum cross-validated balanced accuracy required before proposed labels can be applied downstream. |
| `optional.grouping_diagnostics.power.enabled` | bool | `false` | Whether this module or nested analysis is scheduled. |
| `optional.grouping_diagnostics.power.sample_sizes` | list | `[3,5,10,15,20]` | Ordered values used for sample sizes by the grouping diagnostics module. |
| `optional.grouping_diagnostics.power.simulations` | int | `100` | Numeric setting for simulations in the grouping diagnostics module. |
| `optional.grouping_diagnostics.power.permutations` | int | `99` | Number of permutations used by the statistical test. |
| `optional.grouping_diagnostics.power.alpha` | float | `0.05` | Numeric setting for alpha in the grouping diagnostics module. |
| `optional.grouping_diagnostics.power.min_groups` | int | `2` | Minimum accepted groups for the grouping diagnostics module. |

### `optional.asv_mag_network`

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `optional.asv_mag_network.enabled` | bool | `false` | Whether this module or nested analysis is scheduled. |
| `optional.asv_mag_network.output_dir` | str | `asv_mag_network` | Directory or published subdirectory used for output dir. |
| `optional.asv_mag_network.prefix` | str | `asv_mag_network` | Value selecting or naming prefix for the ASV mag network module. |
| `optional.asv_mag_network.graph_variant` | str | `all` | Network graph to annotate: all or the thresholded subgraph. |
| `optional.asv_mag_network.min_pident` | float | `99.5` | Minimum accepted pident for the ASV mag network module. |
| `optional.asv_mag_network.min_qcov` | float | `100.0` | Minimum accepted qcov for the ASV mag network module. |
| `optional.asv_mag_network.asv_taxonomy_source` | str | `ncbi` | Value selecting or naming ASV taxonomy source for the ASV mag network module. |
| `optional.asv_mag_network.mag_taxonomy_source` | str | `gtdb` | Value selecting or naming mag taxonomy source for the ASV mag network module. |
| `optional.asv_mag_network.mag_id_mode` | str | `exact` | Value selecting or naming mag id mode for the ASV mag network module. |
| `optional.asv_mag_network.mag_abundance` | null | `null` | Value selecting or naming mag abundance for the ASV mag network module. |
| `optional.asv_mag_network.mag_abundance_format` | str | `long` | Value selecting or naming mag abundance format for the ASV mag network module. |
| `optional.asv_mag_network.mag_abundance_genome_col` | str | `genome_id` | Input-table column containing mag abundance genome. |
| `optional.asv_mag_network.mag_abundance_sample_col` | str | `sample_id` | Input-table column containing mag abundance sample. |
| `optional.asv_mag_network.mag_abundance_value_col` | str | `read_count` | Input-table column containing mag abundance value. |
| `optional.asv_mag_network.min_shared_samples` | int | `5` | Minimum accepted shared samples for the ASV mag network module. |
| `optional.asv_mag_network.abundance_transform` | str | `log1p` | Value selecting or naming abundance transform for the ASV mag network module. |
| `optional.asv_mag_network.functional_module_min_fraction` | float | `0.5` | Numeric setting for functional module min fraction in the ASV mag network module. |
| `optional.asv_mag_network.functional_annotations` | list | `[]` | Optional functional-annotation tables to summarize across linked MAGs and network modules. |

## `environments`

### Process environment definitions

| Parameter | Type | Template value | Definition |
|---|---|---|---|
| `environments.main` | str | `processes/fastp_qc/env.yml` | Conda environment YAML used by the main process family. |
| `environments.sina` | str | `processes/sina_trim/env.yml` | Conda environment YAML used by the sina process family. |
| `environments.taxonomy` | str | `processes/taxonomy/env.yml` | Conda environment YAML used by the taxonomy process family. |
| `environments.mitomaster` | str | `processes/mitomaster/env.yml` | Conda environment YAML used by the mitomaster process family. |
| `environments.mito_checker` | str | `processes/mito_decontam/env.yml` | Conda environment YAML used by the mito checker process family. |
| `environments.filter_counts` | str | `processes/shared_envs/asv_pipeline.yml` | Conda environment YAML used by the filter counts process family. |
| `environments.three_tier_decontam` | str | `processes/three_tier_decontam/env.yml` | Conda environment YAML used by the three tier decontam process family. |
| `environments.general_stats` | str | `processes/general_stats/env.yml` | Conda environment YAML used by the general stats process family. |
| `environments.sankey` | str | `processes/sankey/env.yml` | Conda environment YAML used by the sankey process family. |
| `environments.plot_metadata` | str | `processes/plot_metadata/env.yml` | Conda environment YAML used by the plot metadata process family. |
| `environments.batch_correction` | str | `processes/asv_batch_correction/env.yml` | Conda environment YAML used by the batch correction process family. |
| `environments.outlier_checker` | str | `processes/outlier_checker/env.yml` | Conda environment YAML used by the outlier checker process family. |
| `environments.collectors_curve` | str | `processes/collectors_curve/env.yml` | Conda environment YAML used by the collectors curve process family. |
| `environments.plot_upset` | str | `processes/plot_upset/env.yml` | Conda environment YAML used by the plot upset process family. |
| `environments.bubbleplotter` | str | `processes/bubbleplotter/env.yml` | Conda environment YAML used by the bubbleplotter process family. |
| `environments.umap_clustering` | str | `processes/umap_clustering/env.yml` | Conda environment YAML used by the umap clustering process family. |
| `environments.diversity` | str | `processes/diversity_analysis/env.yml` | Conda environment YAML used by the diversity process family. |
| `environments.indicspecies` | str | `processes/indicspecies/env.yml` | Conda environment YAML used by the indicspecies process family. |
| `environments.voc_correlation` | str | `processes/voc_correlation/env.yml` | Conda environment YAML used by the VOC correlation process family. |
| `environments.clustermaps` | str | `processes/clustermaps/env.yml` | Conda environment YAML used by the clustermaps process family. |
| `environments.power_analysis` | str | `processes/power_analysis_pipeline/env.yml` | Conda environment YAML used by the power analysis process family. |
| `environments.taxonomy_patient_aware` | str | `processes/taxonomy_patient_aware/env.yml` | Conda environment YAML used by the taxonomy patient aware process family. |
| `environments.lung_status_analysis` | str | `processes/lung_status_analysis/env.yml` | Conda environment YAML used by the lung status analysis process family. |
| `environments.spieceasi` | str | `processes/spieceasi/env.yml` | Conda environment YAML used by the spieceasi process family. |
| `environments.network_modules` | str | `processes/network_modules/env.yml` | Conda environment YAML used by the network modules process family. |
| `environments.network` | str | `processes/graph_network/env.yml` | Conda environment YAML used by the network process family. |
| `environments.network_topology` | str | `processes/network_topology/env.yml` | Conda environment YAML used by the network topology process family. |
| `environments.master_summary` | str | `processes/master_summary/env.yml` | Conda environment YAML used by the master summary process family. |
| `environments.asv_mag_link` | str | `processes/asv_mag_link/env.yml` | Conda environment YAML used by the ASV mag link process family. |
| `environments.measurement_association` | str | `processes/measurement_association/env.yml` | Conda environment YAML used by the measurement association process family. |
| `environments.grouping_diagnostics` | str | `processes/grouping_diagnostics/env.yml` | Conda environment YAML used by the grouping diagnostics process family. |
| `environments.group_label_augmentation` | str | `processes/group_label_augmentation/env.yml` | Conda environment YAML used by the group label augmentation process family. |
| `environments.asv_mag_network` | str | `processes/asv_mag_network/env.yml` | Conda environment YAML used by the ASV mag network process family. |
