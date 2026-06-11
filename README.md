# ASPIRE

ASPIRE is the current Nextflow-based ASV analysis workflow in this repository.

The primary entrypoint is `run_asv_pipeline.sh`, which wraps `asv_pipeline.nf` and manages:
- controller environment bootstrapping
- resume behavior
- stage-aware reruns with `--rerun-from`

Key files:
- `run_asv_pipeline.sh`: main wrapper used for runs
- `asv_pipeline.nf`: current workflow definition
- `asv_pipeline_nextflow.yml`: canonical config template
- `examples/`: copy-ready starter configs and input templates

# ASPIRE Nextflow Pipeline

## TL;DR Quick Start

This is the shortest reliable path to run the pipeline on a new dataset.

### 1. Install runtime prerequisites

```bash
# one-time: mamba is required by the wrapper
command -v mamba
```

Notes:
- Main entrypoint is `run_asv_pipeline.sh`, which bootstraps a controller env at `ASPIRE/.controller_env` from `processes/controller/env.yml` (contains `nextflow` + `yq`).
- `run_asv_pipeline.sh` reads `paths.work_dir` and `paths.conda_cache_dir` from your YAML and exports `NXF_WORK` / `NXF_CONDA_CACHEDIR`.

### 2. Copy and edit a pipeline YAML template

```bash
cd ASPIRE
# Choose one:
cp examples/my_run.full.yml my_run.yml   # metadata + mito workflow
# or
cp examples/my_run.min.yml my_run.yml    # core-only starter
```

You must edit at least:
- `paths.input_dir`
- `paths.output_dir`
- `paths.manifest` (recommended; required for default metadata stage)
- Any placeholder paths under `/abs/path/...`

You can still use `asv_pipeline_nextflow.yml` directly if you prefer starting from the full canonical config.

### 3. Prepare required inputs

#### 3a) FASTQ inputs
- Put FASTQs under `paths.input_dir`, or reference them in `paths.manifest`.
- Supported extensions: `.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`.

#### 3b) Manifest (recommended; required by default metadata plotting)
- Tab-separated, **no header**.
- Columns:
  1. `sample_id`
  2. `fastq_r1`
  3. `fastq_r2` (optional for single-end)
- `#` comment lines are allowed.
- Relative file paths are resolved relative to the manifest file directory.
- Starter template: `examples/manifest.example.tsv`

Example:

```tsv
#sample_id	fastq_r1	fastq_r2
S01	/path/reads/S01_R1.fastq.gz	/path/reads/S01_R2.fastq.gz
S02	/path/reads/S02_R1.fastq.gz	/path/reads/S02_R2.fastq.gz
SE03	/path/reads/SE03_R1.fastq.gz
```

#### 3c) Metadata (required for default config sections)
With the shipped `asv_pipeline_nextflow.yml`, these sections are enabled and require metadata:
- `filter_counts.enabled: true` (metadata optional, but used when provided)
- `metadata_plots.enabled: true` (metadata required)
- `sankey.enabled: true` (metadata required)

Minimal columns for the default example settings:
- `sampleID` (must match manifest `sample_id`)
- `Depth`
- `Color`

Starter template:
- `examples/metadata.example.tsv`

### 4. Set custom BLAST databases (required if `mito.enabled: true`)

Default config has `mito.enabled: true`, so set:
- `mito.mito_db`
- `mito.biof_db`

These should point to BLAST DB base paths (or compatible fasta path). Recommended:

```bash
makeblastdb -in mito_ncbi.fasta -dbtype nucl -out /path/db/mito_ncbi
makeblastdb -in contaminants.fasta -dbtype nucl -out /path/db/ssu_pipeline_contaminants
```

Then in YAML:

```yaml
mito:
  enabled: true
  mito_db: /path/db/mito_ncbi
  biof_db: /path/db/ssu_pipeline_contaminants
```

### 5. (Optional) Core-only mode for first test run

If you want to run only the core ASV+taxonomy flow without metadata/mito extras:

```yaml
mito:
  enabled: false
filter_counts:
  enabled: false
metadata_plots:
  enabled: false
sankey:
  enabled: false
batch_correction:
  enabled: false
outlier_detection:
  enabled: false
collectors_curve:
  enabled: false
diversity:
  enabled: false
indicspecies:
  enabled: false
clustermaps:
  enabled: false
spieceasi:
  enabled: false
```

### 6. Validate, then run (main entrypoint)

```bash
# syntax/help check
./run_asv_pipeline.sh --help

# real run
./run_asv_pipeline.sh my_run.yml

# resume
./run_asv_pipeline.sh my_run.yml

# force rerun from a stage onward
./run_asv_pipeline.sh my_run.yml --rerun-from FILTER_COUNTS
```

Direct Nextflow invocation is still supported:

```bash
nextflow run asv_pipeline.nf -params-file my_run.yml
```

If using `-params-file`, set `config_root` in YAML when you need deterministic relative-path resolution.

---

## Provided Templates

Copy-ready files are available in `ASPIRE/examples/`:

- `examples/my_run.min.yml`: core ASV + taxonomy only.
- `examples/my_run.full.yml`: metadata + mito workflow.
- `examples/my_run.patient_aware_test.yml`: minimal test config for the new patient-aware branches.
- `examples/my_run.asv_mag_test.yml`: minimal test config for the ASV-to-genome barrnap linkage branch.
- `examples/manifest.example.tsv`: manifest format template.
- `examples/metadata.example.tsv`: minimal metadata schema template.

Template usage and field descriptions are documented in:
- `examples/README.md`

---

## What This Pipeline Runs

`asv_pipeline.nf` executes a core ASV workflow plus optional downstream modules, all controlled by `asv_pipeline_nextflow.yml`.

### Core ASV path
1. `FASTP_QC`
2. `MERGE_READS`
3. `FILTER_READS`
4. `RELABEL_FILTERED` (if `concat.relabel`)
5. `DEREPLICATE`
6. `DENOISE`
7. `CHIMERA_CHECK`
8. `CREATE_COUNT_MATRIX`
9. `FILTER_TABLE`
10. `SINA_TRIM`
11. `TAXONOMY`

### Optional modules (YAML-controlled)
- Mito decontamination and count filtering: `MITOMASTER`, `MITO_DECONTAM`, `FILTER_COUNTS`
- Metadata/reporting: `GENERAL_STATS`, `PLOT_METADATA`, `SANKEY`, `PLOT_UPSET`, `BUBBLEPLOTTER`, `UMAP_CLUSTERING`
- Batch/outlier/diversity/ecology/network: `ASV_BATCH_CORRECTION`, `OUTLIER_CHECKER`, `COLLECTORS_CURVE`, `DIVERSITY_ANALYSIS`, `INDICSPECIES`, `CLUSTERMAPS`, `SPIECEASI`, `NETWORK_MODULES`, `GRAPH_NETWORK`

---

## Input Requirements

### Always required
- `paths.input_dir`
- `paths.output_dir`
- Valid FASTQ files (discovered from manifest or filename patterns)

### Required depending on enabled sections

- `paths.manifest`
  - Strongly recommended.
  - Effectively required when `metadata_plots.enabled: true` (default in example YAML), because metadata plotting consumes `--sample-manifest`.

- `filter_counts.metadata` (if `filter_counts.enabled` and you want group/sample filtering)
  - Defaults in example YAML use this file.

- `metadata_plots.metadata` (if `metadata_plots.enabled`)
  - Must contain configured columns (`sample_col`, `type_col`, `color_col`; defaults often `sampleID`, `Depth`, `Color`).

- `sankey.metadata` (if `sankey.enabled`)
  - Also needs `filter_counts.enabled: true`, `filter_counts.save_intermediates: true`, and `general_stats.enabled: true`.

- `mito.mito_db` and `mito.biof_db` (if `mito.enabled`)
  - Must exist and be usable by BLAST.

### Optional references (auto-download fallback exists)
- `sina.reference` or `sina.reference_url`
- `taxonomy.ref_taxonomy` or `taxonomy.ref_taxonomy_url`
- `taxonomy.ref_sequences` or `taxonomy.ref_sequences_url`

If local files are not present, the pipeline attempts downloads to output subdirectories.

---

## Config Sections You Should Edit First

### `paths`
- `input_dir`
- `output_dir`
- `manifest` (recommended)
- `work_dir` (optional override; otherwise Nextflow uses `nextflow.config` default)
- `conda_cache_dir` (optional override)

### `resources`
- `threads`
- `single_end`

### `environments`
Map logical stage environments to YAML files under `ASPIRE/processes/`.

Common keys:
- `main`, `sina`, `taxonomy`, `mitomaster`, `mito_checker`, `filter_counts`, `general_stats`, `plot_metadata`, `sankey`
- Optional advanced keys: `diversity`, `indicspecies`, `power_analysis`, `clustermaps`, `spieceasi`, `asv_mag_link`, `network`, `network_modules`, `plot_upset`, `bubbleplotter`, `umap_clustering`

### Feature toggles to review
- `mito.enabled`
- `filter_counts.enabled`
- `metadata_plots.enabled`
- `sankey.enabled`
- `batch_correction.enabled`
- `outlier_detection.enabled`
- `collectors_curve.enabled`
- `diversity.enabled`
- `indicspecies.enabled`
- `power_analysis.enabled`
- `clustermaps.enabled`
- `spieceasi.enabled`
- `asv_mag_link.enabled`

### `power_analysis`
Optional patient-aware power-analysis branch. This stage is off by default and requires `metadata_plots.enabled: true`, because it builds a long-format master table from `ASV_meta_micro*.tsv` plus the final micro count table. If `indicspecies.enabled: true`, the power-analysis plotting wrapper will also reuse the run's indicspecies outputs for aligned ISA power figures.

### `diversity.patient_aware`
Optional patient-aware Bray-Curtis sub-analysis nested under `diversity`. When enabled, `DIVERSITY_ANALYSIS` still produces the standard Shannon/Bray/Jaccard outputs, then also runs the patient-aware Bray PERMANOVA workflow into `diversity/<patient_aware.output_dir>/`. This covers the old standalone Bray patient-aware branch without requiring a separate workflow module.

Common keys:
- `enabled`
- `output_dir`
- `sample_col`
- `patient_col`
- `case_col`
- `type_col`
- `sample_types`
- `exclude_contralateral_in_cancer`
- `contralateral_col`
- `cancer_site_col`
- `lung_side_col`
- `contralateral_value`
- `contralateral_sample_types`
- `transform` (`none` or `rclr`)
- `permutations`
- `seed`
- `require_complete_types`

### Flexible ISA Groups
`indicspecies.group_cols` can contain more than two grouping columns. The ISA analysis itself already runs all listed groups, and the downstream network layers support `groupN_*` ISA modes based on that configured order.

Useful keys:
- `indicspecies.group_cols`
- `indicspecies.group_palettes`
- `indicspecies.group_orders`
- `indicspecies.focus_labels`
- `spieceasi.isa_overlay_groups`
- `spieceasi.group_palettes`
- `spieceasi.group_orders`
- `spieceasi.focus_labels`

Example:

```yaml
indicspecies:
  group_cols:
    - anomaly_type
    - o2_subcompartment_final
    - o2_compartment
  group_orders:
    anomaly_type: [oxic, dysoxic, anoxic]

spieceasi:
  isa_overlay_groups:
    - anomaly_type
    - o2_subcompartment_final
    - o2_compartment
  network_modes:
    - group1_isa_all
    - group2_isa_all
    - group3_isa_all
```

Common keys:
- `enabled`
- `output_dir`
- `sample_col`
- `patient_col`
- `case_col`
- `type_col`
- `sample_sizes_cancer`
- `sample_sizes_stype`
- `n_control`
- `n_simulations`
- `n_perm`
- `alpha`
- `seed`
- `skip_estimate`
- `skip_plot`
- `transform` (`none` or `rclr`)
- `keep_contralateral_in_cancer`
- `contralateral_sample_types`

### `taxonomy_patient_aware`
Optional patient-aware taxonomy branch. This stage is off by default and requires `metadata_plots.enabled: true`. It builds a fresh `ASV_master_long.tsv` from the current run outputs, runs `run_taxonomic_abundance_analysis.py` for case/control comparisons and `run_taxonomic_sample_type_analysis.py` for paired sample-type comparisons, then renders the combined figures with `plot_taxonomic_observed_analysis.py`. When batch correction is enabled, the branch uses the corrected count matrix and corrected ASV metadata so the taxonomy summaries reflect the same processed data used by the rest of the run.

Common keys:
- `enabled`
- `output_dir`
- `sample_col`
- `patient_col`
- `case_col`
- `type_col`
- `count_col`
- `tax_levels`
- `sample_types`
- `min_prevalence`
- `exclude_contralateral_in_cancer`
- `contralateral_col`
- `cancer_site_col`
- `lung_side_col`
- `contralateral_value`
- `contralateral_sample_types`
- `skip_omnibus`
- `transform` (`none` or `rclr`)
- `alpha`
- `top_n`

### `lung_status_analysis`
Optional lung-status contrast branch. This stage is off by default and requires `metadata_plots.enabled: true`. It builds a fresh `ASV_master_long.tsv`, then for each configured sample type runs `prepare_lung_status_data.py`, `run_lung_status_analysis.R`, and `plot_lung_status_analysis.py`. The prep step prefers explicit metadata columns such as `TumorSide`, `Contralateral`, `Healthy`, and `lung_status` when present, and otherwise falls back to deriving the three-way status from `Case`, `Cancer_Site`, and `lung_code`. When batch correction is enabled, the branch uses the corrected count matrix and corrected ASV metadata.

Common keys:
- `enabled`
- `output_dir`
- `sample_col`
- `type_col`
- `sample_types`
- `case_col`
- `patient_col`
- `cancer_site_col`
- `lung_code_col`
- `tumor_side_col`
- `contralateral_col`
- `healthy_col`
- `lung_status_col`

### `asv_mag_link`
Optional ASV-to-genome linkage branch. This stage is off by default and runs after the ASV feature table is finalized but before `MASTER_SUMMARY`, so its outputs can be pulled into the summary tables. It aligns the filtered ASV FASTA against barrnap-derived 16S/SSU references recovered from one or more genome/MAG sources, then writes per-hit tables, best-hit pairings, pairing-status summaries, genome-level count tables, MAG-enriched summary tables, and downstream plots. The preferred input is `genome_qc_dir` or `genome_qc_dirs`, because the linker can autodetect `barrnap/`, prefer final QC genomes from `dedupe/fasta` (falling back to `genomes_subset`), and restrict the eligible MAGs to those that pass the genome-QC barrnap check. If you are not using a genome_qc result directory, you can still provide `barrnap_dir` directly, and optionally `genome_fasta_dir` for GFF-based interval extraction.

Common keys:
- `enabled`
- `genome_qc_dir`
- `genome_qc_dirs`
- `barrnap_dir`
- `genome_fasta_dir`
- `output_dir`
- `threads`
- `min_pident`
- `min_qcov`
- `top_n`
- `plot_top_n`

### `indicspecies.aligned_*`
Optional aligned ISA summary/plot stage built on top of the existing `INDICSPECIES` outputs. This stage is off by default and does not replace the standard `INDICSPECIES_PLOTS` step; it adds a second summary/visualization pass using `plot_indicspecies_aligned.py`.

Common keys:
- `aligned_plot_enabled`
- `aligned_plot_output_dir`
- `aligned_alpha`
- `aligned_min_stat`
- `aligned_top_n`

---

## Manifest vs Auto-Discovery

If `paths.manifest` is set, the pipeline uses it.

If not set, the pipeline scans `paths.input_dir` and detects pairs using:
- `filename_patterns.r1_tokens`
- `filename_patterns.r2_tokens`
- `filename_patterns.ext_patterns`
- `filename_patterns.sample_strip_regex`

This works for common Illumina names, but manifest mode is safer and reproducible.

---

## Output Structure (Top-Level)

Under `paths.output_dir`, typical directories include:
- `fastp/`, `merged/`, `filtered/`, `concat/`, `derep/`, `denoise/`, `nochimeras/`, `ASVs/`
- `sina/`, `taxonomy/`, `mito/`, `stats/`, `logs/`
- optional: `metadata/`, `batch_correction/`, `outliers_corrected/`, `diversity/`, `indicspecies/`, `clustermaps/`, `spieceasi/`, `network/`
- optional: `power_analysis/`
- optional: `taxonomy_patient_aware/`
- optional: `lung_status_analysis/`
- optional: `asv_mag_link/`

Nextflow runtime artifacts:
- Work dir: `${HOME}/.nextflow/work` unless overridden
- Local run log: `.nextflow.log`

---

## Practical Troubleshooting

- `No usable entries detected in manifest ...`
  - Check tab format, sample IDs, and file paths in manifest.

- `metadata_plots metadata file not found ...`
  - Provide `metadata_plots.metadata` or disable `metadata_plots.enabled`.

- `power_analysis.enabled requires metadata_plots.enabled to be true`
  - Enable `metadata_plots`, or disable `power_analysis`.

- `... BLAST database not found ...`
  - Set `mito.mito_db` / `mito.biof_db` to valid DB prefixes and ensure files exist.

- `Sankey requires filter_counts.save_intermediates to be true`
  - Set `filter_counts.save_intermediates: true` when `sankey.enabled: true`.

- Conda solve issues
  - Confirm `mamba` is installed in your runner environment.
  - Pin/override env YAMLs via `environments.*` keys.

---

## Example Run Commands

```bash
# Standard run using custom YAML
nextflow run asv_pipeline.nf --config my_run.yml

# Resume
nextflow run asv_pipeline.nf --config my_run.yml -resume

# Run with explicit Nextflow reports
nextflow run asv_pipeline.nf --config my_run.yml \
  -with-report reports/nf_report.html \
  -with-trace reports/nf_trace.tsv \
  -with-timeline reports/nf_timeline.html
```
