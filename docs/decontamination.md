# Decontamination and non-target filtering

## Cohort-level control prevalence filtering

Enable `core.control_decontam.enabled` to run independent **TECH** and
**BIO** prevalence analyses with R `decontam`. The integrated algorithm uses
these two control classes.
BIO means biological/sample controls (for example scope flushes), distinct from
TECH water/extraction/library blanks. Real biological samples retain the
`biological` class; controls use `bio_control` and are exempt from the depth cutoff.
The module runs in the core workflow after ASV construction and full taxonomy
assignment, **before** reference screening, `FILTER_ASVS`, and `PLOT_METADATA`.

All supplied biological samples, biological/sample controls, technical-negative and positive-control FASTQs
share the upstream trimming, quality filtering, denoising and chimera-removal
path. ASPIRE accepts already-demultiplexed FASTQs. Configure the primers and
SINA region for the assay; for a V4 assay set `core.sina.trim_to: V4` and the
appropriate read trimming. Taxonomy now covers all inferred ASVs, including
control-only ASVs. Original `ASV_counts.tsv` is preserved.

The **≥5,000 post-QC reads cutoff is a biological-sample inclusion criterion**.
It deliberately does **not** apply to scope flushes or technical negatives.
Depth means the sum of counts mapped to non-chimeric ASVs before feature or
host filtering. Failing biological samples enter neither prevalence test nor
downstream biological analyses. Nonzero controls are retained regardless of
depth. Zero-depth controls, including metadata samples missing from the count
matrix, remain in QC reporting but do not enter the prevalence denominators.

Configure explicit metadata class values, for example:

```yaml
core:
  control_decontam:
    enabled: true
    metadata: /path/to/metadata.tsv  # null falls back to filter_counts.metadata
    metadata_sample_col: Sample
    class_col: sample_class
    biological_labels: [biological]
    technical_labels: [water_blank, extraction_blank, library_blank]
    bio_control_labels: [scope_flush]
    positive_labels: [positive]
    min_biological_reads: 5000
    technical_enabled: true
    bio_control_enabled: true
    technical_score_threshold: 0.1
    bio_control_score_threshold: 0.1
```

Values are read from `class_col`, not inferred from filenames or sample IDs.
Each value must map to exactly one class. Missing/ambiguous IDs, unknown class
values and overlapping label mappings fail validation. Hyphen/underscore ID
normalization follows the existing convention, with collisions rejected.
Keep the metadata scoped to the sequencing run; metadata-only rows are reported
as having zero usable reads. All original metadata columns, including patient,
procedure and custom pairing identifiers, are preserved in `sample_qc.tsv`.
Pairing does not replace the cohort-level statistical test.

- TECH compares QC-passing biological samples with usable technical negatives;
  biological/sample and positive controls are excluded.
- BIO compares the same biological samples with usable biological/sample controls;
  technical negatives and positive controls are excluded.
- Either arm can be disabled independently with `technical_enabled: false` or
  `bio_control_enabled: false`. An enabled arm with no usable
  controls fails explicitly; it never silently produces a clean result.
- The final contaminant set is the union of the enabled arms' `decontam` calls.
  An ASV is flagged when its decontam prevalence score is below the configured
  threshold (default 0.1), not merely because it occurs in a control. These
  scores are distinct from relative-abundance percentages and BH-adjusted q-values.
  Retained ASV counts carry forward unchanged.

The biological depth criterion is applied once. `FILTER_ASVS` and
`PLOT_METADATA` skip secondary depth cutoffs when biological QC has already run.
Relative-abundance/non-target filters then operate on the cleaned biological
counts, followed by existing analysis-specific prevalence/variance filters.
The depth cutoff is not reapplied after contaminant removal. Optional batch
correction remains a separately configured downstream operation.

## Audit and QC outputs

Published under `modules/contamination_filtering/` using the existing output
layout (tables, plots and logs are inventoried in the module manifest):

| Artifact | Content |
|---|---|
| `three_tier_results/sample_qc.tsv` | All sample metadata, class, post-QC depth and inclusion status |
| `three_tier_results/read_depth_summary.tsv` | Class-specific sample counts and depth summaries |
| `three_tier_results/read_depth_by_class.svg` | Depth distributions including low/zero-depth controls |
| `three_tier_results/biological_qc_counts.tsv` | Biological counts after sample inclusion, before ASV removal |
| `three_tier_results/TECH_counts.tsv`, `BIO_counts.tsv` | Exact input cohort for each enabled arm |
| `three_tier_results/TECH_metadata.tsv`, `BIO_metadata.tsv` | Sample order and negative flag used by R |
| `three_tier_results/TECH_scores.tsv`, `BIO_scores.tsv` | Separate decontam scores, calls and test status |
| `three_tier_results/contamination_calls.tsv` | Every original ASV, taxonomy, prevalence and denominators for all three classes, both scores/calls/statuses, final category |
| `three_tier_results/removed_asvs_with_taxonomy.tsv` | Removed ASVs and TECH, BIO or TECH+BIO reason codes |
| `three_tier_results/positive_qc_counts.tsv` | Positive-control counts kept separately for recovery/cross-talk inspection |
| `three_tier_results/step_counts.tsv` | Sample, ASV, nonzero-ASV and read counts across cohort preparation and union removal |
| `three_tier_results/settings.json` | Exact class mappings, enable flags and thresholds |
| `ASV_final_three_tier.tsv` | Cleaned biological counts entering downstream feature filters |

Final categories are CLEAN, TECH, BIO and TECH+BIO. CLEAN means neither
**enabled** arm called contamination; inspect the arm status to distinguish
`disabled`, `absent_in_arm`, `not_estimable` and `tested`. Missing scores are
never invented. Prevalence uses QC-passing biological samples and nonzero usable
controls. Positive-control counts support QC inspection; automatic expected-taxon
or index-hopping estimates require study-specific expectations and are not
calculated here.

## Migration and historical reproduction

`core.control_decontam` is the current configuration section. The previous
`optional.three_tier_decontam` (or flat `three_tier_decontam`) name remains an
alias for configurations already using the TECH/BIO settings. Do not supply
both names. The same alias applies to `environments.control_decontam`.
Historical audit filenames remain stable for existing downstream consumers.

Old `negative_control_labels`, `positive_control_labels`, frequency thresholds
and plausibility settings do not describe this algorithm. Replace them with the
class mappings and settings above. Legacy enabled configurations are rejected
with migration guidance rather than silently reinterpreted. DNA concentration
is no longer required. The CAMI mock uses `Type_Group`: Airways/Oral are biological samples, Skin is
the BIO control, and Control identifies TECH blanks. Both arms run. The older
paired-airway mock has no BIO controls and runs TECH only.

The standalone historical SPARK scripts in `processes/three_tier_decontam/`
remain available for manuscript reproduction and do not implement this new
procedure. `RUN_METADATA_ANALYSES` cannot rerun contamination removal from
already-filtered tables; use the full workflow with resume. The core-stage
rerun identifier is `core:CONTROL_DECONTAM`.

The launcher accepts old stage names as aliases: `THREE_TIER_DECONTAM` maps to
`CONTROL_DECONTAM`; `FILTER_TABLE` and `FILTER_COUNTS` map to `FILTER_ASVS`.
Historical audit filenames remain stable; use `core.control_decontam` in new YAML.

Parameter aliases preserve existing run values without changing units:

| Previous name (within its section) | Current name |
|---|---|
| `table_filter.min_sample_sum` | `table_filter.min_sample_reads` |
| `table_filter.min_asv_sum` | `table_filter.min_relative_abundance_pct` |
| `filter_counts.abundance_threshold` | `filter_counts.min_relative_abundance_pct` |
| `control_decontam.technical_threshold` | `control_decontam.technical_score_threshold` |
| `control_decontam.bio_control_threshold` | `control_decontam.bio_control_score_threshold` |
| `spieceasi.min_rel_abund` | `spieceasi.min_relative_abundance_fraction` |
| `spieceasi.min_prevalence` | `spieceasi.min_prevalence_fraction` |
| `taxonomy_patient_aware.min_prevalence` | `taxonomy_patient_aware.min_prevalence_fraction` |
| `measurement_association.min_prevalence` | `measurement_association.min_prevalence_fraction` |

The old network keys additionally accept their historical percentage shorthand
(e.g. `min_prevalence: 5`); migration converts it to fraction 0.05. New `_fraction`
keys accept only 0–1. Do not define both the previous and current key.

## Reference-based screening

`PREPARE_BLAST_DATABASES`, `MITOMASTER` and `MITO_DECONTAM` combine local
mitochondrial/contaminant BLAST evidence, taxonomy and optional MITOMASTER
results. `FILTER_ASVS` applies these separately from the TECH/BIO calls,
with configured abundance, taxonomy-quality and explicit taxon-exclusion rules.
Configure these under `standard.mito` and `standard.filter_counts`.

Set `standard.mito.run_mitomaster: false` to retain local reference screening
without an external lookup. Review the separate non-target removal tables and
`ASV_target.tsv` as well as the contamination audit. Keep parameter experiments
in separate output directories and record the selected downstream table.

## Combined feature filtering

The main path is `TAXONOMY` → `CONTROL_DECONTAM` → reference/mitochondrial
screening → `FILTER_ASVS` → `PLOT_METADATA`. `GENERAL_STATS` is a parallel QC
branch.

`FILTER_ASVS` applies the
`core.table_filter` abundance/nonzero rules first, followed by
`standard.filter_counts` group-size, non-target, microbial abundance and taxonomy
rules, then a final nonzero check. The two abundance thresholds have different
count denominators and remain independently configurable. Group sizes count
only samples still in the input count table. No depth threshold is reapplied
when control decontamination has already selected the biological cohort.

`ASV_filtered.tsv` and `ASVs_filtered.fasta.gz` retain the initial checkpoint.
The existing `.decon.tsv`, `.micro.tsv` and `.mito.tsv` outputs remain available.
`ASV_target.tsv` and `ASVs_target.fasta.gz` contain the final matching ASVs;
genome linkage now uses this final FASTA. `filter_audit.tsv` records sample/ASV/read
counts and `filter_removed_asvs.tsv` records the removal checkpoint.
`PLOT_METADATA` always receives the fully filtered microbial table.

## Relative-abundance threshold in the mock

The general and mock templates set `standard.filter_counts.min_relative_abundance_pct: 0.1`, meaning
**0.1%**, or a fraction of 0.001. This is independent of the TECH/BIO prevalence
score thresholds, which also default to 0.1 but are not percentages.

Relative abundance is calculated separately in each retained biological sample
using its microbial counts after control-ASV and reference-based removals, before
taxonomy-quality filtering. An ASV passes if it reaches **at least 0.1% in any one
biological sample**. Passing ASVs retain their original counts in every remaining
sample; individual low-abundance cells are not zeroed. Skin BIO controls, TECH
blanks and positive controls are absent from this calculation.

The initial `core.table_filter.min_relative_abundance_pct` setting is a per-sample
percentage; both templates set it to zero and apply the 0.1% gate in the final
filter. Choose the threshold in the YAML used for the run.

## Biological sample prevalence

Both templates set `standard.filter_counts.min_prevalence_fraction: 0.05`.
`FILTER_ASVS` retains ASVs with a nonzero count in at least 5% of the biological
samples remaining after group-size filtering. With 50 samples, this requires
at least three nonzero observations. Zero-depth columns created by earlier
feature removal remain in that denominator because sample inclusion was decided
upstream. Presence means any positive count, independently of the 0.1% RA gate.

Both gates are calculated from the same microbial input table after reference
removal and before taxonomy-quality filtering. `ASV_target.feature_qc.tsv`
records each ASV's biological sample denominator, nonzero sample count, prevalence,
maximum per-sample RA, configured thresholds, gate results and final retention.
The retained ASVs and original counts then feed metadata and downstream analyses.
A positive biological prevalence threshold requires control decontamination to be
enabled; the general template enables decontamination and the 5% gate together.
If you disable decontamination, set this gate to zero as required by the pipeline.
