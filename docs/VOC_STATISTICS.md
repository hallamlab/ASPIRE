# VOC analysis: exploratory displays and patient-level inference

The optional `VOC_CORRELATION` process preserves existing matching, QC, taxonomy,
ISA membership and sample-level displays. Patient inference is enabled by default
when VOC analysis is enabled. It does not change upstream ASVs or any non-VOC
analysis. Output files are published under `modules/voc_correlation/`.

## Interpretation

- Existing `asv_voc_*` and `isa_*asv_voc_*` outputs use sample-level raw-count
  Spearman correlations. Their historical p/q values ignore repeated patients;
  treat these outputs as exploratory, not evidence of independent associations.
- New `patient_*asv_voc_*` outputs test **between-patient** relationships. They
  do not estimate within-patient effects or adjust for cancer, batch, environment,
  detection limits or other confounders. ISA-based selection uses the same study
  and does not turn these analyses into independent validation.
- Neither normalization estimates absolute bacterial abundance. An effect-size
  display cutoff such as `isa_min_abs_rho: 0.35` is not a significance threshold.
  ISA significance and ASV–VOC significance are separate questions.

## Patient inference

1. Use the already matched samples and the full supplied QC ASV count matrix.
   Exclude zero-depth samples from ASV inference. Missing patient IDs or conflicting
   patient case labels raise an error; unknown case labels are not controls.
2. **Primary:** divide each sample's ASV counts by its full retained ASV total.
   **Sensitivity:** add `clr_pseudocount` to all counts, log-transform, and subtract
   the sample mean log-count. Normalize before selecting candidate ASVs.
   CLR is sensitive to zeros, the pseudocount and the retained ASV universe; it
   is a sensitivity check, not a guarantee against compositional bias.
3. For each VOC, average normalized ASVs and VOC values within each patient using
   only cognate samples with a finite measurement for that VOC. This gives equal
   weight to patients, not reads or numbers of repeated samples. These are means
   on the input VOC scale; ASPIRE does not undo upstream log transformations.
4. Test the union of globally eligible ASVs and significant Type_Group ISA ASVs.
   Defaults require at least six patients, three with detected ASV counts, and
   variation in both variables. These gates do not guarantee adequate power.
   Untestable pairs remain in the long table with a reason and missing p/q values.
5. Compute two-sided Spearman statistics. Shuffle VOC patient ranks (not sample
   rows) using a fixed seed and `patient_permutations` draws. Use
   `(extreme + 1)/(draws + 1)` with an absolute-statistic comparison. This tests
   exchangeable patient-level independence, not causality or covariate-adjusted
   association. Missing VOC values are not replaced by zero; censored values
   require appropriate preparation and are not modeled as censoring here.
   Patient averages are rounded to 14 decimal places before ranking to stabilize
   floating-point ties; identical repeated values must not gain distinct ranks.
6. Apply BH correction over all tested union ASV–VOC pairs separately for primary
   and CLR analyses, before any plot selection. ISA/global/bronchial subsets
   inherit these q-values; they do not receive separate, smaller FDR families.
   CLR is a sensitivity analysis, not an additional route to claiming discovery.

Both signs are retained in patient tables/plots regardless of legacy display
direction. Focused plots use `isa_min_abs_rho` on rows AND columns. Full tested
and untested pairs remain in the long table. When a matrix has missing cells,
the table is saved but its clustered plot is skipped rather than inventing zero
correlations. Patient ASV matrix exports average all matched positive-depth
samples; with missing VOCs the pair-specific averages used in tests may differ.

## Cancer–control VOC comparison

Average repeated VOC measurements within each patient in **both** groups. Use
a two-sided Mann–Whitney rank-sum permutation test, retaining ties. Enumerate
all group allocations if their number is no greater than `patient_permutations`;
otherwise use seeded Monte Carlo allocations with the +1 correction. Require
at least three observed patients per group. BH correction covers all tested VOCs.
This compares distributions, not specifically medians unless distributional
shape assumptions hold. Median differences are descriptive effect summaries.
Bars show mean ± SD of patient-level z-scores; tests use unstandardized patient
means on the supplied VOC scale. Unknown cases are excluded from case comparisons.

## Configuration and outputs

All options are under `optional.voc_correlation` in the complete template:
`patient_inference`, `patient_permutations`, `patient_seed`,
`patient_min_patients`, `patient_min_nonzero`, `clr_pseudocount`.
The mock template uses 999 permutations for demonstration; the full template and
SPARK-compatible local configuration use 9999. These settings never change SPARK
upstream analyses. `patient_inference: false` disables additional ASV inference,
but cancer–control tests still use the corrected permutation method.

New files:

- `patient_asv_voc_permutation_long.tsv`: ASV, VOC, normalization, patient counts,
  rho, p, q, test status, permutations, family membership, taxonomy labels.
- `patient_asv_voc_normalization_comparison.tsv`: both normalizations side by
  side; agreement, insufficient observations and normalization-sensitive support.
- `patient_asv_relative_abundance.tsv`, `patient_asv_clr.tsv`: patient averages
  across the full ASV universe, before pair-specific missing-VOC exclusions.
- `patient_{all_asv,isa_all_sample_types,isa_bronchial_brush}_asv_voc_spearman.tsv`
  and corresponding `*_clustermap.{pdf,png,svg}`: focused primary displays.
- `patient_inference_summary.json`: denominators, gates, settings and test counts.
- `patient_voc_case_tests_brush.tsv`: updated case tests, including method and
  number of permutations. With no testable VOCs, no case-test table is produced;
  any case barplot is descriptive and has no inferential annotation.

See [the complete parameter catalogue](CONFIG_PARAMETERS.md) for defaults.
