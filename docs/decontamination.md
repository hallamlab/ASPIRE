# Decontamination and non-target filtering

ASPIRE has two distinct filtering layers. Reference-based screens identify host,
mitochondrial and other configured non-target sequences. The optional three-tier
screen uses extraction controls, DNA concentrations and taxonomy to flag potential
contamination in the microbial tables.

## Reference-based screening

`PREPARE_BLAST_DATABASES`, `MITOMASTER` and `MITO_DECONTAM` combine local
mitochondrial/contaminant BLAST evidence, taxonomy and optional MITOMASTER
results. `FILTER_COUNTS` applies the calls together with configured abundance,
prevalence, taxonomy-quality and explicit taxon-exclusion rules.

Configure these under `standard.mito` and `standard.filter_counts`.
`standard.mito.run_mitomaster: false` disables the external lookup while retaining
local reference screens. Supply the required local databases or FASTAs for an
offline run. Review the intermediate removal tables as well as `ASV_target.tsv`.

## Three-tier control-based filtering

Enable `optional.three_tier_decontam.enabled` only after configuring the required
control-bearing metadata and raw counts. The process runs after metadata-linked
long/wide table construction and before optional group diagnostics, batch
correction and downstream statistical analyses.

| Tier | Evidence | Implementation |
|---|---|---|
| Pooled prevalence | Occurrence across extraction controls and biological samples | R `decontam` |
| Within-type frequency | Abundance versus DNA concentration within configured sample types | R `decontam` |
| Biological plausibility | Taxonomic rules for implausible airway residents | ASPIRE filtering script |

The plausibility rules were designed for airway studies; they are not a generic
list of contaminants for environmental microbiomes. Set
`optional.three_tier_decontam.biological_plausibility: false` when those rules do
not fit your study. This keeps the statistical tiers available.

The scorer uses the **raw count matrix containing controls**, rather than a table
from which controls have already been removed. The resulting flags are applied
to the final microbial long and wide tables, producing
`ASV_meta_three_tier.tsv` and `ASV_final_three_tier.tsv`. When enabled, these
replace the base tables entering downstream analysis preparation. Existing
host/non-target filtering is not undone.

Review the pooled and within-type score tables, removed-ASV flags and filtering
summary in the module's published audit. The [configuration guide](CONFIGURATION.md)
and [parameter catalogue](CONFIG_PARAMETERS.md) describe control identifiers,
concentration/sample-type columns, thresholds and combination settings.

## Workflow test and study interpretation

The [public test fixture](test.md) supplies synthetic extraction
controls and DNA concentrations to exercise the implementation. For a real
study, use the actual control design and measured concentrations. A technical
validation pass does not establish that filtering choices suit every study.

When comparing runs, record the selected final table, enabled tiers, thresholds,
removed features and read loss. Keep each parameter experiment in a separate
output directory. [Expert operations](EXPERT_GUIDE.md) explains reproducible
reruns.
