# ASPIRE Process and I/O Reference

This catalogue describes the processes available on `main`. Exact filenames
are inventoried after each run in
`<output_dir>/summary/tables/module_output_manifest.tsv`; the paths below show
the stable public destination or the principal intermediate family.

Every process has a tier-qualified user-facing identifier. Use these identifiers
with `./run_asv_pipeline.sh --rerun-from`; unqualified names remain accepted for
compatibility:

- `core:PROCESS`: required ASV construction and taxonomy stages.
- `standard:PROCESS`: final-table preparation used by the usual analysis graph.
- `optional:PROCESS`: selectable analyses and reporting stages.

The launcher prints the same three groups before execution, so a reviewer can
distinguish the pipeline backbone from optional branches in the terminal.

## Core ASV Workflow

| Process | Principal input | Principal output | Purpose |
|---|---|---|---|
| `FASTP_QC` | manifest FASTQ pairs | `intermediates/fastp/` | Fixed trimming and read-level QC. |
| `MERGE_READS` | trimmed paired reads | merged FASTQ/FASTA intermediates | Merge overlapping read pairs. |
| `FILTER_READS` | merged reads | quality-filtered FASTA | Expected-error and length filtering. |
| `RELABEL_FILTERED` | per-sample filtered FASTA | sample-labelled FASTA | Preserve sample identity in sequence headers. |
| `CONCAT_FASTAS` | labelled sample FASTAs | concatenated read FASTA | Build study-wide denoising and mapping inputs. |
| `DEREPLICATE` | concatenated reads | dereplicated sequences | Collapse identical sequences and retain abundance. |
| `DENOISE` | dereplicated sequences | denoised ASV candidates | Remove inferred sequencing errors. |
| `CHIMERA_CHECK` | ASV candidates | non-chimeric ASV FASTA | Remove chimeric sequences. |
| `CREATE_COUNT_MATRIX` | labelled reads and ASV FASTA | raw ASV-by-sample matrix | Map reads back to ASVs. |
| `FILTER_TABLE` | raw matrix and ASV FASTA | technically filtered counts/FASTA | Enforce sample-depth and ASV-total thresholds. |
| `SINA_TRIM` | filtered ASV FASTA | aligned/region-trimmed ASVs | Identify and trim configured 16S variable regions. |
| `TAXONOMY` | trimmed ASVs and SILVA artifacts | `modules/taxonomy/tables/` | Assign taxonomy and classification statistics. |

## Non-Target Filtering and Final Tables

| Process | Principal input | Principal output | Purpose |
|---|---|---|---|
| `PREPARE_BLAST_DATABASES` | configured BLAST databases or FASTAs | `references/reference/blast_databases/` | Archive/rebuild reproducible run databases. |
| `MITOMASTER` | filtered ASVs | MITOMASTER result/chunk tables | Query candidate mitochondrial sequences when enabled. |
| `MITO_DECONTAM` | taxonomy, MITOMASTER and local BLAST evidence | `modules/non_target_filtering/` | Combine non-target evidence and audit calls. |
| `FILTER_COUNTS` | filtered matrix, taxonomy, non-target calls, metadata | `ASV_target.tsv` plus audit tables | Produce final host/non-target/abundance-filtered microbial counts. |
| `PLOT_METADATA` | final microbial and mitochondrial counts, metadata | `modules/metadata_plots/` | Build `ASV_meta_micro.tsv`, `ASV_final.micro.tsv`, updated metadata and master tables. |
| `THREE_TIER_DECONTAM` | raw control-bearing counts, final long/wide tables, control metadata | `modules/contamination_filtering/` | Score control prevalence/frequency/plausibility and filter downstream ASV tables. |

`ASV_target.tsv` is the final microbial table emitted by `FILTER_COUNTS`.
`ASV_final.micro.tsv` is its downstream wide sample-by-ASV representation after
metadata/control handling. When three-tier or batch correction is enabled, the
corresponding filtered/corrected wide table becomes the downstream count input.

## Descriptive and Metadata Processes

| Process | Input | Public output |
|---|---|---|
| `GENERAL_STATS` | concatenated reads | `modules/general_stats/` |
| `SANKEY` | stage count summaries and metadata | `modules/sankey/` |
| `PLOT_UPSET` | metadata-linked ASV tables | `modules/upset/` |
| `BUBBLEPLOTTER` | long ASV/metadata table | `modules/bubbleplotter/` |
| `UMAP_CLUSTERING` | long ASV/metadata table | `modules/umap_clustering/` |
| `OUTLIER_CHECKER` | ASV values and metadata | `modules/outlier_detection/` |
| `COLLECTORS_CURVE` | wide ASV table and metadata | `modules/collectors_curve/` |
| `CLUSTERMAPS` | long/wide ASV tables, metadata and ISA completion | `modules/clustermaps/` |

## Ecological and Study-Design Processes

| Process | Principal input | Principal output |
|---|---|---|
| `ASV_BATCH_CORRECTION` | metadata-linked ASV tables | `modules/batch_correction/` corrected counts and diagnostics |
| `ASV_META_FROM_CORRECTED` | corrected wide counts and original long annotations | regenerated corrected long table for downstream modules |
| `DIVERSITY_ANALYSIS` | downstream counts and metadata | `modules/diversity/` Shannon, Bray-Curtis, Jaccard, PERMANOVA/PERMDISP and figures |
| `INDICSPECIES` | downstream counts and groups | `modules/indicator_analysis/tables/` complete and summary ISA tables |
| `INDICSPECIES_PLOTS` | ISA tables and metadata | `modules/indicator_analysis/plots/` |
| `INDICSPECIES_ALIGNED_PLOTS` | collected ISA tables | aligned ISA tables and figures |
| `POWER_ANALYSIS_PIPELINE` | downstream ASV summaries and ISA completion | `modules/power_analysis/` simulation results and plots |
| `TAXONOMY_PATIENT_AWARE` | downstream ASV summaries | `modules/taxonomy/` patient-aware abundance tests and plots |
| `LUNG_STATUS_ANALYSIS` | downstream ASV summaries and lung metadata | `modules/lung_status_analysis/` prepared cohorts, tests and plots |

## Generalized Optional Processes and Aliases

The generalized stage names below retain the older names as launcher aliases.
Both flat legacy YAML and tiered YAML remain accepted; define a section only once.

| Process | Inputs | Outputs / configuration |
|---|---|---|
| `GROUP_POWER_ANALYSIS` | downstream ASV tables and ISA completion | Same power analysis as `POWER_ANALYSIS_PIPELINE`; `optional.power_analysis` remains supported. |
| `TAXONOMY_GROUP_ASSOCIATION` | downstream long/wide tables and grouping metadata | Taxonomic abundance tests; alias `TAXONOMY_PATIENT_AWARE`, configured with `optional.taxonomy_patient_aware`. |
| `PAIRED_GROUP_CONTRAST` | downstream tables and paired-group metadata | Paired comparisons; alias `LUNG_STATUS_ANALYSIS`, configured with `optional.lung_status_analysis`. |
| `GROUPING_DIAGNOSTICS` | downstream counts and sample metadata | Group separation, ordinations, optional balanced-resampling support, and optional proposed labels in `modules/grouping_diagnostics/`. |
| `GROUP_LABEL_AUGMENTATION` | diagnostic assignments, validation summary, metadata and downstream long table | Audited augmented metadata/long table; only changes downstream labels when `optional.grouping_diagnostics.soft_labeling.apply_downstream` is enabled and validation gates pass. |
| `MEASUREMENT_ASSOCIATION` | downstream counts/long table, metadata and configured measurements | Spearman associations and configured CCA/RDA/dbRDA analyses in `modules/measurement_association/`; complementary to the dedicated VOC process. |
| `ASV_MAG_NETWORK` | inferred network, taxonomy, ASV-MAG links and optional MAG abundance/functions | Taxonomy-filtered network mappings and functional/abundance summaries in `modules/asv_mag_network/`. Requires ASV-MAG linking and network analysis. |

These additional modules remain optional. Their full settings appear in the
[parameter catalogue](CONFIG_PARAMETERS.md). Grouping diagnostics and label
augmentation use the three-tier-filtered counts/annotations when that filter is enabled.

## VOC Process Details

`VOC_CORRELATION` receives the downstream long and wide ASV tables, metadata,
the configured VOC table, and collected ISA outputs. It publishes under
`modules/voc_correlation/`:

- matched sample and participant audit tables;
- baseline ASV-VOC Spearman coefficients, p-values and FDR values;
- direction- and magnitude-filtered ASV-VOC matrices and clustermaps;
- ISA-focused all-sample-type and bronchial/lung subset tables and figures;
- sample VOC z-score heatmaps and participant/case VOC displays.
- patient-level relative-abundance permutation tests and CLR sensitivity tables;
- patient-level focused clustermaps, testability reasons and a settings summary;
- tie-aware permutation cancer–control VOC tests on patient averages.

Baseline sample-level p-values are exploratory because they ignore repeated
patients. Use the new patient-level inference tables for statistical reporting.
See [VOC statistics](VOC_STATISTICS.md) for assumptions, FDR families and outputs.

The sample VOC heatmaps display standardized VOC abundance, whereas ASV-VOC
heatmaps display Spearman correlations. They are complementary outputs and
must not be interpreted as the same statistic.

## Network and Genome Processes

| Process | Principal input | Principal output |
|---|---|---|
| `SPIECEASI` | downstream wide ASV table and optional ISA retention table | `modules/network_analysis/` matrices, GraphML, edge/node tables and cohort audit |
| `NETWORK_TOPOLOGY` | thresholded SPIEC-EASI GraphML | observed topology summary and degree-preserving null draws |
| `NETWORK_MODULES` | full and thresholded graphs | consensus module assignments and stability summaries |
| `GRAPH_NETWORK` | graphs, ASV counts, metadata, taxonomy, ISA and modules | editable network figures and annotated node tables |
| `ASV_MAG_LINK` | filtered ASV FASTA and genome/barrnap inputs | `modules/asv_mag_link/` match tables and plots |
| `MODULE_MAG_ANCHORS` | module, node, taxonomy, abundance and ASV-MAG products | module anchor, sample score and summary tables |

## Integrated Publication

`MASTER_SUMMARY` joins the final ASV annotations with configured indicator,
network, VOC, clustermap and optional MAG outputs, then writes integrated tables
under `summary/tables/`. After Nextflow succeeds,
`processes/output_layout/organize_outputs.py` atomically publishes the module
tree, builds file/checksum inventories, archives logs and creates
`summary/report/ASPIRE_run_report.html`.

The report is the first review page; the TSV inventories and individual module
tables remain the authoritative machine-readable results.
