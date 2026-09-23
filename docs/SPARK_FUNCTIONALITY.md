# SPARK and ASPIRE manuscript functionality

This crosswalk identifies where a reviewer can find and run the analytical
functionality reported in the paired SPARK data paper and ASPIRE application
note. It describes capabilities, not patient-data availability. The public mock
dataset is the executable demonstration dataset.

| Reported analysis | YAML switch | ASPIRE implementation | Principal outputs |
|---|---|---|---|
| Read QC, merging, filtering, denoising, chimera removal, count matrix | core workflow | `asv_pipeline.nf`; `FASTP_QC` through `CREATE_COUNT_MATRIX` | FASTQ QC, filtered reads, ASV FASTA and count matrix |
| SILVA/SINA taxonomy and non-target screening | `sina`, `taxonomy`, `mito`, `filter_counts` | `processes/sina_trim`, `processes/taxonomy`, `processes/mitomaster`, `processes/mito_decontam`, `processes/filter_counts` | taxonomy tables, mitochondrial/non-target calls, final microbial table |
| Three-tier SPARK contamination screen | `optional.three_tier_decontam.enabled` | `processes/three_tier_decontam/pipeline` | pooled and within-type decontam scores, biological-plausibility flags, filtered long/wide tables |
| Read/sample retention and ASV overlap | `general_stats`, `sankey`, `plot_upset`, `metadata_plots` | corresponding processes listed in `README.md` | retention summaries, Sankey, UpSet and Venn outputs |
| Collector curves | `collectors_curve.enabled` | `processes/collectors_curve/collectors_curve.py` | curve figures and richness summaries |
| Shannon, Bray-Curtis, PERMANOVA and PERMDISP | `optional.diversity.enabled`; `optional.diversity.patient_aware.enabled` | `processes/diversity_analysis`; `processes/bray_patient_aware` | diversity tables/figures and patient-blocked tests |
| UMAP and hierarchical clustermaps | `umap_clustering.enabled`, `clustermaps.enabled` | `processes/umap_clustering`, `processes/clustermaps` | ordinations and taxonomic/ASV heatmaps |
| Tumour-side, contralateral and healthy-lung comparisons | `optional.lung_status_analysis.enabled` | `processes/lung_status_analysis` | prepared status tables, PERMANOVA/PERMDISP results and figures |
| Patient-aware taxonomic comparisons | `optional.taxonomy_patient_aware.enabled` | `processes/taxonomy_patient_aware` | sample-type and cancer-status taxonomic tests and figures |
| Sample-type and cancer-status indicator species | `optional.indicspecies.enabled` | `processes/indicspecies`, `processes/indicspecies_plots` | complete ISA result, summary and figure files |
| SPIEC-EASI network inference | `optional.spieceasi.enabled` | `processes/spieceasi/run_spieceasi.R` | GraphML, edge list, matrices and node centralities |
| Consensus network modules | `optional.spieceasi.modules_enabled` | `processes/network_modules/network_modules.R` | module assignments, run details and stability summaries |
| Global topology and null-network enrichment | `optional.network_topology.enabled` | `processes/network_topology/network_topology_stats.py` | `network_topology_summary.tsv`, `network_topology_null_draws.tsv` |
| Network figures and metadata/ISA overlays | `optional.spieceasi.network_enabled` | `processes/graph_network/graph_network.py` | editable SVG/PDF network panels |
| Simulation-based statistical power | `optional.power_analysis.enabled` | `processes/power_analysis_pipeline` | effect sizes, simulated power tables and figures |
| ASV-VOC correlations and participant/sample VOC displays | `optional.voc_correlation.enabled` | `processes/voc_correlation/plot_voc_corr.py` | global and ISA-focused correlations, FDR values, heatmaps and participant comparisons |

## Network topology settings

```yaml
network_topology:
  enabled: true
  output_dir: spieceasi
  n_null: 1000
  seed: 42
  skip_null: false
```

The implementation follows the deposited SPARK supplementary method: the
thresholded positive SPIEC-EASI GraphML is simplified to an undirected graph,
self-loops are removed, and configuration-model graphs preserving the observed
degree sequence are generated with a fixed seed. Observed transitivity and mean
local clustering are compared with the null distributions. Every null draw is
written so the reported fold enrichment and empirical probabilities can be
independently checked.

## Public mock-data scope

`examples/mock.local.yml` enables the manuscript analysis modules and runs the
topology method with 100 null draws. `examples/validate_mock_run.sh` verifies
that the topology tables exist, contain a non-empty graph, and contain the
completed null ensemble.

The public fixture intentionally contains no protected patient data. DECOI adds
four labelled synthetic extraction blanks, synthetic DNA concentrations, and
the corresponding control truth table so the optional three-tier prevalence and
within-type frequency models can also be exercised end to end.
