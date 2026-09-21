# ASV_MAG_NETWORK

`ASV_MAG_NETWORK` exports Kieft-style ASV association networks annotated with
sequence-derived MAG links. It intentionally reuses existing ASPIRE outputs
instead of rerunning network inference or ASV-MAG alignment.

## Inputs

- SPIEC-EASI ASV graph GraphML from `spieceasi`.
- SPIEC-EASI node feature table from `spieceasi`.
- ASV taxonomy table from `taxonomy`.
- ASV count table from metadata/batch-corrected ASPIRE outputs.
- ASV-MAG pairing table and barrnap reference catalog from `asv_mag_link`.
- Optional MAG abundance table.
- Optional functional annotation tables keyed by `genome_id`.

Ecological anchors use the shared `asv_mag_network.anchor_top_n` setting. ASVs
are ranked independently within their module by degree, eigenvector centrality,
and betweenness centrality and are retained when they qualify for at least one
metric. The default top rank is 1; tied ranks are retained.

Preferred MAG abundance input is a long-form table with an explicit value
column. For an upstream-normalized metagenome table, for example:

```tsv
genome_id	sample_id	fpkm
MAG_001	SampleA	1.234
MAG_001	SampleB	0
MAG_002	SampleA	0.055
```

Set `mag_abundance_value_col: fpkm` and
`mag_abundance_normalization: provided_fpkm` to preserve these values, or
use `tpm` and `provided_tpm` for transcript abundance. For raw paired-end
recruitment counts, provide the matching SeqKit table and use
`mag_abundance_normalization: auto` or `input_fragment_fpm`. Both settings
require the SeqKit table; legacy `median_ratio` must be selected explicitly.
R1 and R2 counts
must agree, and one mate's `num_seqs` becomes the number of input fragments.
MAG recruitment is then reported as fragments per million input fragments.
ASV counts remain sample-relative, and agreement is computed per accepted
ASV-MAG pair using only samples shared by that pair. Separate abundance and
SeqKit inputs are supported for metagenomes and metatranscriptomes.

For studies where MAG identifiers differ by source prefix/order, set
`mag_id_mode: suffix_after_double_underscore`. This keeps capitalization and
uses only the final token after the last `__`, so both
`100m__xPG_SAGs__AB-746_B09_AB-901` and
`xPG_SAGs__100m__AB-746_B09_AB-901` match as `AB-746_B09_AB-901`.

Preferred functional input is module-level genome summaries:

```tsv
bin_id	feature_id	feature_name	observed_ko_count	total_feature_kos	fraction_covered
MAG_001	M00529	Denitrification nitrate nitrogen	5	6	0.83
```

The module table is automatically subset to accepted MAGs in the analysis.

## Outputs

Outputs are written under `asv_mag_network` by default:

- `validation/`: ASV-MAG taxonomy validation and paper-threshold status.
- `mapping/`: accepted and ambiguous ASV-MAG mappings.
- `rrna/`: copied MAG 16S reference catalog and MAGs lacking accepted ASVs.
- `network/`: ASV-only paper network and heterogeneous ASV-MAG network as TSV,
  GraphML, and Cytoscape JSON.
- `abundance/`: optional ASV-MAG abundance agreement tables.
- `functional/`: optional MAG functional summaries.
- `qc/`: mapping-class summaries, plots, and run parameters.

SPIEC-EASI edges are statistical ASV-ASV associations. MAG edges are
sequence-match annotation edges and should not be interpreted as direct
biological interactions.

When ASV and MAG taxonomy sources differ, for example NCBI/SILVA ASV taxonomy
versus GTDB MAG taxonomy, conflicts are labeled
`taxonomy_unvalidated_crossdb` instead of being rejected unless a source-aware
crosswalk is added later.
