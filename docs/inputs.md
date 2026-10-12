# Inputs and references



## Inputs

FASTQs can be discovered from `core.paths.input_dir`. Every run writes a normalized,
reusable manifest to `<output_dir>/summary/tables/run_manifest.tsv`, regardless of
whether discovery or `core.paths.manifest` supplied the inputs.

Manifest format:

- Tab-separated; the `sample_id`, `fastq_r1`, `fastq_r2` header is optional.
- Column 1: `sample_id`.
- Column 2: R1 FASTQ.
- Column 3: R2 FASTQ, optional for single-end data.
- Lines starting with `#` are ignored.
- Relative FASTQ paths are resolved relative to the manifest file.

See `examples/manifest.template.tsv` for a reusable template.

Metadata is required by enabled metadata-aware branches such as metadata plots, Sankey, diversity, indicator species, VOC correlation, power analysis, taxonomy patient-aware analysis, lung-status analysis, and several network overlays. The configured sample column must match the manifest sample IDs.

Start with the {download}`minimal metadata workbook <../examples/metadata.template.xlsx>`.
It contains a **Definitions** sheet and an empty **Metadata** table with `Sample`,
`sample_class`, `Type_Group`, `Participant_ID`, `Case` and `Set`. The sample-class
column has a dropdown for biological, technical, biological-control and positive
roles. Export only the completed Metadata sheet as a UTF-8 tab-separated `.tsv`
file, then configure ASPIRE's column and group mappings as described in the workbook.

Metadata column names are configured per run (`sample_col`, `type_col`,
`case_col`, `patient_col`, and related settings); ASPIRE does not require fixed
study-specific names. If the configured color column is absent, ASPIRE assigns
deterministic colors and writes both an augmented metadata table and a reusable
two-column palette under `<output_dir>/modules/metadata_plots/tables`. Set `palette_file` in
`standard.metadata_plots` or `optional.sankey` to override it. See
`examples/metadata_palette.template.tsv`.

Reference inputs depend on enabled branches:

- `core.sina.reference` and `core.taxonomy` references are used for SINA alignment and taxonomy assignment. The config can point at local files or URLs.
- Mitochondrial/contaminant decontamination accepts existing database prefixes
  through `standard.mito.mito_db` and `standard.mito.biof_db`, or FASTA inputs through
  `standard.mito.mito_fasta` and `standard.mito.contaminant_fasta`. ASPIRE exports/copies and
  rebuilds both as run references under `<output_dir>/references/reference/blast_databases`.
- `standard.mito.run_mitomaster: false` skips the external MITOMASTER service while
  retaining taxonomy and local BLAST screens, which is useful for offline tests.
- `optional.voc_correlation.voc_table` is required when VOC correlation is enabled.
- `optional.asv_mag_link.*` inputs are required only when ASV-to-MAG linkage is enabled.
