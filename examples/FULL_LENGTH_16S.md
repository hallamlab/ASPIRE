# Full-length 16S example

`full_length_16s.local.yml` is a conservative core-processing example for
single-end, high-accuracy, near-full-length 16S amplicons. It is best matched
to PacBio CCS/HiFi consensus reads that have already been demultiplexed, primer
trimmed, and consistently oriented. Copy the manifest template, replace all
absolute paths, and review the length and expected-error settings against the
actual primer pair and sequencing run.

Run the complete core workflow with:

```bash
./run_asv_pipeline.sh examples/full_length_16s.local.yml --phase all
```

The current ASPIRE inference path uses VSEARCH UNOISE. It does not implement
PacBio's run-specific error model and should therefore be validated with a mock
community before exact variants are interpreted as ASVs. For production PacBio
CCS data, DADA2/QIIME 2 `denoise-ccs` is the usual platform-aware alternative.
The example is not appropriate unchanged for raw Oxford Nanopore reads.

Short shotgun metagenomic reads are a different input type. Reconstructing
near-full-length 16S genes from recruited shotgun reads requires a targeted
assembly program such as phyloFlash or MATAM. Those reconstructed marker
sequences are contigs or population-level reconstructions, not amplicon ASVs,
and should enter a separately labeled comparison rather than this ASV workflow.
