# SPARK three-tier decontamination

This directory supports two related modes. The integrated Nextflow
`THREE_TIER_DECONTAM` process is optional, YAML-configurable, and runs between
`PLOT_METADATA` and every downstream analysis. It learns statistical flags from
the raw control-bearing count matrix, applies them to the final microbial long
and wide tables, and publishes its full audit under the configured
`three_tier_decontam.output_dir`.

The standalone wrapper below is the stricter manuscript-resumption process. It
starts from an existing authoritative ASPIRE output and never modifies that
directory. It reconstructs the historical
1,425-ASV long and 1,377-ASV corrected-wide inputs, executes the exact decontam
scripts deposited in `SPARK_analysis_code.zip`, and requires exact SD3 table
equality plus SD1 flag/score agreement before downstream analysis can begin.

The process belongs after ASPIRE's corrected final-table construction and before
diversity, indicator-species, taxonomy, contralateral, or network analysis.

The 1,388-ASV filtered long table is retained as the tier audit artifact. A
separate downstream long projection is required to match validated SD3 exactly
(1,342 ASVs, 125 samples, 2,872,255 reads); all analyses consume that projection.

For analyses that intentionally improve on SPARK by excluding residual host
reads, run `filter_host_taxa.py` after SD3 validation and route every downstream
branch from its paired host-filtered long and wide tables. This remains a
separate checkpoint so the exact published SPARK input is preserved.

The extracted analysis copy receives one compatibility-only edit for current
pandas/scikit-bio: PCoA variance uses positional `.iloc` indexing. Statistical
methods and values are unchanged.

Run the wrapper with a new, nonexistent output directory. It refuses to overwrite
an earlier run.

```bash
processes/three_tier_decontam/run_spark_manuscript_pipeline.sh \
  --final-dir /path/to/spark_set1-2_output_final \
  --metadata /path/to/spark_set1-2_metadata.tsv \
  --submission-package /path/to/SPARK_submission_package_final \
  --output-dir /path/to/new_SPARK_manuscript_output \
  --rscript /path/to/decontam/environment/bin/Rscript \
  --downstream-rscript /path/to/manuscript-R/environment/bin/Rscript \
  --downstream-python /path/to/manuscript-python/environment/bin/python
```

The deposited indicator-species script does not set an RNG seed, so its ASV
memberships and indicator statistics reproduce but permutation p-values need not
be byte-identical. The published SD9 consensus Leiden assignments are not
recomputable from the submission package: its README identifies them as an
upstream-facility result and the exact network invocation was not recorded.
