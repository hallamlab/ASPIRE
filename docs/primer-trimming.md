# Primer removal before fastp

The optional `core:PRIMER_TRIM` module uses Cutadapt 5.2 to remove paired primers from raw reads. `core:PRIMER_TRIM_CHECK` then verifies that every sample belongs to the same amplicon family before releasing reads to fastp, merging and ASV construction. Biological samples and all control classes use this same preprocessing.

## Configuration

For a known 515F/926R assay, add this to your existing `core` mapping:

```yaml
core:
  primer_trimming:
    enabled: true
    primers:
      - name: 515F_926R
        forward: GTGYCAGCMGCCGCGGTAA
        reverse: CCGYCAATTYMTTTRAGTTT
    sample_reads: 5000
    max_prefix: 12
    min_pair_fraction: 0.5
    min_family_fraction: 0.01
    min_family_reads: 5
    screen_error_rate: 0.0
    error_rate: 0.1
    end_overlap: 12
    minimum_length: 1
  fastp:
    trim_front_r1: 0
    trim_tail_r1: 0
    trim_front_r2: 0
    trim_tail_r2: 0
```

The template enables primer trimming, sets all four fixed fastp clipping values to zero, and includes both 515F/806R and 515F/926R candidates. For a known V4 assay, use only `515F_806R` with forward `GTGYCAGCMGCCGCGGTAA` and reverse `GGACTACNVGGGTWTCTAAT`. Custom families require a unique name and forward/reverse IUPAC DNA sequences. Candidate screening identifies support among the supplied sequences; it cannot discover unknown primers or establish the exact oligo formulation used in the laboratory.

All four fixed fastp clipping values must explicitly be zero when this module is enabled. Fastp still performs quality filtering and its usual adapter processing. Configure merged-read length limits and SINA regions for your assay: V4 and V4–V5 inserts need different length settings. Already primer-trimmed reads should use this module disabled and appropriate fastp clipping values.

## Detection and filtering

The first 5,000 pairs per sample are screened for exact IUPAC-compatible paired primer matches, in either mate orientation and allowing up to 12 leading bases. At least 50% of screened pairs must support the catalogue. Each selected family/orientation requires at least five supporting pairs and 1% of screened pairs. These are screening thresholds, distinct from Cutadapt's mismatch rate. `screen_error_rate` defaults to zero; for noisy primer reads it can allow substitutions (for example, 0.1 permits one mismatch in a 19-base primer). It requires a complete primer, permits no indels and does not change Cutadapt's separate `error_rate`. Inspect primer evidence before changing either rate.

Cutadapt removes the full anchored 5′ primer plus detected prefix, permits a 0.1 mismatch rate without indels, and removes optional opposite-primer read-through at the 3′ end with a minimum overlap of 12 bases. Both mates must match a selected primer pair. Unmatched pairs and pairs with a mate shorter than `minimum_length` are retained in separate audit bins and excluded from fastp input. Either detected mate orientation is supported; original mate ordering is retained.

Unsupported samples and samples containing multiple selected amplicon families stop the workflow. Different families across samples also stop the cohort before fastp. Inspect the failed task's `audit/` reports and confirm the assay, primer catalogue and screening thresholds. Split genuinely different amplicon families into separate runs. This module requires paired-end inputs; low-depth controls must also supply enough supported pairs for the configured screen.

## Outputs and reruns

After successful publication:

- `modules/primer_trimming/tables/<sample>/`: primer detection, input/trimmed manifests, summary TSV, run settings, Cutadapt JSON and command logs.
- `modules/primer_trimming/tables/<sample>.primer_summary.json`: sample counts and portable paths to retained reads.
- `modules/primer_trimming/tables/primer_cohort.json`: cohort family and all sample summaries.
- `intermediates/primer_trimmed/<sample>/`: `trimmed/`, `unmatched/` and `too_short/` paired FASTQs.

The detailed utility manifests record task-local paths; the JSON summaries provide public retained-read paths. Original FASTQs remain unchanged. Primer audit counts explain losses before fastp; fastp's own input counts describe primer-retained reads.

Run with `./run_asv_pipeline.sh study.yml`. To rerun this stage and its dependents, use `./run_asv_pipeline.sh --rerun-from core:PRIMER_TRIM study.yml`. The module has its own managed environment in `processes/primer_trimming/env.yml`; an existing environment can be selected with `environments.primer_trimming`.

The underlying utility also supports [standalone use](https://github.com/hallamlab/ASPIRE/blob/main/scripts/PRIMER_TRIMMING.md). The integrated workflow adds manifest IDs, strict family checks, publication and dependency handling.
