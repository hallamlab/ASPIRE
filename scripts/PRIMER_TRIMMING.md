# Standalone paired-read primer trimming

`trim_amplicon_primers.py` is an optional preprocessing utility, not an ASPIRE
workflow step. It leaves input FASTQs unchanged and requires a new or empty output
directory. Python standard library and Cutadapt >=4 are required; tested with 5.2.

```bash
mamba create -n aspire_primer_qc -c conda-forge -c bioconda python=3.11 cutadapt=5.2
mamba activate aspire_primer_qc
python ~/repos/ASPIRE/scripts/trim_amplicon_primers.py \
  --fastq-dir ~/data/ASV/SI_fastq \
  --output-dir ~/data/ASV/SI_fastq_primer_trimmed \
  --threads 8
```

Use the actual raw FASTQ directory if different. Filenames must end in `_R1`/`_R2`,
`_R1_001`/`_R2_001`, `_1`/`_2`, or equivalent dot-separated suffixes, followed by
`.fastq` or `.fq`, optionally `.gz`. Lane/chunk prefixes are preserved as separate
sample IDs; files are not concatenated. Use `--recursive` for nested directories.
Duplicate sample IDs, missing mates and unrecognized FASTQ filenames are errors.

## Detection and trimming

The first 5,000 pairs per sample are screened for exact IUPAC-compatible matches
to 515F/806R and 515F/926R in either mate orientation, allowing 0–12 prefix bases.
These are broad primer families: the 515F-Y/806RB degenerate sequences also match
the original formulations. A match cannot conclusively identify the oligo version
used in the laboratory. This is candidate screening, not de novo primer discovery.
Already-trimmed or unknown-primer data may have no supported match.

At least 50% of screened pairs must support a candidate pair. Each selected
family/orientation needs at least five pairs and 1% of screened pairs. Adjust with
`--sample-reads`, `--max-prefix`, `--min-pair-fraction`, `--min-family-reads`, and
`--min-family-fraction`. Screening the first N records can miss rare or later
families; it is not a random sample. Both forward and swapped orientations are
reported; output mates remain in their original order.

Cutadapt processes the full sample with selected candidates, requiring full 5′
primer matches at allowed prefix lengths. Prefix and primer are removed together.
Trimming allows a 0.1 mismatch rate (`--error-rate`) and no indels. It also removes
optional 3′ read-through of the opposite primer's reverse complement, requiring
12 bases of overlap (`--end-overlap`). Both mates must have a matched 5′ primer
and meet `--minimum-length` (default 1) to enter the main output. Matching on each
mate is independent: this is not per-read primer-family demultiplexing.

Mixed families are flagged and trimmed with all supported candidates, but remain
in the same sample FASTQs. **Do not pool different amplified regions for ASV
inference.** Review or split mixed-region libraries before downstream analysis.
Unsupported samples produce reports but no trimmed FASTQs; processing continues
for other samples. Exit status 2 indicates unsupported or failed samples; 0 means
all samples were processed (including flagged mixed-family samples).

To provide another candidate catalogue, use `--primers-json primers.json`, where
both primers are written 5′→3′ as synthesized (not reverse-complemented):

```json
[
  {"name": "515F_806R", "forward": "GTGYCAGCMGCCGCGGTAA", "reverse": "GGACTACNVGGGTWTCTAAT"},
  {"name": "515F_926R", "forward": "GTGYCAGCMGCCGCGGTAA", "reverse": "CCGYCAATTYMTTTRAGTTT"}
]
```

## Outputs and provenance

- `trimmed/`: paired primer-trimmed FASTQs.
- `unmatched/`: pairs lacking a matched primer on either mate. A matched mate may
  already be trimmed; these are not copies of the original reads.
- `too_short/`: processed pairs failing the length threshold.
- `summary.tsv`: per-sample status, selected families, paired-read counts and paths.
- `primer_detection.tsv`: screened counts, orientations and prefix offsets.
- `trimmed_manifest.tsv`: paths for successfully processed samples, including
  mixed-family samples that still require review.
- `input_manifest.tsv`: source paths, byte sizes and modification timestamps.
- `run_config.json`: parameters, primer catalogue and Cutadapt version.
- `reports/`: per-sample Cutadapt JSON, console log and exact command argument list.
- `run.log`: progress log; override with `--log-file PATH`.

Samples run sequentially, with `--threads` used by Cutadapt within each sample.
Progress prints at screening, trimming and completion; detailed Cutadapt output is
saved in each sample's report. Successful outputs are moved into place only after
Cutadapt finishes and read-pair accounting passes. Input provenance records are
not content checksums. No general quality filtering is performed here.

**Before using these FASTQs in ASPIRE, disable any existing fixed primer clipping
for these inputs.** Otherwise, a second round of clipping will remove biological
sequence. This utility does not change ASPIRE's fastp settings or launch ASPIRE.

Focused tests (integration test needs Cutadapt on PATH):

```bash
python -m unittest discover -s ~/repos/ASPIRE/tests -p test_trim_amplicon_primers.py
```

Cutadapt adapter and paired-read behavior:
<https://cutadapt.readthedocs.io/en/stable/guide.html>.
