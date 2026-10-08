# Manuscript demonstration: aspire-cami-mock

`aspire-cami-mock` is the full CAMI2-derived demonstration for the ASPIRE
manuscript. The validated run contains **50 synthetic patients and 179 libraries**.
For a small installation test, use the [quickstart](test.md).

## Dataset and study design

DECOI generated paired V4 amplicon reads from CAMI2 reference genomes and community
compositions, with seed 42 and 50,000 initial read pairs per non-blank library.
The simulation adds sequence errors, differential abundance, microbial network
signals, chemical associations, batch effects, PCR bias, contaminants,
mitochondrial reads and chimeras. Patient, disease-group and lung-side labels are
synthetic. The resulting associations demonstrate analysis behavior on known
simulation inputs.

| Cohort | Patients | Libraries |
| --- | ---: | --- |
| Control | 25 | One healthy airway, oral and skin sample per patient |
| Cancer | 25 | Tumor-side airway, contralateral airway, oral and skin per patient |
| Technical controls | — | Four extraction blanks |

This gives 75 Airways and 50 Oral biological samples, 50 Skin biological controls
and four technical controls. Each patient belongs to one of two batches. Chemical
measurements are supplied for the non-blank libraries. Original source-community
identifiers and patient pairing are preserved in metadata.

## Reproduce the demonstration

Prepare the 50-patient synthetic clinical study using the DECOI
[CAMI body-site study instructions](https://hallamlab-decoi.readthedocs.io/en/latest/cami-body-sites.html).
Use `aspire-cami-mock-inputs` for preparation and `aspire-cami-mock` for simulation
output. The completed input to ASPIRE is `dataset/mock_dataset/` beneath that
simulation directory. Its manifest has 179 libraries. ASPIRE's configurator
supplies the control classes and analysis settings described below.

From the ASPIRE repository:

```bash
DATASET=/absolute/path/to/aspire-cami-mock/dataset/mock_dataset
RESULTS=/absolute/path/to/aspire-cami-mock-results

./examples/configure_mock_run.sh \
  --dataset "$DATASET" --output "$RESULTS" \
  --config-out aspire-cami-mock.yml --threads 32
./run_asv_pipeline.sh aspire-cami-mock.yml
./examples/validate_mock_run.sh --dataset "$DATASET" --results "$RESULTS"
```

The public name describes the dataset and benchmark. Existing completed runs
can be validated in their original directories; published provenance retains
the paths used for execution. The repository records the validated settings and
results in `examples/benchmarks/aspire-cami-mock/`.

## Analysis settings

Full-length taxonomy precedes control decontamination. Biological samples need
at least 5,000 post-QC ASV reads. Nonzero TECH and BIO controls enter their
respective prevalence tests. TECH uses extraction blanks and BIO uses Skin,
each against the same eligible Airways/Oral cohort. Both decontam score thresholds
are 0.1, and the union of flagged ASVs is removed. Local contaminant and
mitochondrial references provide the subsequent reference screen.

Before metadata construction, ASVs must meet **both** of these criteria:

- At least **0.1% relative abundance in one biological sample**.
- A nonzero count in at least **5% of biological samples**: seven of the 125
  retained samples in this run.

Patient-aware diversity, taxonomy, lung-side comparisons, power analysis,
indicators, VOC associations, network inference and topology run alongside the
other benchmarked modules. Genome linkage is enabled when separate genome
inputs are supplied. See [decontamination](decontamination.md) for exact
threshold definitions and [configuration](CONFIGURATION.md) for units.

## Validated results

The completed run was independently revalidated on **7 October 2026**. All
**29 validator checks passed**, and all 751 recorded Nextflow tasks completed
or were reused from cache. Raw and final counts also matched the earlier
completed run exactly when ASVs were matched by sequence identity.

The full rebuild changed five of 26,875 selected batch-corrected pseudocount
cells by one count. The largest absolute change in VOC Spearman rho was 0.00095
for sample-level correlations and 0.00510 for patient-level correlations; no
VOC significance calls changed. Network edges changed from 1,569 to 1,570 and
the largest module from 30 to 27 ASVs; recovery remained nine planted within-module
edges. The table below reports the refreshed run, and the provenance record
includes the count and VOC comparisons.

| Check or output | Observed result |
| --- | ---: |
| Published libraries | 179 / 179 |
| Eligible biological samples | 125 / 125 |
| Raw inferred ASVs | 568 |
| Raw ASV counts | 7,309,193 |
| TECH calls | 108 |
| BIO calls | 97 |
| TECH/BIO overlap | 15 |
| Unique contaminant ASVs removed | 190 |
| Final ASVs | 221 |
| Final biological counts | 3,912,205 |
| Non-target truth identified | 4 / 4 |
| ASVs with assigned taxonomy | 536 |
| Significant indicator entries | 141 |
| Planted indicator ASVs recovered | 23 |
| Planted positive VOC associations recovered | 7 |
| Patient-aware VOC tests | 2,340; 50 patients; zero insufficient tests |
| Network edges | 1,570 |
| Largest inferred module | 27 ASVs |
| Planted within-module edges recovered | 9 |
| Topology null draws | 100 |

The validator also checks control membership, low-depth control inclusion,
class-specific prevalence denominators, score-based calls, pairing, exact retained
counts, FASTA/table agreement, biological feature-filter audits, SVG outputs and
publication-manifest hashes. Patient-aware VOC validation checks valid inference
outputs; it does not require statistical significance. The network benchmark
requires at least three planted within-module edges, ten inferred edges and an
inferred module with at least two members.

Download the recorded evidence:

- {download}`Validation log <../examples/benchmarks/aspire-cami-mock/validation.txt>`
- {download}`Parameter record <../examples/benchmarks/aspire-cami-mock/parameters.yml>`
- {download}`Provenance and hashes <../examples/benchmarks/aspire-cami-mock/provenance.json>`

These files preserve the
check results, portable parameter record and hashes identifying the input and
output tables. The parameter record uses path placeholders; use the configurator
above to produce a runnable YAML. The manuscript benchmark measures recovery for
this fixed synthetic design. The [installation fixture](test.md) uses a separate,
smaller community with deliberately strong planted signals.

## Primer trimming in current runs

The configurator enables Cutadapt before fastp for both mock datasets, using the
515F/806R V4 primer family and a 64-base prefix search window to accommodate
the large mock's R2 adapter prefixes. Screening and trimming each allow a 0.1 substitution
rate; the minimum paired-primer screening support remains 50%. All four fixed
fastp clipping settings are zero. The validator additionally checks primer-family
membership, read-pair accounting and the handoff to fastp.

The recorded validation logs and numerical results on this page describe the
previous fixed-clipping run. They remain historical evidence; the Cutadapt
variation requires a new full run and validation. Re-running the configuration
command updates an existing YAML, and Nextflow reruns affected tasks.

### Power settings for new runs

Newly generated large-mock configurations use the production power settings:
1,000 simulations, 999 permutations, patient counts
`4,6,8,10,15,20,30,40,50`, and sample-type counts `10,15,20,25,30,40,50`.
The committed run record retains the settings used for that completed run.
The small quickstart uses 10 simulations, 49 permutations, and counts
`4,6,8,10` to keep installation checks short.
