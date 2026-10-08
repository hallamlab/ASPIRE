# Quickstart and installation test

Run the bundled **aspire-quickstart** fixture to check your installation and try
ASPIRE's full benchmarked workflow. It has 218,000 paired reads in 46 libraries
and occupies approximately **2.5 MB** when generated. Sequence seeds and the
reproducible generator are included in the repository.

The larger **aspire-cami-mock** is the separate
[manuscript demonstration](https://hallamlab-aspire.readthedocs.io/en/latest/cami-mock.html),
with 50 synthetic patients and 179 libraries.

## Run the quickstart

From the ASPIRE repository, with Python 3 and Mamba available:

```bash
./examples/run_quickstart.sh --output "$PWD/aspire-quickstart" --threads 4
```

This builds the fixture, configures ASPIRE, runs the pipeline and executes the
**same full validator used for the manuscript demonstration**. Success ends with
`All mock-run checks passed.` The command retains a validation log and resumes
existing runs when repeated.

```text
aspire-quickstart/
├── dataset/                 # FASTQs, metadata, references, truth and checksums
├── run.yml                  # Machine-local pipeline configuration
├── run.manifest.tsv
├── run.metadata.tsv
├── results/                 # Published ASPIRE output and .aspire runtime cache
└── validation.txt
```

Open `aspire-quickstart/results/summary/report/ASPIRE_run_report.html` to inspect
sample accounting, plots, analysis outputs and execution logs. The
[output guide](https://hallamlab-aspire.readthedocs.io/en/latest/outputs.html)
explains how to view reports over SSH.

The first run installs analysis environments and downloads full SINA and
QIIME2/SILVA references. Those dependencies are much larger than the input
fixture. Runtime depends on environment and reference caches, hardware and the
analysis modules. The quickstart uses 10 power simulations and 49 permutations per power test to
exercise that module efficiently. These coarse power curves serve as execution
checks. The manuscript configuration uses 100 simulations and 199 permutations.
Both fixtures retain 999 indicator permutations, the full network settings and
identical validator thresholds.

## What the fixture exercises

Twelve synthetic patients provide 18 airway and 12 oral libraries, with 12 Skin
BIO controls and four extraction blanks. Clinical labels, paired lung-side
samples and two batches exercise patient-aware, clinical and batch modules.
The fixture uses 32 V4 sequences and deliberately strong body-site, microbial
network and chemical signals. FASTQs contain error-free amplicon reads with the
configured V4 primer prefixes. The generator's seed is fixed at 20261007.

Biological libraries contain 6,000 read pairs. Skin libraries contain 1,800;
blanks contain 800, 1,600, 6,000 and 8,000. Nonzero controls participate in their
own decontam arms. Only biological samples face the 5,000 post-QC ASV-read cutoff.
TECH uses blanks and BIO uses Skin; each score cutoff is 0.1. Their union removal
precedes reference screening and combined biological feature filtering.

The final feature gates require both **0.1% RA in at least one biological sample**
and **nonzero counts in at least 5% of biological samples**. For the quickstart's
30 biological samples, the prevalence gate requires two nonzero samples.
Retained counts and patient pairing feed the downstream analyses.

The full validator checks:

- Publication layout, sample accounting and output-manifest hashes.
- Control classes, biological depth QC, low-depth control inclusion, prevalence
  denominators, TECH/BIO scores, union removal, retained counts and pairing.
- Biological prevalence/RA audit results and FASTA/table agreement.
- Taxonomy, known non-target removal and sequence-based truth mapping.
- Planted indicators, positive VOC associations and patient-inference outputs.
- Network truth recovery, topology null draws and SVG outputs.

The network contract requires at least ten inferred edges, a module with two
members, and three edges between members of the same planted module. Failed
checks remain failures. Genome linkage uses separate genome inputs and is outside
both amplicon fixtures. The installation fixture demonstrates successful
execution on a small, designed input; manuscript recovery results are reported
for `aspire-cami-mock`.

## Validated quickstart result

The recorded run passed all **29 checks** on 7 October 2026: 30 biological samples,
24 final ASVs, two recovered planted indicators, all three planted positive VOC
associations, and 30 planted within-module network edges. It used eight threads
and the short power grid; the wrapper defaults to four threads. The
`examples/quickstart/validation.txt` and `validation.json` files record the checks
and fixture identity. Package/reference caches were reused for this validation.

## Run the steps separately

```bash
python3 examples/mock_test/build_quickstart.py --output "$PWD/aspire-quickstart/dataset"
./examples/configure_mock_run.sh \
  --dataset "$PWD/aspire-quickstart/dataset" \
  --output "$PWD/aspire-quickstart/results" \
  --config-out "$PWD/aspire-quickstart/run.yml" --threads 4
./run_asv_pipeline.sh "$PWD/aspire-quickstart/run.yml"
./examples/validate_mock_run.sh \
  --dataset "$PWD/aspire-quickstart/dataset" \
  --results "$PWD/aspire-quickstart/results"
```

Run the builder once into a new directory. The configurator checks input
checksums, writes absolute local paths and builds the controller environment as
needed. The launcher and validator reuse that environment. Change `--threads`
to suit the host; the quickstart wrapper defaults to four.

## Supplied datasets and validation contract

The same configurator and validator accept a supplied DECOI mock dataset. It
must contain paired FASTQs, `fastq_manifest.tsv`, `sample_metadata.tsv`,
`chemistry.tsv`, and these truth/reference files:

```text
ground_truth_feature_registry.tsv
ground_truth_group_effects.tsv
ground_truth_asv_chem.tsv
ground_truth_network_modules.tsv
ground_truth_extraction_controls.tsv
ground_truth_reference_filters.tsv
references/mitochondria.fasta
references/contaminants.fasta
```

Manifest IDs must match metadata and FASTQs. Required clinical metadata columns
are `sample_id`, `Participant_ID`, `Case`, `Type_Group`, `lung_status`, `batch`,
`DNA_conc`, `is_negative_control` and `is_positive_control`. CAMI-profile metadata
also records `source_sample`, `body_site` and `study_role`; `synthetic_clinical`
enables the patient and clinical analysis paths. Chemistry has `sample_id` and
numeric compound columns. Truth tables use stable ASV identifiers and V4
sequences so inferred features can be mapped back to the planted signals.

The optional `checksums.tsv` contains `relative_path` and `sha256` columns;
the configurator verifies it before writing a configuration. The full validation
contract is recorded in `examples/mock_test/expected_results.tsv`.

## Resume and diagnostics

Rerun the quickstart wrapper or `./run_asv_pipeline.sh aspire-quickstart/run.yml`
after correcting an interrupted download or environment solve. The retained
`.aspire/` directory supplies the resume cache. To rerun only the affected stages:

```bash
./run_asv_pipeline.sh --list-stages
./run_asv_pipeline.sh aspire-quickstart/run.yml --rerun-from standard:FILTER_ASVS
```

Changes to biological abundance/prevalence thresholds start at `FILTER_ASVS`.
Changes to control labels, biological depth or decontam score cutoffs start at
`core:CONTROL_DECONTAM`. Successful runs publish an atomic, consistent output
snapshot with configuration and task provenance.

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
