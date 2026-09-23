# Mock Dataset Test

This guide runs every ASPIRE analysis module against a separately supplied
mock dataset.

## Obtain the Mock Dataset

Download the ASPIRE mock dataset from Zenodo:

- DOI: [10.5281/zenodo.22906294](https://doi.org/10.5281/zenodo.22906294)
- Zenodo record: [https://zenodo.org/records/22906294](https://zenodo.org/records/22906294)

Extract the downloaded archive and use the extracted `mock_dataset/` directory
as the `--dataset` argument in the commands below.

## Reviewer Fast Path

From a clone of the `main` branch, edit the first two paths and run the block:

```bash
DATASET=/absolute/path/to/mock_dataset
RESULTS=/absolute/path/to/aspire_mock_output

./examples/configure_mock_run.sh \
  --dataset "$DATASET" \
  --output "$RESULTS" \
  --config-out mock_run.generated.yml
./run_asv_pipeline.sh mock_run.generated.yml --no-resume
./examples/validate_mock_run.sh --dataset "$DATASET" --results "$RESULTS"
```

Success ends with `All mock-run checks passed.` Continue below for the dataset
contract, expected outputs, interpretation, resume behavior, and diagnostics.

## Required Dataset Layout

```text
mock_dataset/
├── fastq/
│   ├── SAMPLE_001_R1.fastq.gz
│   └── SAMPLE_001_R2.fastq.gz
├── fastq_manifest.tsv
├── sample_metadata.tsv
├── chemistry.tsv
├── checksums.tsv                  # recommended
├── ground_truth_feature_registry.tsv
├── ground_truth_group_effects.tsv
├── ground_truth_asv_chem.tsv
├── ground_truth_network_modules.tsv
├── ground_truth_extraction_controls.tsv
├── ground_truth_reference_filters.tsv
└── references/
    ├── mitochondria.fasta
    └── contaminants.fasta
```

Additional truth tables may be included for interpretation. All files above
except the recommended `checksums.tsv` are required for configuration and
validation. The manifest sample
IDs must exactly match `sample_metadata.tsv` and the FASTQ filenames.
When `checksums.tsv` is present, the generator verifies it before writing a
configuration. It contains `relative_path` and `sha256` columns with paths
relative to `mock_dataset/`.

Required metadata columns for the committed mock configuration are
`sample_id`, `Participant_ID`, `Case`, `Type_Group`, `lung_status`, `batch`,
`DNA_conc`, `is_negative_control`, and `is_positive_control`.
The chemistry table is tab-separated, has one row per `sample_id`, and contains
one or more numeric VOC columns. Truth tables retain stable ASV IDs;
the validator maps ASPIRE's inferred ASV numbering back to those IDs by sequence.

The current control-bearing mock dataset contains 60 biological samples from 20 participants:
10 control participants contribute bronchial-brush and BAL samples (20 samples),
and 10 cancer participants contribute tumor-side and contralateral samples for
both types (40 samples). Four explicitly labelled synthetic extraction controls
exercise the three-tier decontamination process and are removed before the
biological downstream analyses.

## Configure

From the ASPIRE repository root:
The user does not create or activate a controller environment. This command
automatically creates or updates `<ASPIRE>/.controller_env` from the committed
controller environment definition, runs script inside it, and leaves it
in place. `run_asv_pipeline.sh` and `validate_mock_run.sh` automatically reuse
the same environment for the rest of the test. Only `mamba` must already be
available on `PATH`.

```bash
./examples/configure_mock_run.sh \
  --dataset /absolute/path/to/mock_dataset \
  --output /absolute/path/to/aspire_mock_output \
  --config-out mock_run.generated.yml
```

The script checks the fixture, writes absolute paths for the local machine,
sets the persistent runtime directory to `<output>/.aspire`, and prints the exact run command. The committed
`examples/mock.local.yml` is the template; users should not edit its
developer-specific paths manually. The generated configuration retains the
runtime directory so interrupted and completed runs remain resumable. It also
preserves the explicit `core`, `standard`, and `optional` namespaces so a
reviewer can distinguish required construction stages from selectable analyses.

## Run

```bash
./run_asv_pipeline.sh mock_run.generated.yml --no-resume
```

The first run downloads the configured SINA and QIIME2/SILVA references and
solves per-process Conda environments. Internet access is therefore required
unless those references and package caches were provided separately. The mock
configuration disables the external MITOMASTER service and uses the mock dataset's local mitochondrial and contaminant FASTA files.

Runtime depends on hardware, package-cache state, and FASTQ depth. A cold run
can take some time because references and process environments are prepared
before all 64 input samples (60 biological samples plus four controls) are
analyzed. The configured output contains the hidden
`.aspire/` resume cache plus `modules/`, `intermediates/`, `references/`,
`summary/`, and `logs/`. The non-interpretive inventory and provenance report at
`summary/report/ASPIRE_run_report.html` begins with sample, group, sequence, and
ASV accounting plus the Sankey, swarmplot, UpSet, and collector's-curve summaries.
It then inventories module outputs and links the Nextflow execution report, timeline, trace table,
workflow DAG, launch command, runtime version, controller log, and per-task logs.

## Validate

```bash
./examples/validate_mock_run.sh \
  --dataset /absolute/path/to/mock_dataset \
  --results /absolute/path/to/aspire_mock_output
```

The validator checks sample accounting, output structure, checksums,
taxonomy output, known non-target removal, indicator recovery, positive VOC
associations, a non-degenerate microbial network, network topology with a
degree-preserving null ensemble, and SVG production. Its
machine-readable contract is `examples/mock_test/expected_results.tsv`.

A zero exit status and `All mock-run checks passed.` indicate that the test completed the expected end-to-end benchmark. A failed statistical
signal check should be investigated rather than bypassed by relaxing ASPIRE's
analysis thresholds.

The public fixture contains synthetic extraction blanks and synthetic DNA
concentrations generated by DECOI. These values are benchmark inputs only and
carry no patient information.

## Resume And Diagnostics

After resolving an interrupted external download or environment solve, rerun:

```bash
./run_asv_pipeline.sh mock_run.generated.yml
```

ASPIRE resumes from the runtime cache retained after a failed or interrupted
run. For a targeted rerun, list valid stages and select the earliest affected
stage:

```bash
./run_asv_pipeline.sh --list-stages
./run_asv_pipeline.sh mock_run.generated.yml --rerun-from standard:FILTER_COUNTS
```

Targeted reruns remain available after successful runs because
`core.paths.keep_runtime_dir` defaults to `true`. Set it to `false` only when you
explicitly want to discard the work and Conda caches after publication.

The public output is replaced atomically only after a successful workflow, so a
failed rerun does not merge partial files into the previous completed result.
