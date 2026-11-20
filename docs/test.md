# Quickstart installation test

From the repository root, after [installation](installation.md), run the bundled 46-library fixture with 218,000 paired reads.

```bash
./examples/run_quickstart.sh --output "$PWD/aspire-quickstart" --threads 4
```

The wrapper generates the fixture, configures ASPIRE, processes the reads and checks the published results against planted signals. The first run downloads environments and taxonomy references. Repeat the command after an interruption to resume.

Look for `All mock-run checks passed.` in `aspire-quickstart/validation.txt`. Open `aspire-quickstart/results/summary/report/ASPIRE_run_report.html` to browse the outputs. The recorded Cutadapt quickstart passed 31 checks. {download}`Validation log <../examples/quickstart/validation.txt>` · {download}`Fixture identity and settings <../examples/quickstart/validation.json>`.

Continue to the [respiratory microbiome vignette](cami-mock.md) for a worked study.
