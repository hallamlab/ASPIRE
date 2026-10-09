# ASPIRE quickstart fixture

`features.tsv` contains 32 V4 seed sequences, their original DECOI identifiers,
sequence hashes and taxonomy. `source.json` records the source registry hash and
selection. The source sequences come from the larger CAMI2-derived fixture.
The deterministic builder creates a separate synthetic community, counts,
metadata, chemistry, truth tables and error-free paired FASTQs.

From the repository root:

```bash
./examples/run_quickstart.sh --output "$PWD/aspire-quickstart" --threads 4
```

The fixture has 12 synthetic patients, 30 biological libraries, 12 Skin BIO
controls and four TECH blanks. Biological libraries contain 6,000 read pairs,
Skin libraries 1,800, and blanks 800, 1,600, 6,000 and 8,000. Total: 218,000 pairs.
There are 24 biological ASVs, four BIO-control ASVs, three TECH-control ASVs and
one mitochondrial fixture. Three five-ASV modules share strong latent abundance
signals; six ASVs have body-site effects; three chemical measurements track their
specified ASV drivers. Sparse control carry-over exercises TECH/BIO removal.

The seed is fixed at 20261007. FASTQ gzip timestamps are fixed, and the generated
`checksums.tsv` covers all fixture files. Generation requires only Python 3's
standard library. The input directory occupies approximately 2.5 MB. Package
installation, full taxonomy references and analysis outputs require additional
space; cached installations reuse those dependencies.

This fixture tests installation and execution with deliberately strong signals.
Use the separate `aspire-cami-mock` demonstration for a larger synthetic study.
The same full mock-run validator checks both fixtures, with the same thresholds.

The configurator recognizes `name: aspire-quickstart` in the dataset manifest and
sets the power grid to 10 simulations and 49 permutations per test. Every power
analysis path is retained. Other analysis settings and all truth-recovery checks
are shared with the larger-demonstration configuration, which uses 1,000
simulations, 999 permutations and the expanded patient-count grids.

`validation.txt` and `validation.json` record the successful pre-Cutadapt 29-check validation,
including the fixture checksum, generator hash and run settings. The recorded run
used eight threads and retained 24 ASVs in 30 biological samples. It recovered
two planted indicators, three positive VOC pairs and 30 planted network edges.

The current generated configuration enables Cutadapt and its additional primer audit checks. Re-run and validate this configuration to obtain current-run evidence.
