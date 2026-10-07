# aspire-cami-mock: manuscript run record

This record documents the completed 50-patient, 179-library demonstration.
The dedicated manuscript section is `docs/cami-mock.md`.

- `validation.txt`: all 29 checks passed on 7 October 2026. Only the local output
  directory in the log is replaced with the public benchmark name.
- `parameters.yml`: the published run configuration, with local filesystem paths
  replaced by explicit placeholders. Generate a runnable configuration with
  `examples/configure_mock_run.sh`.
- `provenance.json`: dataset/table/configuration hashes, sample and ASV counts,
  task status counts and the recorded Nextflow version.

The result tables and FASTQs remain in the completed run and dataset directories.
The provenance hashes identify those exact inputs and outputs independently of
how a local directory is named. For the compact installation fixture, use
`examples/run_quickstart.sh`.
