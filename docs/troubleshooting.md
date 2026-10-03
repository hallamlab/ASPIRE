# Troubleshooting

For task-level command and log inspection, use the [expert guide](EXPERT_GUIDE.md#reading-a-failed-task). Fix the reported failure, then use the same wrapper command to resume.

## Troubleshooting

If environment creation fails while collecting repodata with
`BlockingIOError: [Errno 11] Resource temporarily unavailable`, the failure is
a Conda cache-lock collision rather than evidence of an unavailable package.
Current wrapper runs use a run-local package cache and serialize environment
creation through `flock`. Restart the failed run normally; its retained
Nextflow work directory remains resumable. DNS, connection-timeout, or HTTP
errors instead indicate network or repository availability and are retried by
the controller before the run fails.

Interrupted Nextflow sessions can also leave `.env-*.lock` marker files even
when no `mamba` process remains. The wrapper takes an exclusive lock on the
configured Conda cache, rejects concurrent runs that share that cache, and then
removes these orphaned markers automatically. Messages about waiting for the
serialized mamba slot indicate active environment creation, not a deadlock.

- `No usable entries detected in manifest`: check tab separation, sample IDs, and FASTQ paths.
- `metadata file not found`: set the branch-specific metadata path or disable that branch.
- Host taxa still appear downstream: confirm `standard.filter_counts.exclude_taxa` is set and rerun from `standard:FILTER_COUNTS` or at least from `standard:PLOT_METADATA` if the final `ASV_target.tsv` is already corrected.
- VOC direction did not change outputs: rerun from `VOC_CORRELATION`.
- Sankey complains about intermediates: set `standard.filter_counts.save_intermediates: true`.
- BLAST database errors: set `standard.mito.mito_db` and `standard.mito.biof_db` to valid database prefixes or compatible FASTA paths.
- Conda solve errors: confirm `mamba` is available and review the relevant `environments.*` config entry.
