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
- Host taxa still appear downstream: confirm `standard.filter_counts.exclude_taxa` is set and rerun from `standard:FILTER_ASVS` or at least from `standard:PLOT_METADATA` if the final `ASV_target.tsv` is already corrected.
- VOC direction did not change outputs: rerun from `VOC_CORRELATION`.
- Sankey complains about intermediates: set `standard.filter_counts.save_intermediates: true`.
- BLAST database errors: set `standard.mito.mito_db` and `standard.mito.biof_db` to valid database prefixes or compatible FASTA paths.
- Conda solve errors: confirm `mamba` is available and review the relevant `environments.*` config entry.

## Full mock stops at CLUSTERMAPS with exit 137

Exit 137 with `Killed` means the process received SIGKILL. Memory pressure is one possible cause. Earlier versions exported even the largest clustermaps at 600 DPI. Current versions scale the raster resolution to the figure dimensions and preserve vector cells and text in PDF and SVG.

For the vignette, update the checkout and rerun the same configuration. Keep the existing results and hidden runtime directory so Nextflow can reuse completed tasks.

```bash
cd ~/aspire-vignette/ASPIRE
git pull --ff-only origin main
./run_asv_pipeline.sh ../cami.yml
```

If the updated process is killed again, check the kernel or scheduler memory report for the time of the failure. It can distinguish an out-of-memory kill from another external termination.

## MITOMASTER API availability

Requests have bounded retries and timeouts. `standard.mito.mitomaster_failure_policy: fail` stops after an API failure; `continue` preserves successful responses and allows local BLAST screening to finish even if all requests fail. Inspect the published status JSON and failure TSV before interpreting mitochondrial evidence. See the [configuration reference](CONFIGURATION.md).

## Unavailable outlier classifications

A detector that fails or cannot form valid clusters is unavailable. If too few detectors remain for the vote threshold, consensus is unknown. Review detector availability and group diagnostics rather than treating unknown calls as inliers.

## Taxonomy stops while importing `Self` from `typing`

A traceback ending in `flufl/lock/_lockfile.py` during `core:TAXONOMY` indicates an incompatible locking-library build in the QIIME 2 Python 3.10 environment. ASPIRE main pins the compatible `flufl.lock=9.2.0` release. From the checkout used for the failed quickstart, update main and repeat the original quickstart command with the same output directory and thread count:

```bash
git pull --ff-only origin main
./examples/run_quickstart.sh --output "$PWD/aspire-quickstart" --threads 4
```

Keep the existing results and hidden runtime directory so Nextflow can reuse completed tasks and build the updated taxonomy environment.
