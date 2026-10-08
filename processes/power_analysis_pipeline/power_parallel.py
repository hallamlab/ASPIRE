"""Ordered, seed-stable simulation workers for the Linux/macOS ASPIRE runner.

Fork shares read-only input matrices without serializing them for each replicate.
Only the parent writes checkpoints; scientific functions own their indexed seeds.
"""
import hashlib
import importlib.metadata
import json
import multiprocessing as mp
import os
from pathlib import Path
import pickle
import time

_WORKERS = 1
_CHECKPOINT_DIR = None
_REPLICATE = None


def add_parallel_arguments(parser):
    parser.add_argument('--workers', type=int, default=1,
                        help='Maximum independent simulation processes')
    parser.add_argument('--checkpoint-dir', default=None,
                        help='Persistent replicate checkpoints (default: outdir/checkpoints)')


def configure_parallel(workers=1, checkpoint_dir=None):
    global _WORKERS, _CHECKPOINT_DIR
    if isinstance(workers, bool) or not isinstance(workers, int) or workers < 1:
        raise ValueError('workers must be a positive integer')
    _WORKERS, _CHECKPOINT_DIR = workers, checkpoint_dir


def _execute(index):
    return _REPLICATE(index)


def _json_default(value):
    if hasattr(value, 'tolist'):
        return value.tolist()
    raise TypeError(f'Unsupported checkpoint value: {type(value)}')


def run_replicates(function, count, context):
    """Yield results in replicate order, independent of worker count and resume.

    Context includes the simulation inputs/parameters before any replicate runs.
    Checkpoints are invalidated by inputs, parameters, source or package versions.
    """
    global _REPLICATE
    if count < 1:
        raise ValueError('simulation count must be positive')
    state = {k: v for k, v in context.items() if not callable(v)}
    digest = hashlib.sha256(pickle.dumps((function.__module__, function.__qualname__, count, state), protocol=5))
    for source in sorted(Path(__file__).parent.glob('*.py')):
        digest.update(source.name.encode())
        digest.update(source.read_bytes())
    for package in ('numpy', 'scipy', 'pandas', 'scikit-learn', 'statsmodels'):
        digest.update(importlib.metadata.version(package).encode())
    checkpoint = None
    completed = []
    if _CHECKPOINT_DIR:
        directory = Path(_CHECKPOINT_DIR)
        directory.mkdir(parents=True, exist_ok=True)
        checkpoint = directory / (digest.hexdigest() + '.json')
        if checkpoint.exists():
            completed = json.loads(checkpoint.read_text())
            if not isinstance(completed, list) or len(completed) > count:
                raise ValueError(f'Invalid checkpoint: {checkpoint}')
    def save():
        if checkpoint:
            temporary = checkpoint.with_suffix(f'.{os.getpid()}.tmp')
            temporary.write_text(json.dumps(completed, default=_json_default))
            temporary.replace(checkpoint)
    start = time.monotonic()
    label = function.__qualname__.split('.')[0]
    print(f'    {label}: {len(completed)}/{count} cached; workers={_WORKERS}', flush=True)
    for result in completed:
        yield result
    remaining = range(len(completed), count)
    if not remaining:
        return
    pool = None
    _REPLICATE = function
    try:
        if _WORKERS > 1:
            if 'fork' not in mp.get_all_start_methods():
                raise RuntimeError('Parallel power analysis requires a POSIX fork platform; use workers=1')
            pool = mp.get_context('fork').Pool(min(_WORKERS, len(remaining)))
            results = pool.imap(_execute, remaining, chunksize=1)
        else:
            results = map(function, remaining)
        for result in results:
            completed.append(result)
            if len(completed) % 25 == 0 or len(completed) == count:
                save()
                print(f'    {label}: {len(completed)}/{count}; elapsed {time.monotonic()-start:.1f}s', flush=True)
            yield result
    finally:
        save()
        if pool:
            pool.terminate()
            pool.join()
        _REPLICATE = None
