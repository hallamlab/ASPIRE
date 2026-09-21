from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_cache_management_is_dry_run_first_and_phase_safe() -> None:
    launcher = (ROOT / "run_asv_pipeline.sh").read_text(encoding="utf-8")
    assert "cache prune [--keep N] [--force]" in launcher
    assert 'CACHE_FORCE=0' in launcher
    assert "cache_pointer_runs" in launcher
    assert "cache_managed_runs" in launcher
    assert 'index($7, expected_work) > 0' in launcher
    assert "adopt_legacy_cache_baseline" in launcher
    assert "legacy_successful_pipeline_run" in launcher
    assert "valid legacy successful ASPIRE run exists" in launcher
    assert 'command !~ /(^|[[:space:]])-preview([[:space:]]|$)/' in launcher
    assert 'nextflow clean "$run_name" -f -q' in launcher


def test_cache_usage_and_explicit_conda_clear_are_exposed() -> None:
    launcher = (ROOT / "run_asv_pipeline.sh").read_text(encoding="utf-8")
    assert "cache usage" in launcher
    assert "--include-conda" in launcher
    assert "nextflow_work" in launcher
    assert "conda_cache" in launcher
    assert "publication_staging" in launcher


def test_post_run_pruning_is_opt_in() -> None:
    launcher = (ROOT / "run_asv_pipeline.sh").read_text(encoding="utf-8")
    assert 'PRUNE_CACHE_AFTER_RUN=0' in launcher
    assert "--prune-cache-after-run" in launcher
    assert 'if [[ "$PRUNE_CACHE_AFTER_RUN" -gt 0 ]]' in launcher


def test_console_process_names_are_shortened_without_renaming_workflows() -> None:
    launcher = (ROOT / "run_asv_pipeline.sh").read_text(encoding="utf-8")
    pipeline = (ROOT / "asv_pipeline.nf").read_text(encoding="utf-8")
    assert "format_nextflow_console" in launcher
    assert "s/RUN_METADATA_ANALYSES:/COMMUNITY:/g" in launcher
    assert "s/ASPIRE_ANALYSIS:/ANALYSIS:/g" in launcher
    assert "workflow RUN_METADATA_ANALYSES" in pipeline
    assert "workflow ASPIRE_ANALYSIS" in pipeline
