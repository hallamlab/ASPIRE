from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_default_phase_is_all_and_supported_phases_are_explicit() -> None:
    launcher = (ROOT / "run_asv_pipeline.sh").read_text(encoding="utf-8")
    assert 'ASPIRE_PHASE="all"' in launcher
    assert 'all|preprocess|analysis)' in launcher
    assert 'ASPIRE_PHASE="$ASPIRE_PHASE"' in launcher


def test_phase_resume_pointers_are_independent() -> None:
    launcher = (ROOT / "run_asv_pipeline.sh").read_text(encoding="utf-8")
    assert 'if [[ "$ASPIRE_PHASE" != "all" ]]; then' in launcher
    assert 'last_successful_run.${ASPIRE_PHASE}' in launcher
    assert 'last_successful_run.preprocess"' not in launcher
    assert 'last_successful_run.analysis"' not in launcher


def test_analysis_phase_uses_canonical_dataset_without_raw_reads() -> None:
    main = (ROOT / "asv_pipeline.nf").read_text(encoding="utf-8")
    preprocess_start = main.index("workflow ASPIRE_PREPROCESS")
    analysis_start = main.index("workflow ASPIRE_ANALYSIS")
    entry_start = main.index("workflow {", analysis_start)
    preprocess = main[preprocess_start:analysis_start]
    analysis = main[analysis_start:entry_start]
    assert "FASTP_QC" in preprocess
    assert "PREPROCESS_DATASET" in preprocess
    assert "FASTP_QC" not in analysis
    assert "rawReadsChannel" not in analysis
    assert "VALIDATE_PREPROCESS_DATASET" in analysis
    assert "preprocessing_dataset" in analysis
    assert "datasetFile('sample_manifest.tsv')" in analysis


def test_validation_deep_hashes_canonical_dataset_across_phase_paths() -> None:
    preprocessing_module = (
        ROOT / "workflow/modules/preprocessing.nf"
    ).read_text(encoding="utf-8")
    validation = preprocessing_module[
        preprocessing_module.index("process VALIDATE_PREPROCESS_DATASET"):
    ]
    assert "cache 'deep'" in validation


def test_analysis_manifest_is_an_explicit_process_input() -> None:
    metadata_module = (
        ROOT / "workflow/modules/metadata_ecology.nf"
    ).read_text(encoding="utf-8")
    assert metadata_module.count("path(sample_manifest)") == 2
    assert '${sample_manifest}' in metadata_module
    assert '${manifestPath}' not in metadata_module


def test_canonical_dataset_is_not_reorganized_as_a_module() -> None:
    organizer = (
        ROOT / "processes/output_layout/organize_outputs.py"
    ).read_text(encoding="utf-8")
    assert '"preprocessing_dataset",' in organizer
