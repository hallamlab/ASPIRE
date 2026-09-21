from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]


def test_every_modular_process_references_scoped_cache_generation():
    modules = list((ROOT / "workflow" / "modules").glob("*.nf"))
    source = "\n".join(path.read_text() for path in modules)
    controller = (ROOT / "run_asv_pipeline.sh").read_text()
    process_count = len(re.findall(r"^process\s+\w+\s*\{", source, re.MULTILINE))
    hook_count = source.count(
        "def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0"
    )
    assert process_count >= 48
    assert hook_count == process_count
    assert "ext.aspire_cache_generation" in controller


def test_main_is_orchestration_only_and_controller_composes_all_modules():
    main = (ROOT / "asv_pipeline.nf").read_text()
    controller = (ROOT / "run_asv_pipeline.sh").read_text()
    assert not re.search(r"^process\s+\w+\s*\{", main, re.MULTILINE)
    for module in (
        "preprocessing.nf",
        "taxonomy_decontamination.nf",
        "metadata_ecology.nf",
        "indicators_diagnostics.nf",
        "networks_reporting.nf",
    ):
        assert f'workflow/modules/{module}' in controller


def test_metadata_and_downstream_python_use_activated_conda_interpreter():
    # These stages can be entered directly with --rerun-from while ASPIRE
    # itself is launched under an outer `conda run`.
    modules = [
        ROOT / "workflow" / "modules" / name
        for name in (
            "metadata_ecology.nf",
            "indicators_diagnostics.nf",
            "networks_reporting.nf",
        )
    ]
    source = "\n".join(path.read_text() for path in modules)
    # A surrounding `conda run` can leave an outer Python ahead of the
    # Nextflow process environment on PATH. Resolve Python from CONDA_PREFIX so
    # imports always use the dependencies declared by the process environment.
    assert not re.search(r"^\s*python\s", source, re.MULTILINE)
    assert '"\\${CONDA_PREFIX}/bin/python"' in source


def test_controller_uses_persistent_generations_not_workdir_invalidation():
    source = (ROOT / "run_asv_pipeline.sh").read_text()
    assert 'CACHE_GENERATION_FILE="${RUNTIME_DIR}/cache_generations.tsv"' in source
    assert 'RERUN_CONFIG="${RUNTIME_DIR}/cache_generations.config"' in source
    assert "candidate_workdirs" not in source
    assert 'rm -rf "$workdir_real"' not in source


def test_generation_is_advanced_only_for_explicit_target():
    source = (ROOT / "run_asv_pipeline.sh").read_text()
    assert 'CACHE_GENERATIONS["$RERUN_FROM_CANONICAL"]=$((' in source
    assert "RERUN_STAGE_SET" not in source
    assert "true dependency descendants" in source


def test_clustermaps_hashes_selected_isa_table_as_a_path_input():
    main = (ROOT / "asv_pipeline.nf").read_text()
    module = (ROOT / "workflow" / "modules" / "networks_reporting.nf").read_text()

    assert "path(isa_table)" in module
    assert "val(indicspecies_ready)" not in module
    assert "clustermapsIsaSearchDir" not in main
    assert "clustermapsIsaSearchDir" not in module
    assert "indicspecies_stage.all_tables.collect().map" in main
    assert "indicspeciesIsaForClustermaps" in main
    assert "indicspeciesReadyForClustermaps" not in main
    assert "ISA_SOURCE_DIR" not in module


def test_si_network_filter_and_review_defaults_are_explicit():
    config = (ROOT / "examples" / "si.local.yml").read_text()
    assert "min_rel_abund: 0.005" in config
    assert "min_prevalence: 0.005" in config
    assert "force_keep_isa_asvs: false" in config
    assert "module_best_top_n: 8" in config
    assert "degree_size_mode: linear" in config
    assert "degree_scale: 3.5" in config


def test_network_and_guild_review_outputs_are_wired():
    main = (ROOT / "asv_pipeline.nf").read_text()
    module = (ROOT / "workflow" / "modules" / "networks_reporting.nf").read_text()
    graph = (ROOT / "processes" / "graph_network" / "graph_network.py").read_text()
    guild = (
        ROOT / "processes" / "group_guild_function" / "group_guild_function.py"
    ).read_text()

    assert "--module-best-top-n" in module
    assert "--graph" in module
    assert "--max-modules" in module
    guild_wrapper = main[main.index("workflow RUN_GROUP_GUILD_FUNCTION"):]
    assert re.search(
        r"take:\s+metadata_table\s+asv_counts\s+modules_all\s+node_features\s+isa_tables",
        guild_wrapper,
    )
    assert "weighted_modularity_contribution" in graph
    assert "module_quality_score" in graph
    assert graph.index("G_sub = load_graph(graph_sub)") < graph.index("graph=G_sub")
    assert graph.index("G_all = load_graph(graph_all)") < graph.index("graph=G_all")
    assert 'args.group1_summary or not args.isa_group_cols' in graph
    assert "ecological_module_asv_pca_overlay_asv_level" in guild
    assert "ecological_module_asv_pca_overlay_asv_level_abundant_asvs_labeled" in guild
    assert "ecological_module_asv_membership_all.tsv" in guild
    assert "ecological_module_asv_membership_mag_paired.tsv" in guild


def test_asv_clustermap_plot_can_be_skipped_without_dropping_table():
    config = (ROOT / "examples" / "si.local.yml").read_text()
    main = (ROOT / "asv_pipeline.nf").read_text()
    module = (ROOT / "workflow" / "modules" / "networks_reporting.nf").read_text()
    script = (ROOT / "processes" / "clustermaps" / "plot_clustermaps.py").read_text()

    assert "plot_asv_level: false" in config
    assert "clustermapsPlotAsvLevel" in main
    assert "--skip-asv-plot" in module
    assert "args.skip_asv_plot and rank == args.asv_id_col" in script
    assert 'pivot.to_csv(outdir / f"clustermap_{colname}.tsv"' in script


def test_corrected_asv_metadata_runs_for_every_enabled_consumer():
    module = (ROOT / "workflow" / "modules" / "metadata_ecology.nf").read_text()
    block = module[
        module.index("process ASV_META_FROM_CORRECTED"):
        module.index("process OUTLIER_CHECKER")
    ]
    for consumer_flag in (
        "bubbleplotterEnabled",
        "umapClusteringEnabled",
        "clustermapsEnabled",
        "vocCorrelationEnabled",
        "measurementAssociationEnabled",
        "masterSummaryEnabled",
        "powerAnalysisEnabled",
        "taxonomyPatientAwareEnabled",
        "lungStatusAnalysisEnabled",
    ):
        assert consumer_flag in block
