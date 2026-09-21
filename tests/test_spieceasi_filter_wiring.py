from pathlib import Path


ROOT = Path(__file__).parents[1]


def test_isa_force_retention_is_explicit_and_disabled_for_si():
    pipeline = (ROOT / "asv_pipeline.nf").read_text()
    network_module = (ROOT / "workflow" / "modules" / "networks_reporting.nf").read_text()
    config = (ROOT / "examples" / "si.local.yml").read_text()
    assert "spieceasiForceKeepIsaAsvs" in pipeline
    assert "indicspeciesEnabled && spieceasiForceKeepIsaAsvs" in network_module
    assert "force_keep_isa_asvs: false" in config


def test_asv_mag_force_retention_is_explicit_and_enabled_for_si():
    pipeline = (ROOT / "asv_pipeline.nf").read_text()
    network_module = (ROOT / "workflow" / "modules" / "networks_reporting.nf").read_text()
    config = (ROOT / "examples" / "si.local.yml").read_text()
    script = (ROOT / "processes" / "spieceasi" / "run_spieceasi.R").read_text()
    assert "spieceasiForceKeepAsvMagAsvs" in pipeline
    assert "asv_mag_link_stage.pairing" in pipeline
    assert "--force-keep-asv-mag-links" in network_module
    assert "force_keep_asv_mag_asvs: true" in config
    assert "read_force_keep_asv_mag_asvs" in script
    assert "retained_for_asv_mag_link" in script


def test_si_network_filters_are_the_requested_values():
    config = (ROOT / "examples" / "si.local.yml").read_text()
    block = config[config.index("spieceasi:"):config.index("asv_mag_link:")]
    assert "min_rel_abund: 0.005" in block
    assert "min_prevalence: 0.005" in block
    assert "lambda_min_ratio: 0.01" in block
    assert "thresh: 0.05" in block
    assert "edge_threshold: 0.0" in block


def test_spieceasi_standard_filter_is_a_measurement_association_dependency():
    pipeline = (ROOT / "asv_pipeline.nf").read_text()
    network_module = (ROOT / "workflow" / "modules" / "networks_reporting.nf").read_text()
    measurement_module = (
        ROOT / "workflow" / "modules" / "indicators_diagnostics.nf"
    ).read_text()
    measurement_script = (
        ROOT
        / "processes"
        / "measurement_association"
        / "measurement_association.py"
    ).read_text()
    config = (ROOT / "examples" / "si.local.yml").read_text()

    assert 'path("spieceasi_filtering_audit.csv"), emit: filter_audit' in network_module
    assert 'ln -sf "\\${FILTER_AUDIT}" spieceasi_filtering_audit.csv' in network_module
    assert "path(asv_subset_audit)" in measurement_module
    assert '--asv-subset-source "${measurementAssociationSubsetSource}"' in measurement_module
    assert "spieceasi_stage.filter_audit.map" in pipeline
    assert pipeline.index("spieceasi_stage = SPIECEASI(") < pipeline.index(
        "MEASUREMENT_ASSOCIATION("
    )
    assert "passes_standard_filter" in measurement_script
    assert "retained_final" not in measurement_script
    assert "asv_subset_source: spieceasi_standard_filter" in config
    assert "max_asvs: 0" in config

    controller = (ROOT / "run_asv_pipeline.sh").read_text()
    order = controller[
        controller.index("PROCESS_ORDER=(") :
        controller.index("declare -A PROCESS_ALIASES")
    ]
    assert order.index("ASV_MAG_LINK") < order.index("SPIECEASI")
    assert order.index("SPIECEASI") < order.index("MEASUREMENT_ASSOCIATION")


def test_spieceasi_workers_disable_nested_numeric_threads():
    section = (ROOT / "workflow" / "modules" / "networks_reporting.nf").read_text()
    for variable in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
        "RCPP_PARALLEL_NUM_THREADS",
    ):
        assert f"export {variable}=1" in section


def test_module_subnetworks_use_auditable_consortium_prominence_encoding():
    pipeline = (ROOT / "asv_pipeline.nf").read_text()
    workflow = (ROOT / "workflow" / "modules" / "networks_reporting.nf").read_text()
    script = (ROOT / "processes" / "graph_network" / "graph_network.py").read_text()
    config = (ROOT / "examples" / "si.local.yml").read_text()

    assert "networkModuleSubnetworkProminenceMetrics" in pipeline
    assert "--module-subnetwork-prominence-metrics" in workflow
    assert "--module-subnetwork-prominence-threshold" in workflow
    for metric in ("max_relative_abundance", "eigenvector", "participation"):
        assert f"    - {metric}" in config
    assert 'metrics["consortium_prominence_score"]' in script
    assert '"consortium_prominence_rank"' in script
    assert 'for rank, attr in (("phylum", "Phylum"),):' in script
    assert "def pfg_taxonomy" in script
    assert "network_ecological_module_renewal_phase_POS_ALL" in script
    assert "network_ecological_module_renewal_post_renewal_vs_stagnation_POS_ALL" in script
    assert "Post-renewal vs stagnation" in script
    assert '"is_high_prominence_asv"' in script
    assert 'module_high_prominence_asvs.tsv' in script
    assert 'set(mag_nodes + high_prominence_nodes + top_nodes)' in script
    assert '"plot_label_number"' in script
    assert "Top-ranked prominence ASVs" in script
    assert "High-prominence ASVs (score >=" in script
    assert "Genome-linked ASVs (bold border)" in script
    assert "role_band_heights" in script
    assert "asv_key_axes = [fig.add_subplot(grid[index, :]) for index in range(1, 4)]" in script
    assert "asv_key_axes, role_specs, role_column_counts" in script
    assert "uniform module zoom without node displacement" in script
    assert "external_label_positions" in script
    assert "ordered exterior lanes" in script
    assert "module_subnetwork_prominence_min_area: 30.0" in config
    assert "module_subnetwork_prominence_max_area: 520.0" in config
    assert 'CompleteTaxonomy' in script
