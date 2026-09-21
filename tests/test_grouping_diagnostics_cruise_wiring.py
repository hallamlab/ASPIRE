from pathlib import Path


PROJECT = Path(__file__).parents[1]


def _yaml_block(text: str, key: str) -> str:
    marker = f"{key}:\n"
    start = text.index(marker) + len(marker)
    remainder = text[start:]
    lines = []
    for line in remainder.splitlines():
        if line and not line.startswith((" ", "\t", "#")):
            break
        lines.append(line)
    return "\n".join(lines)


def test_active_grouping_diagnostics_are_limited_to_depth_hybrid_and_mc():
    config = (PROJECT / "examples" / "si.local.yml").read_text()
    block = _yaml_block(config, "grouping_diagnostics")
    assert "- Depth" in block
    assert "- o2_subcompartment_final" in block
    assert "- microbial_compartment" in block
    assert "cruise_level_group_cols: []" in block
    assert "require_cruise_level_groups: false" in block
    assert "primary_group: microbial_compartment" in block


def test_si_analyses_prioritize_hybrid_and_mc_with_focused_renewal_posthoc():
    config = (PROJECT / "examples" / "si.local.yml").read_text()
    isa_block = _yaml_block(config, "indicspecies")
    assert "renewal_phase" not in isa_block
    assert "- o2_subcompartment_final" in isa_block
    assert "- microbial_compartment" in isa_block
    assert "within_col: o2_subcompartment_final" in isa_block
    assert "group_col: microbial_compartment" in isa_block
    assert "within_col: microbial_compartment" in isa_block
    assert "group_col: o2_subcompartment_final" in isa_block

    master_summary_block = _yaml_block(config, "master_summary")
    assert "microbial_compartment_indicator_species_summary.tsv" in master_summary_block
    assert "stratified_microbial_compartment_within_o2_subcompartment_final" in master_summary_block
    assert "stratified_o2_subcompartment_final_within_microbial_compartment" in master_summary_block

    state_block = _yaml_block(config, "microbial_state_interpretation")
    assert "renewal_col: renewal_phase" in state_block
    assert "renewal_levels: [baseline, renewal, post-renewal]" in state_block
    assert "- o2_subcompartment_final" in state_block

    spieceasi_block = _yaml_block(config, "spieceasi")
    assert "- o2_subcompartment_final" in spieceasi_block
    assert "- microbial_compartment" in spieceasi_block
    assert "- renewal_phase" not in spieceasi_block


def test_wrapper_and_renderer_enforce_required_cruise_groups():
    wrapper = (PROJECT / "run_asv_pipeline.sh").read_text()
    renderer = (PROJECT / "processes" / "grouping_diagnostics" / "grouping_diagnostics.py").read_text()
    assert "GROUPING_REQUIRE_CRUISE_LEVEL" in wrapper
    assert "GROUPING_CRUISE_LEVEL_COUNT" in wrapper
    assert '--require-cruise-level-groups' in renderer
    assert "missing configured metadata columns" in renderer
    assert "grouping_hierarchical_within_cruise_states.tsv" in renderer
    assert "grouping_hierarchical_across_cruise_states.tsv" in renderer
    assert "grouping_hierarchical_interactions.tsv" in renderer
    assert "restricted_label_permutation" in renderer
    assert '"whole_cruise"' in renderer


def test_wrapper_uses_one_standard_nextflow_rerun_path():
    wrapper = (PROJECT / "run_asv_pipeline.sh").read_text()
    workflow = (PROJECT / "asv_pipeline.nf").read_text()
    assert "--rerun-from PROCESS_NAME" in wrapper
    assert "cache = false" not in wrapper
    assert "--start-from" not in wrapper
    assert "START_FROM_GROUPING_DIAGNOSTICS" not in workflow
