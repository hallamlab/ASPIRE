from pathlib import Path


PIPELINE = (
    Path(__file__).parents[1]
    / "workflow"
    / "modules"
    / "indicators_diagnostics.nf"
)


def test_isa_plot_plan_is_bounded_to_configured_primary_groups():
    source = PIPELINE.read_text()
    section = source[
        source.index("process INDICSPECIES_PLOTS"):
        source.index("process INDICSPECIES_ALIGNED_PLOTS")
    ]
    assert "configured_groups = json.loads" in section
    assert 'summary_files = [Path(f"{group}{selected_suffix}") for group in configured_groups]' in section
    assert 'glob("*_indicator_species*_summary.tsv")' not in section
    assert "isa_plot_plan.tsv" in section


def test_isa_plot_output_replaces_old_tree_only_after_rendering():
    source = PIPELINE.read_text()
    section = source[
        source.index("process INDICSPECIES_PLOTS"):
        source.index("process INDICSPECIES_ALIGNED_PLOTS")
    ]
    render = section.index("subprocess.run(cmd, check=True)")
    replace = section.index("shutil.rmtree(final_out_root)")
    assert render < replace


def test_isa_plot_subprocess_reuses_active_conda_python():
    source = PIPELINE.read_text()
    section = source[
        source.index("process INDICSPECIES_PLOTS"):
        source.index("process INDICSPECIES_ALIGNED_PLOTS")
    ]
    assert "sys.executable" in section
    assert '"python",\n        str(plot_script)' not in section
