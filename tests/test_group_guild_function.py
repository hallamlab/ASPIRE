import importlib.util
from pathlib import Path

import pandas as pd


SCRIPT = Path(__file__).parents[1] / "processes" / "group_guild_function" / "group_guild_function.py"
SPEC = importlib.util.spec_from_file_location("group_guild_function", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_evidence_chain_retains_all_eligible_mags():
    isa = pd.DataFrame({
        "ASV_ID": ["ASV1"], "grouping": ["cruise_group"], "group_level": ["G1"],
        "isa_stat": [0.8], "isa_p_value": [0.01], "isa_q_value": [0.02],
    })
    modules = pd.DataFrame({"Taxon": ["ASV1"], "module_label": ["Module 1"]})
    nodes = pd.DataFrame({"Taxon": ["ASV1"], "Degree": [3]})
    mappings = pd.DataFrame([
        {"ASV_ID": "ASV1", "genome_id": "MAG1", "mag_match_id": "MAG1", "eligible_mag_count": 2, "asv_link_weight": 0.5},
        {"ASV_ID": "ASV1", "genome_id": "MAG2", "mag_match_id": "MAG2", "eligible_mag_count": 2, "asv_link_weight": 0.5},
    ])
    functions = pd.DataFrame([
        {"mag_match_id": "MAG1", "module_id": "M00529", "module_name": "Denitrification", "fraction_covered": 1.0, "present": True},
        {"mag_match_id": "MAG2", "module_id": "M00529", "module_name": "Denitrification", "fraction_covered": 0.2, "present": False},
    ])

    result = MODULE.evidence_chain(isa, modules, nodes, mappings, functions, ["M00529"])

    assert set(result["mag_match_id"]) == {"MAG1", "MAG2"}
    assert set(result["asv_link_weight"]) == {0.5}


def test_load_isa_retains_stratified_provenance(tmp_path):
    pd.DataFrame({
        "ASV": ["ASV1"], "s.oxic__gmm0": [1], "stat": [0.8], "q.value": [0.02],
        "stratified_within_col": ["cruise_group"],
        "stratified_within_value": ["cruise_group_1"],
        "stratified_group_col": ["o2_subcompartment_final"],
    }).to_csv(
        tmp_path / "stratified_o2_subcompartment_final_within_cruise_group_indicator_species_summary.tsv",
        sep="\t", index=False,
    )

    result = MODULE.load_isa(tmp_path, ["o2_subcompartment_final"])

    assert result.loc[0, "grouping"] == "o2_subcompartment_final"
    assert result.loc[0, "group_level"] == "oxic__gmm0"
    assert result.loc[0, "stratified_within_col"] == "cruise_group"
    assert result.loc[0, "stratified_within_value"] == "cruise_group_1"


def test_load_isa_retains_multigroup_assignment_as_one_combination(tmp_path):
    pd.DataFrame({
        "ASV": ["ASV1"], "s.group_a": [1], "s.group_b": [1],
        "stat": [0.9], "q.value": [0.01],
    }).to_csv(
        tmp_path / "cruise_group_indicator_species_summary.tsv",
        sep="\t", index=False,
    )

    result = MODULE.load_isa(tmp_path, ["cruise_group"])

    assert len(result) == 1
    assert result.loc[0, "group_level"] == "group_a+group_b"
    assert result.loc[0, "association_group_count"] == 2
    assert result.loc[0, "association_scope"] == "multigroup"


def test_matched_compartment_module_tests_uses_common_observed_cohort():
    metadata = pd.DataFrame({
        "sampleID": [f"S{i}" for i in range(9)],
        "Cruise": [f"C{i // 3}" for i in range(9)],
        "o2_compartment": ["oxic"] * 3 + ["dysoxic"] * 3 + ["suboxic"] * 3,
        "gmm_component": [0] * 3 + [1] * 3 + [2] * 3,
        "o2_subcompartment_final": ["oxic__gmm0"] * 3
        + ["dysoxic__gmm1"] * 3 + ["suboxic__gmm2"] * 3,
        "o2_subcompartment_final_assignment_source": ["observed"] * 8 + ["imputed"],
    })
    scores = pd.DataFrame({
        "sampleID": [f"S{i}" for i in range(9)] * 2,
        "ecological_module": ["M1"] * 9 + ["M2"] * 9,
        "module_relative_abundance": list(range(9)) + list(reversed(range(9))),
    })

    tests, summary, cohort = MODULE.matched_compartment_module_tests(
        scores, metadata, "sampleID", "Cruise", permutations=9, seed=42,
        minimum_group_n=2,
    )

    assert set(tests["grouping"]) == {
        "o2_compartment", "gmm_component", "o2_subcompartment_final",
    }
    assert set(tests["common_cohort_samples"]) == {8}
    assert set(summary["classification"]) == {"Legacy O2", "GMM", "Hybrid O2-GMM"}
    assert "adjusted_r_squared" in tests
    assert len(cohort) == 8
    assert set(cohort["o2_subcompartment_final_assignment_source"]) == {"observed"}


def test_seasonal_module_summaries_are_ordered_and_normalized():
    associations = pd.DataFrame({
        "grouping": ["Season", "Season"],
        "ecological_module": ["M10", "M2"],
        "effect_eta_squared": [0.4, 0.2],
        "p_value": [0.01, 0.2],
        "q_value": [0.02, 0.2],
        "highest_group_level": ["Fall", "Winter"],
        "n_units": [20, 20],
        "n_levels": [4, 4],
    })
    levels = pd.DataFrame([
        {
            "grouping": "Season", "ecological_module": module,
            "group_level": season, "sample_unit": "cruise", "n_units": 5,
            "mean_module_relative_abundance": value,
            "median_module_relative_abundance": value,
            "overall_mean_module_relative_abundance": 2.5,
            "mean_difference_from_other_levels": 0.0,
            "mean_ratio_to_other_levels": 1.0,
        }
        for module in ("M10", "M2")
        for season, value in zip(MODULE.SEASON_ORDER, (1.0, 2.0, 3.0, 4.0))
    ])

    association, abundance = MODULE.seasonal_module_summaries(
        associations, levels
    )

    assert association.ecological_module.tolist() == ["M2", "M10"]
    assert abundance.ecological_module.drop_duplicates().tolist() == ["M2", "M10"]
    assert abundance.loc[abundance.ecological_module.eq("M2"), "group_level"].tolist() == MODULE.SEASON_ORDER
    assert abundance.groupby("ecological_module")[
        "within_module_season_fraction"
    ].sum().round(12).eq(1).all()


def test_summarize_compartment_function_tests_reports_corrected_support():
    rows = []
    for grouping in ("o2_compartment", "gmm_component", "o2_subcompartment_final"):
        for level in ("a", "b"):
            rows.append({
                "grouping": grouping,
                "group_level": level,
                "stratified_within_col": "",
                "functional_module_id": "M00529",
                "presence_q_value": 0.01 if grouping != "o2_compartment" else 0.2,
                "completeness_q_value": 0.03 if grouping == "o2_subcompartment_final" else 0.2,
                "function_prevalence": 0.5,
                "odds_ratio": 2.0,
            })
    result = MODULE.summarize_compartment_function_tests(pd.DataFrame(rows))

    hybrid = result.loc[result["classification"].eq("Hybrid O2-GMM")].iloc[0]
    legacy = result.loc[result["classification"].eq("Legacy O2")].iloc[0]
    assert hybrid["presence_enrichments_q05"] == 2
    assert hybrid["completeness_enrichments_q05"] == 2
    assert legacy["presence_enrichments_q05"] == 0


def test_prepare_functional_indicator_evidence_filters_to_significant_single_network_asvs():
    isa = pd.DataFrame([
        {"grouping": "o2_compartment", "group_level": "oxic",
         "association_scope": "single_group", "stratified_within_col": "",
         "stratified_within_value": "", "ASV_ID": "ASV1", "isa_q_value": 0.01},
        {"grouping": "o2_compartment", "group_level": "oxic+suboxic",
         "association_scope": "multigroup", "stratified_within_col": "",
         "stratified_within_value": "", "ASV_ID": "ASV2", "isa_q_value": 0.01},
        {"grouping": "o2_compartment", "group_level": "suboxic",
         "association_scope": "single_group", "stratified_within_col": "",
         "stratified_within_value": "", "ASV_ID": "ASV3", "isa_q_value": 0.20},
        {"grouping": "gmm_component", "group_level": "1",
         "association_scope": "single_group", "stratified_within_col": "",
         "stratified_within_value": "", "ASV_ID": "ASV4", "isa_q_value": 0.01},
    ])
    mapping = pd.DataFrame({"ASV_ID": ["ASV1"], "mag_match_id": ["MAG1"]})

    evidence, audit = MODULE.prepare_functional_indicator_evidence(
        isa, {"ASV1", "ASV2", "ASV3"}, mapping,
        ["o2_compartment", "gmm_component", "o2_subcompartment_final"],
    )

    assert list(evidence["ASV_ID"]) == ["ASV1"]
    o2 = audit.loc[audit["grouping"].eq("o2_compartment")].iloc[0]
    hybrid = audit.loc[
        audit["grouping"].eq("o2_subcompartment_final")
    ].iloc[0]
    assert o2["mag_linked_functional_indicator_asvs"] == 1
    assert bool(o2["functional_test_available"])
    assert hybrid["mag_linked_functional_indicator_asvs"] == 0
    assert not bool(hybrid["functional_test_available"])


def test_asv_pca_positions_collapse_duplicate_cruise_depth_libraries():
    counts = pd.DataFrame(
        {
            "S1": [10.0],
            "S2": [30.0],
            "S3": [20.0],
        },
        index=["ASV1"],
    )
    modules = pd.DataFrame({"Taxon": ["ASV1"], "module_label": ["M1"]})
    metadata = pd.DataFrame({
        "sampleID": ["S1", "S2", "S3"],
        "Cruise": [1, 1, 1],
        "Depth": [10, 10, 20],
    })
    scores = pd.DataFrame({
        "Cruise": [1, 1],
        "Depth_anchored": [10, 20],
        "PC1": [0.0, 2.0],
        "PC2": [0.0, 4.0],
    })

    result = MODULE.build_asv_pca_positions(
        counts, modules, metadata, scores, "sampleID"
    )

    assert len(result) == 1
    # S1/S2 are one location after median aggregation; the two environmental
    # positions then contribute equal relative-abundance weights.
    assert result.loc[0, "matched_environmental_positions_n"] == 2
    assert result.loc[0, "observed_PC1"] == 1.0
    assert result.loc[0, "observed_PC2"] == 2.0


def test_module_phylum_positions_pool_asvs_by_abundance_weight():
    positions = pd.DataFrame({
        "ASV_ID": ["ASV1", "ASV2", "ASV3"],
        "module_label": ["M1", "M1", "M1"],
        "observed_PC1": [0.0, 2.0, 10.0],
        "observed_PC2": [0.0, 4.0, 10.0],
        "relative_abundance_weight_sum": [1.0, 3.0, 1.0],
    })
    taxonomy = pd.DataFrame({
        "Feature ID": ["ASV1;size=10", "ASV2;size=20", "ASV3;size=5"],
        "Taxon": [
            "d__Bacteria; p__Pseudomonadota",
            "d__Bacteria; p__Pseudomonadota",
            "d__Bacteria; p__Bacteroidota",
        ],
    })

    result = MODULE.build_module_phylum_pca_positions(positions, taxonomy)

    pseudo = result.loc[result["phylum"].eq("Pseudomonadota")].iloc[0]
    assert pseudo["asvs_n"] == 2
    assert pseudo["observed_PC1"] == 1.5
    assert pseudo["observed_PC2"] == 3.0


def test_pca_label_selection_uses_full_denominator_and_strict_threshold():
    counts = pd.DataFrame({'s1': [1., 2., 97.], 's2': [1., 2., 97.]},
                          index=['boundary', 'selected', 'unmapped'])
    positions = pd.DataFrame({'ASV_ID': ['boundary', 'selected'],
                              'observed_PC1': [0., 1.], 'observed_PC2': [0., 1.],
                              'module_label': ['M1', 'M2']})
    result = MODULE.select_pca_label_asvs(counts, positions).set_index('ASV_ID')
    assert set(result.index) == {'selected', 'unmapped'}
    assert result.loc['selected', 'mean_relative_abundance'] == .02
    assert result.loc['selected', 'plotted']
    assert not result.loc['unmapped', 'plotted']


def test_local_labels_avoid_points_and_crossing_leaders():
    import numpy as np
    anchors = np.array([[90., 100.], [100., 105.], [110., 100.]])
    half = np.array([[15., 8.]] * 3)
    points = np.vstack([anchors, [[90., 125.], [110., 80.]]])
    radii = np.full(len(points), 5.)
    layout = MODULE.local_label_layout(anchors, half, points, radii,
                                       np.array([0., 0., 220., 220.]))
    for i, (center, end) in enumerate(layout):
        assert np.all(center - half[i] >= 0) and np.all(center + half[i] <= 220)
        gap = np.maximum(np.abs(points - center) - half[i], 0)
        assert np.all((gap ** 2).sum(axis=1) > radii ** 2)
        assert np.linalg.norm(center - anchors[i]) < 100
        for j, (other, other_end) in enumerate(layout[:i]):
            assert not np.all(np.abs(center - other) <= half[i] + half[j])
            assert not MODULE.leader_segments_intersect(anchors[i], end, anchors[j], other_end)


def test_taxonomy_legend_preserves_family_only_assignment():
    selected = pd.DataFrame({'ASV_ID': ['ASV4']})
    taxonomy = pd.DataFrame({'Feature ID': ['ASV4;size=20'],
                             'Taxon': ['d__Bacteria; p__Campilobacterota; f__Arcobacteraceae; g__uncultured']})
    result = MODULE.attach_label_taxonomy(selected, taxonomy)
    assert result.loc[0, 'family'] == 'Arcobacteraceae'
    assert result.loc[0, 'genus'] == 'uncultured'


def test_label_module_peak_uses_matched_hybrid_not_parent_classification():
    selected = pd.DataFrame({'ASV_ID': ['ASV1', 'ASV3'], 'module_label': ['M1', 'M1']})
    stats = pd.DataFrame({'grouping': ['o2_compartment', 'o2_subcompartment_final'],
                          'ecological_module': ['M1', 'M1'],
                          'highest_group_level': ['suboxic', 'suboxic__gmm1'],
                          'highest_group_mean': [.65, .66]})
    result = MODULE.attach_label_module_peaks(selected, stats)
    assert len(result) == 2
    assert result.module_peak_hybrid_compartment.eq('suboxic__gmm1').all()
    assert result.module_peak_mean_abundance.eq(.66).all()
