import importlib.util
from pathlib import Path

import pandas as pd
import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "processes" / "ecological_context_atlas" / "ecological_context_atlas.py"
SPEC = importlib.util.spec_from_file_location("ecological_context_atlas", SCRIPT)
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def test_mag_function_join_uses_normalized_match_id():
    mappings = pd.DataFrame({
        "ASV_ID": ["ASV1"], "genome_id": ["source_A"], "mag_match_id": ["MAG1"],
        "link_pident": [100.0], "mag_tax_species": ["Species one"],
    })
    functions = pd.DataFrame({
        "genome_id": ["different_source_spelling"], "mag_match_id": ["MAG1"],
        "module_id": ["M00001"], "module_name": ["test"], "present": [True],
    })
    network = pd.DataFrame({"ASV_ID": ["ASV1"], "module_label": ["M1"]})
    taxonomy = pd.DataFrame({"ASV_ID": ["ASV1"], "taxonomy": ["d__Bacteria"]})
    result = MOD.mag_function_links(mappings, functions, network, taxonomy)
    assert len(result) == 1
    assert result.loc[0, "module_id"] == "M00001"


def test_context_profiles_retain_explicit_zeros():
    sample = pd.DataFrame({
        "ASV_ID": ["ASV1", "ASV1"], "sample_ID": ["S1", "S2"],
        "relative_abundance": [0.0, 0.2], "hybrid": ["A", "A"],
        "Cruise": [1, 2], "Depth": [10, 20], "Date": ["2020-01-01", "2020-02-01"],
    })
    result = MOD.context_profiles(sample, ["hybrid"], "Cruise", "Depth", "Date")
    row = result.iloc[0]
    assert row.number_of_available_samples == 2
    assert row.number_of_detected_samples == 1
    assert row.prevalence == 0.5


def test_paired_multiomics_requires_exact_keys():
    data = pd.DataFrame({
        "ASV_ID": ["ASV1", "ASV1", "ASV1"], "mag_match_id": ["M1"] * 3,
        "species": ["sp"] * 3, "sample_join_key": ["S1", "S1", "S2"],
        "data_modality": ["metagenome", "metatranscriptome", "metatranscriptome"],
        "normalized_recruitment": [2.0, 4.0, 9.0],
    })
    data["normalization"] = "provided_tpm"
    result = MOD.paired_multiomics(data)
    paired = result.loc[result.sample_join_key.eq("S1")].iloc[0]
    assert paired.paired_dna_rna_available
    assert paired.rna_to_dna_recruitment_ratio == 2.0


def test_paired_multiomics_does_not_divide_tpm_by_fpkm():
    data = pd.DataFrame({
        "ASV_ID": ["ASV1", "ASV1"], "mag_match_id": ["M1", "M1"],
        "species": ["sp", "sp"], "sample_join_key": ["S1", "S1"],
        "data_modality": ["metagenome", "metatranscriptome"],
        "normalized_recruitment": [2.0, 4.0],
        "normalization": ["provided_fpkm", "provided_tpm"],
    })
    result = MOD.paired_multiomics(data)
    assert "rna_to_dna_recruitment_ratio" not in result


def test_paired_multiomics_uses_precomputed_tpm_ratio_and_se(tmp_path):
    data = pd.DataFrame({
        "ASV_ID": ["ASV1", "ASV1"], "mag_match_id": ["MAG1", "MAG1"],
        "species": ["sp", "sp"], "sample_join_key": ["S1", "S1"],
        "data_modality": ["metagenome", "metatranscriptome"],
        "normalized_recruitment": [2.0, 4.0],
        "normalization": ["provided_fpkm", "provided_tpm"],
    })
    ratio = tmp_path / "ratio.tsv"
    se = tmp_path / "se.tsv"
    pd.DataFrame({"genome_id": ["100m__MAG1"], "S1": [1.5]}).to_csv(
        ratio, sep="\t", index=False
    )
    pd.DataFrame({"genome_id": ["100m__MAG1"], "S1": [0.2]}).to_csv(
        se, sep="\t", index=False
    )
    result = MOD.paired_multiomics(
        data, log2_ratio_path=ratio, log2_se_path=se,
        mag_id_mode="suffix_after_double_underscore", max_log2_se=0.3,
    )
    row = result.iloc[0]
    assert row.log2_rna_to_dna_recruitment_ratio == 1.5
    assert row.rna_to_dna_recruitment_ratio == pytest.approx(2 ** 1.5)
    assert row.paired_tpm_ratio_passes_se


def test_multiomics_matching_uses_explicit_sample_field(tmp_path):
    abundance = pd.DataFrame({
        "ASV_ID": ["ASV1"], "species": ["sp"], "mag_match_id": ["M1"],
        "sample": ["SI001_100m"], "normalized_recruitment": [3.0],
    })
    abundance.to_csv(tmp_path / "x_metagenome_ASV1_species_recruitment_by_sample.tsv", sep="\t", index=False)
    sample = pd.DataFrame({
        "ASV_ID": ["ASV1"], "sample_ID": ["SI001_100m_plate8"],
        "sample": ["SI001_100m"], "relative_abundance": [0.2],
    })
    result = MOD.multiomics_context(tmp_path, sample)
    assert result.loc[0, "sample_context_matched"]
    assert result.loc[0, "relative_abundance"] == 0.2


def test_titan_network_encoding_is_grayscale_monotonic_and_bounded():
    assert MOD.TITAN_DIRECTION_COLORS == {"z-": "#BDBDBD", "z+": "#000000"}
    areas = [MOD.titan_marker_area(value) for value in [0, 5, 10, 20, 30, 40]]
    assert areas[:5] == sorted(areas[:5])
    assert areas[-1] == areas[-2]


def test_ecological_linkage_crosswalk_retains_statistics_and_mag_bridge():
    network = pd.DataFrame({
        "ASV_ID": ["ASV1", "ASV2"], "module_label": ["M1", "M1"],
        "node_stability": [0.98, 0.92], "Degree": [10, 5],
    })
    taxonomy = pd.DataFrame({
        "ASV_ID": ["ASV1", "ASV2"], "taxonomy": ["d__Bacteria", "d__Bacteria"],
        "phylum": ["P1", "P2"],
    })
    response = pd.DataFrame({
        "ASV_ID": ["ASV1", "ASV2"], "environmental_variable": ["Oxygen", "Oxygen"],
        "spearman_rho": [0.8, 0.1], "spearman_p_value": [0.001, 0.2],
        "spearman_q_value": [0.01, 0.2], "spearman_n": [20, 20],
        "change_point": [5.0, 7.0], "response_direction": ["z+", "z-"],
        "passes_purity_and_reliability": [True, False],
    })
    compartment = pd.DataFrame({
        "grouping": ["o2_subcompartment_final"], "ecological_module": ["M1"],
        "effect_eta_squared": [0.4], "adjusted_r_squared": [0.35],
        "p_value": [0.001], "q_value": [0.01],
        "highest_group_level": ["oxic__gmm1"], "highest_group_mean": [0.2],
        "n_units": [20], "n_levels": [3],
    })
    measurement = pd.DataFrame({
        "ecological_module": ["M1"], "measurement": ["Oxygen"],
        "rho": [0.7], "p_value": [0.001], "q_value": [0.01],
        "samples_measured": [20], "cruises_measured": [5], "depths_measured": [4],
    })
    member = pd.DataFrame({
        "ecological_module": ["M1"], "measurement": ["Oxygen"],
        "enrichment_odds_ratio": [3.0], "enrichment_p_value": [0.002],
        "enrichment_q_value": [0.02], "significant_asvs": [1],
        "significant_fraction": [0.5], "significant_positive_asvs": [1],
        "significant_negative_asvs": [0], "sign_concordance": [1.0],
        "median_rho_all_asvs": [0.45], "median_rho_significant_asvs": [0.8],
    })
    mappings = pd.DataFrame({
        "ASV_ID": ["ASV1", "ASV1"], "mag_match_id": ["MAG1", "MAG2"],
        "genome_id": ["G1", "G2"], "mag_tax_species": ["sp1", "sp2"],
        "mag_tax_genus": ["g1", "g2"], "mag_tax_phylum": ["p1", "p1"],
        "pairing_status": ["paired_ambiguous", "paired_ambiguous"],
    })
    full, supported, mag_linked, definitions = MOD.ecological_linkage_crosswalk(
        network, taxonomy, response, compartment, measurement, member, mappings, 0.05,
    )
    assert len(full) == 2
    assert supported.ASV_ID.tolist() == ["ASV1"]
    assert mag_linked.ASV_ID.tolist() == ["ASV1"]
    assert mag_linked.iloc[0].analysis_eligible_mag_link_count == 2
    assert mag_linked.iloc[0].analysis_eligible_mag_species == "sp1|sp2"
    assert mag_linked.iloc[0].module_compartment_q_value == 0.01
    assert mag_linked.iloc[0].module_measurement_q_value == 0.01
    assert mag_linked.iloc[0].spearman_q_value == 0.01
    membership = definitions.loc[definitions.linkage.eq("ASV to ecological module")].iloc[0]
    assert membership.inferential_statistic == "None"


def test_asv_ecological_module_mag_crosswalk_is_one_row_per_asv():
    network = pd.DataFrame({
        "ASV_ID": ["ASV1", "ASV2", "ASV3"],
        "module_label": ["M1", "M1", "M2"],
        "module_id": [1, 1, 2],
        "Degree": [10, 5, 2],
    })
    taxonomy = pd.DataFrame({
        "ASV_ID": ["ASV1", "ASV2", "ASV3"],
        "taxonomy": ["d__Bacteria"] * 3,
        "phylum": ["P1", "P2", "P3"],
    })
    mappings = pd.DataFrame({
        "ASV_ID": ["ASV1", "ASV1", "ASV3"],
        "mag_match_id": ["MAG1", "MAG2", "MAG3"],
        "genome_id": ["G1", "G2", "G3"],
        "mag_tax_species": ["sp1", "sp2", "sp3"],
        "mag_tax_genus": ["g1", "g2", "g3"],
        "mag_tax_phylum": ["p1", "p1", "p3"],
        "pairing_status": ["paired_ambiguous", "paired_ambiguous", "paired_unique"],
    })
    result = MOD.asv_ecological_module_mag_crosswalk(
        network, taxonomy, mappings,
    )
    assert result.ASV_ID.tolist() == ["ASV1", "ASV2", "ASV3"]
    assert result.loc[result.ASV_ID.eq("ASV1"), "analysis_eligible_mag_link_count"].item() == 2
    assert result.loc[result.ASV_ID.eq("ASV2"), "analysis_eligible_mag_link_count"].item() == 0
    assert result.loc[result.ASV_ID.eq("ASV1"), "module_mag_linked_asv_count"].item() == 1
    assert result.loc[result.ASV_ID.eq("ASV3"), "module_mag_linked_asv_fraction"].item() == 1.0


def test_ecological_linkage_sources_are_explicitly_wired():
    root = Path(__file__).resolve().parents[1]
    workflow = (root / "workflow" / "modules" / "networks_reporting.nf").read_text()
    pipeline = (root / "asv_pipeline.nf").read_text()
    for argument in [
        "--module-compartment-associations",
        "--module-measurement-associations",
        "--module-measurement-member-support",
        "--linkage-q-threshold",
    ]:
        assert argument in workflow
    assert "group_guild_function_stage.done" in pipeline
