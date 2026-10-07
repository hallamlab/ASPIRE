#!/usr/bin/env python3
"""Validate an ASPIRE mock run against its dataset and publication contract."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd
import yaml


class Audit:
    def __init__(self):
        self.failures: list[str] = []

    def check(self, condition: bool, label: str, detail: str) -> None:
        status = "PASS" if condition else "FAIL"
        print(f"[{status}] {label}: {detail}")
        if not condition:
            self.failures.append(f"{label}: {detail}")


def read_table(path: Path, **kwargs) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t", **kwargs)


def sample_ids_from_manifest(path: Path) -> set[str]:
    frame = pd.read_csv(path, sep="\t", comment="#")
    if "sample_id" in frame.columns:
        return set(frame["sample_id"].dropna().astype(str))
    frame = pd.read_csv(path, sep="\t", comment="#", header=None)
    return set(frame.iloc[:, 0].dropna().astype(str))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def read_fasta(path: Path) -> dict[str, str]:
    opener = gzip.open if path.suffix == ".gz" else open
    records: dict[str, str] = {}
    identifier = None
    sequence: list[str] = []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if identifier is not None:
                    records[identifier] = "".join(sequence).upper()
                identifier = line[1:].split(";", 1)[0].split()[0]
                sequence = []
            elif line:
                sequence.append(line)
    if identifier is not None:
        records[identifier] = "".join(sequence).upper()
    return records


def map_inferred_to_truth(dataset: Path, results: Path) -> dict[str, str]:
    registry = read_table(dataset / "ground_truth_feature_registry.tsv")
    truth_sequences = {
        str(row.ASV_ID): str(row.v4_sequence).upper()
        for row in registry[["ASV_ID", "v4_sequence"]].itertuples(index=False)
    }
    inferred = read_fasta(results / "intermediates/ASVs/ASVs.fasta.gz")
    mapping: dict[str, str] = {}
    for inferred_id, sequence in inferred.items():
        oriented = (sequence, reverse_complement(sequence))
        matches = [
            truth_id for truth_id, truth_sequence in truth_sequences.items()
            if any(candidate in truth_sequence or truth_sequence in candidate for candidate in oriented)
        ]
        if len(matches) == 1:
            mapping[inferred_id] = matches[0]
    return mapping


def true_values(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def validate_control_filtering(metadata: pd.DataFrame, results: Path, audit: Audit, cami: bool) -> None:
    """Validate class separation, inclusion, statistical union and unchanged counts.

    Use the raw post-QC matrix as the independent source of expected counts;
    checks never assume that merely detecting an ASV in a control removes it.
    """
    directory = results / "modules/contamination_filtering/tables/three_tier_results"
    settings = json.loads((directory / "settings.json").read_text())
    qc = read_table(directory / "sample_qc.tsv").set_index(settings["sample_col"], drop=False)
    raw = read_table(results / "intermediates/ASVs/ASV_counts.tsv", index_col=0)
    calls = read_table(directory / "contamination_calls.tsv", index_col=0)
    cleaned = read_table(directory / "ASV_cleaned.tsv", index_col=0)
    normalize = lambda values: pd.Index(values).astype(str).str.replace("-", "_", regex=False)
    source = metadata.copy()
    source.index = normalize(source["sample_id"])
    source = source.reindex(normalize(qc.index))
    source.index = qc.index
    expected_classes = source["Type_Group"].map({"Airways": "biological", "Oral": "biological",
                                                "Skin": "bio_control", "Control": "technical"}) if cami else qc.decontam_class
    audit.check(expected_classes.notna().all() and expected_classes.equals(qc.decontam_class),
                "control_classes", f"class counts={qc.decontam_class.value_counts().to_dict()}")
    if "source_sample" in source:
        source["source_sample"] = source["source_sample"].fillna(source["sample_id"])
    pairing = [c for c in metadata if c in {"Participant_ID", "Procedure_ID", "source_sample"}]
    pairing_ok = all(c in qc and source[c].fillna("").astype(str).eq(qc[c].fillna("").astype(str)).all()
                     for c in pairing)
    audit.check(set(normalize(qc.index)) == set(normalize(metadata.sample_id)) and pairing_ok,
                "control_metadata", f"QC samples={len(qc)}, pairing columns={pairing}")
    depths = raw.sum().reindex(qc.index, fill_value=0)
    eligible = (expected_classes == "biological") & (depths >= settings["min_biological_reads"]) & (depths > 0)
    bio_ids = qc.index[eligible].tolist()
    audit.check(depths.eq(qc.post_qc_reads).all() and true_values(qc.biological_pass).eq(eligible).all(),
                "biological_depth_qc", f"cutoff={settings['min_biological_reads']}, passing={len(bio_ids)}")

    def same_counts(actual, expected):
        return (set(actual.index) == set(expected.index) and set(actual.columns) == set(expected.columns)
                and actual.reindex(index=expected.index, columns=expected.columns).eq(expected).all().all())

    audit.check(same_counts(read_table(directory / "biological_qc_counts.tsv", index_col=0), raw.loc[:, bio_ids]),
                "biological_qc_counts", "QC biological counts equal raw counts")
    enabled_calls = {}
    if cami:
        audit.check(settings["technical"]["enabled"] and settings["bio_control"]["enabled"],
                    "cami_control_arms", "TECH and BIO must both be enabled")
    for arm, cls in (("TECH", "technical"), ("BIO", "bio_control")):
        negatives = qc.index[(expected_classes == cls) & (depths > 0)].tolist()
        if settings[cls]["enabled"]:
            arm_counts = read_table(directory / f"{arm}_counts.tsv", index_col=0)
            arm_meta = read_table(directory / f"{arm}_metadata.tsv").set_index("Sample")
            selected = bio_ids + negatives
            audit.check(same_counts(arm_counts, raw.loc[:, selected]) and
                        set(arm_meta.index[true_values(arm_meta.is_negative)]) == set(negatives),
                        f"{arm}_cohort", f"biological={len(bio_ids)}, controls={len(negatives)}, low-depth controls={sum(depths.loc[negatives] < settings['min_biological_reads'])}")
            scores = read_table(directory / f"{arm}_scores.tsv", index_col=0).reindex(calls.index)
            expected = scores.score.lt(settings[cls]["threshold"]).fillna(False)
            enabled_calls[arm] = expected
            audit.check(expected.eq(true_values(scores.contaminant)).all() and
                        expected.eq(true_values(calls[f"{arm}_contaminant"])).all() and
                        calls[f"{arm}_score"].fillna(-1).eq(scores.score.fillna(-1)).all(),
                        f"{arm}_calls", f"threshold={settings[cls]['threshold']}, flagged={int(expected.sum())}")
        else:
            enabled_calls[arm] = pd.Series(False, index=calls.index)
        for klass, ids in ((cls, negatives), ("biological", bio_ids)):
            prevalence = (raw.loc[:, ids] > 0).mean(axis=1).reindex(calls.index) if ids else pd.Series(float('nan'), index=calls.index)
            audit.check((calls[f"prevalence_{klass}"].fillna(-1) - prevalence.fillna(-1)).abs().lt(1e-12).all()
                        and calls[f"n_{klass}"].eq(len(ids)).all(),
                        f"{arm}_{klass}_prevalence", f"denominator={len(ids)}")
    tech, bio = enabled_calls["TECH"], enabled_calls["BIO"]
    category = pd.Series("CLEAN", index=calls.index)
    category.loc[tech] = "TECH"
    category.loc[bio] = "BIO"
    category.loc[tech & bio] = "TECH+BIO"
    removed = set(calls.index[tech | bio])
    audit.check(set(calls.index) == set(raw.index) and category.eq(calls.final_category).all()
                and same_counts(cleaned, raw.loc[~raw.index.isin(removed), bio_ids]),
                "contamination_union", f"categories={category.value_counts().to_dict()}; retained counts unchanged")
    reasons = read_table(directory / "removed_asvs_with_taxonomy.tsv", index_col=0)
    audit.check(set(reasons.index) == removed and reasons.removal_reason.eq(category.reindex(reasons.index)).all(),
                "contamination_reasons", f"removed ASVs={len(reasons)}")
    final = read_table(results / "modules/non_target_filtering/tables/ASV_target.tsv", index_col=0)
    final_sequences = read_fasta(results / "intermediates/ASVs/ASVs_target.fasta.gz")
    audit.check(set(final.index) <= set(cleaned.index) and set(final.columns) <= set(bio_ids)
                and same_counts(final, cleaned.reindex(index=final.index, columns=final.columns))
                and set(final_sequences) == set(final.index),
                "combined_filter_counts_fasta", f"ASVs={len(final)}, biological samples={len(final.columns)}")


    micro = read_table(results / "modules/non_target_filtering/tables/ASV_target.micro.tsv", index_col=0)
    feature_qc = read_table(results / "modules/non_target_filtering/tables/ASV_target.feature_qc.tsv", index_col=0)
    config = yaml.safe_load((results / "summary/tables/run_config.yml").read_text())
    filtering = config.get("standard", {}).get("filter_counts", config.get("filter_counts", {}))
    min_ra = filtering.get("min_relative_abundance_pct", filtering.get("abundance_threshold", 0.005))
    min_prev = filtering.get("min_prevalence_fraction", 0)
    presence = (micro > 0).sum(axis=1)
    prevalence = presence / len(micro.columns)
    ra = micro.div(micro.sum().replace(0, float('nan')), axis=1).fillna(0).mul(100).max(axis=1)
    expected_pass = prevalence.ge(min_prev) & ra.ge(min_ra) & presence.gt(0)
    q = feature_qc.reindex(micro.index)
    audit.check(set(micro.columns) == set(final.columns) and set(micro.columns) <= set(bio_ids)
                and set(micro.index) <= set(cleaned.index)
                and same_counts(micro, cleaned.reindex(index=micro.index, columns=micro.columns))
                and set(feature_qc.index) == set(micro.index)
                and q.n_biological_samples.eq(len(micro.columns)).all()
                and q.n_nonzero_samples.eq(presence).all()
                and (q.prevalence_fraction - prevalence).abs().lt(1e-12).all()
                and (q.max_relative_abundance_pct - ra).abs().lt(1e-10).all()
                and q.min_prevalence_fraction.eq(min_prev).all()
                and q.min_relative_abundance_pct.eq(min_ra).all()
                and true_values(q.passes_prevalence).eq(prevalence.ge(min_prev)).all()
                and true_values(q.passes_relative_abundance).eq(ra.ge(min_ra)).all()
                and true_values(q.passes_feature_filters).eq(expected_pass).all()
                and true_values(q.retained_final).eq(pd.Series(micro.index.isin(final.index), index=micro.index)).all()
                and set(final.index) <= set(micro.index[expected_pass]),
                "biological_prevalence_abundance", f"RA>={min_ra}%; nonzero prevalence>={min_prev}; biological denominator={len(micro.columns)}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--results", required=True, type=Path)
    args = parser.parse_args()
    dataset = args.dataset.expanduser().resolve()
    results = args.results.expanduser().resolve()
    if not dataset.is_dir():
        raise SystemExit(f"Mock dataset directory not found: {dataset}")
    if not results.is_dir():
        raise SystemExit(f"Completed ASPIRE output directory not found: {results}")
    modules = results / "modules"
    audit = Audit()

    required_files = (
        "summary/tables/run_manifest.tsv",
        "summary/tables/module_output_manifest.tsv",
        "summary/tables/nextflow_artifact_manifest.tsv",
        "intermediates/ASVs/ASVs.fasta.gz",
        "intermediates/ASVs/ASV_counts.tsv",
        "intermediates/ASVs/ASVs_target.fasta.gz",
        "intermediates/ASVs/filter_audit.tsv",
        "modules/non_target_filtering/tables/ASV_target.tsv",
        "modules/non_target_filtering/tables/ASV_target.micro.tsv",
        "modules/non_target_filtering/tables/ASV_target.feature_qc.tsv",
        "summary/tables/run_config.yml",
        "modules/non_target_filtering/tables/mitomap/nontarget.master.tsv",
        "modules/taxonomy/tables/ASV_SILVA_tax.full-length.vsearch.tsv",
        "modules/voc_correlation/tables/selected_sample_metadata.tsv",
        "modules/voc_correlation/tables/sample_voc_matrix.tsv",
        "modules/voc_correlation/tables/patient_voc_matrix.tsv",
        "modules/voc_correlation/plots/asv_voc_clustermap.svg",
        "modules/voc_correlation/tables/asv_voc_spearman_long.tsv",
        "modules/voc_correlation/tables/patient_asv_voc_permutation_long.tsv",
        "modules/voc_correlation/tables/patient_inference_summary.json",
        "modules/network_analysis/tables/spieceasi_edge_list.csv",
        "modules/network_analysis/tables/spieceasi_modules_all.tsv",
        "modules/network_analysis/tables/network_topology_summary.tsv",
        "modules/network_analysis/tables/network_topology_null_draws.tsv",
        "modules/contamination_filtering/tables/ASV_final_three_tier.tsv",
        "modules/contamination_filtering/tables/three_tier_results/settings.json",
        "modules/contamination_filtering/tables/three_tier_results/sample_qc.tsv",
        "modules/contamination_filtering/tables/three_tier_results/contamination_calls.tsv",
        "modules/contamination_filtering/tables/three_tier_results/removed_asvs_with_taxonomy.tsv",
        "modules/contamination_filtering/plots/three_tier_results/read_depth_by_class.svg",
        "modules/contamination_filtering/tables/three_tier_results/filtered/filter_summary.txt",
        "modules/non_target_filtering/plots/mitomap/nontarget_non_target_cumulative.svg",
        "modules/outlier_detection/plots/outliers_Case_summary.svg",
        "modules/umap_clustering/plots/umap_clustering_type_group.svg",
        "summary/plots/module_output_summary.svg",
        "summary/report/ASPIRE_run_report.html",
        "logs/controller.log",
        "logs/task_execution.tsv",
        "logs/nextflow_report.html",
        "logs/nextflow_timeline.html",
        "logs/nextflow_trace.tsv",
        "logs/nextflow_dag.html",
        "logs/launch_command.txt",
        "logs/nextflow_version.txt",
    )
    metadata = read_table(dataset / "sample_metadata.tsv")
    cami = {'source_sample', 'body_site', 'study_role'} <= set(metadata.columns) and {
        'Airways','Oral','Skin'} <= set(metadata.body_site.dropna())
    clinical = 'synthetic_clinical' in metadata and metadata['synthetic_clinical'].fillna(False).astype(str).str.lower().isin(['true','1']).any()
    if cami and not clinical:
        required_files = tuple(p.replace('outliers_Case_summary', 'outliers_batch_summary')
                               for p in required_files if '/patient_' not in p)
    missing_files = [relative for relative in required_files if not (results / relative).is_file()]
    if missing_files:
        raise SystemExit("Completed output is missing required files:\n  " + "\n  ".join(missing_files))

    required_modules = {
        "batch_correction", "bubbleplotter", "clustermaps", "collectors_curve", "contamination_filtering",
        "diversity", "general_stats", "indicator_analysis", "lung_status_analysis",
        "metadata_plots", "network_analysis", "non_target_filtering",
        "outlier_detection", "power_analysis", "sankey", "taxonomy",
        "umap_clustering", "upset", "voc_correlation",
    }
    if cami and not clinical:
        required_modules -= {'lung_status_analysis', 'power_analysis'}
    missing_modules = sorted(name for name in required_modules if not (modules / name).is_dir())
    audit.check(not missing_modules, "publication_layout", f"missing modules={missing_modules or 'none'}")

    supplied = sample_ids_from_manifest(dataset / "fastq_manifest.tsv")
    metadata = read_table(dataset / "sample_metadata.tsv")
    published = sample_ids_from_manifest(results / "summary/tables/run_manifest.tsv")
    metadata_ids = set(metadata["sample_id"].astype(str))
    audit.check(
        supplied == metadata_ids == published,
        "sample_accounting",
        f"manifest={len(supplied)}, metadata={len(metadata_ids)}, published={len(published)}",
    )
    validate_control_filtering(metadata, results, audit, cami)

    truth = read_table(dataset / "ground_truth_reference_filters.tsv")
    nontarget = read_table(modules / "non_target_filtering/tables/mitomap/nontarget.master.tsv")
    final = read_table(modules / "non_target_filtering/tables/ASV_target.tsv")
    final_ids = set(final.iloc[:, 0].astype(str))
    identified: set[str] = set()
    for column in ("BF_ID", "MI_Accession"):
        if column in nontarget:
            identified.update(nontarget[column].dropna().astype(str))
    truth_ids = set(truth["ASV_ID"].astype(str))
    identified_rows = nontarget[nontarget.apply(
        lambda row: any(str(row.get(column, "")) in truth_ids for column in ("BF_ID", "MI_Accession")), axis=1
    )]
    removed = set(identified_rows["Sequence_ID"].astype(str)).isdisjoint(final_ids)
    audit.check(
        truth_ids <= identified and removed,
        "non_target_truth",
        f"identified={len(truth_ids & identified)}/{len(truth_ids)}, excluded_from_micro={removed}",
    )

    taxonomy = read_table(modules / "taxonomy/tables/ASV_SILVA_tax.full-length.vsearch.tsv")
    tax_col = "Taxon" if "Taxon" in taxonomy else taxonomy.columns[1]
    n_tax = int(taxonomy[tax_col].fillna("").astype(str).str.strip().ne("").sum())
    audit.check(n_tax > 0, "taxonomy", f"assigned ASVs={n_tax}")

    id_map = map_inferred_to_truth(dataset, results)
    audit.check(bool(id_map), "truth_id_mapping", f"mapped inferred ASVs={len(id_map)}")

    indicator_dir = modules / "indicator_analysis/tables"
    indicator_files = [indicator_dir / "Type_Group_indicator_species_summary.tsv", indicator_dir / ("batch_indicator_species_summary.tsv" if cami and not clinical else "Case_indicator_species_summary.tsv")]
    significant = 0
    significant_ids: set[str] = set()
    for path in indicator_files:
        if path.is_file():
            frame = read_table(path)
            if "significant" in frame:
                frame = frame[true_values(frame["significant"])]
                significant += len(frame)
                significant_ids.update(frame["ASV"].astype(str))
    group_truth = read_table(dataset / "ground_truth_group_effects.tsv")
    significant_truth = {id_map[item] for item in significant_ids if item in id_map}
    recovered_indicators = significant_truth & set(group_truth["ASV_ID"].astype(str))
    audit.check(
        bool(recovered_indicators),
        "indicator_analysis",
        f"significant={significant}, implanted_recovered={len(recovered_indicators)}",
    )

    voc = read_table(modules / "voc_correlation/tables/asv_voc_spearman_long.tsv")
    positive = voc[voc["rho"] >= 0.30].copy() if "rho" in voc else voc.iloc[0:0].copy()
    positive["truth_asv"] = positive["ASV_ID"].astype(str).map(id_map)
    chemistry_truth = read_table(dataset / "ground_truth_asv_chem.tsv")
    chemistry_truth = chemistry_truth[chemistry_truth["direction"].astype(str).str.lower() == "positive"]
    truth_pairs = set(zip(chemistry_truth["ASV_ID"].astype(str), chemistry_truth["compound"].astype(str)))
    recovered_voc = {
        (str(row.truth_asv), str(row.voc)) for row in positive.itertuples(index=False)
        if (str(row.truth_asv), str(row.voc)) in truth_pairs
    }
    audit.check(
        bool(recovered_voc),
        "voc_correlation",
        f"rho>=0.30 pairs={len(positive)}, implanted_positive_recovered={len(recovered_voc)}",
    )

    edge_path = modules / "network_analysis/tables/spieceasi_edge_list.csv"
    if not cami or clinical:
        patient_voc = read_table(modules / "voc_correlation/tables/patient_asv_voc_permutation_long.tsv")
        patient_summary = json.loads((modules / "voc_correlation/tables/patient_inference_summary.json").read_text())
        tested = patient_voc.loc[patient_voc.status.eq("tested")]
        untested = patient_voc.loc[patient_voc.status.ne("tested")]
        audit.check(
            set(patient_voc.normalization) == {"relative_abundance", "clr"}
            and patient_voc.n_patients.le(patient_summary["patients"]).all()
            and tested.p_value.between(0, 1, inclusive="right").all()
            and tested.q_value.between(0, 1).all()
            and tested.rho.abs().le(1 + 1e-12).all()
            and untested.p_value.isna().all(),
            "voc_patient_inference",
            f"patients={patient_summary['patients']}, tested={len(tested)}, insufficient={len(untested)}; significance is not required",
        )

    edges = pd.read_csv(edge_path) if edge_path.is_file() else pd.DataFrame()
    module_path = modules / "network_analysis/tables/spieceasi_modules_all.tsv"
    assignments = read_table(module_path) if module_path.is_file() else pd.DataFrame()
    max_module = int(assignments.groupby("module_label").size().max()) if not assignments.empty else 0
    network_truth = read_table(dataset / "ground_truth_network_modules.tsv")
    truth_module = dict(zip(network_truth["ASV_ID"].astype(str), network_truth["module"].astype(str)))
    within_truth = 0
    for row in edges.itertuples(index=False):
        left = id_map.get(str(row.Taxon1))
        right = id_map.get(str(row.Taxon2))
        if left in truth_module and right in truth_module and truth_module[left] == truth_module[right]:
            within_truth += 1
    audit.check(
        len(edges) >= 10 and max_module >= 2 and within_truth >= 3,
        "network_analysis",
        f"edges={len(edges)}, largest_module={max_module}, implanted_within_module_edges={within_truth}",
    )

    topology = read_table(
        modules / "network_analysis/tables/network_topology_summary.tsv"
    ).set_index("metric")["value"]
    null_draws = read_table(
        modules / "network_analysis/tables/network_topology_null_draws.tsv"
    )
    topology_nodes = int(float(topology["n_nodes"]))
    topology_edges = int(float(topology["n_edges"]))
    null_completed = int(float(topology["n_null_completed"]))
    audit.check(
        topology_nodes > 0
        and topology_edges > 0
        and null_completed == len(null_draws)
        and null_completed > 0,
        "network_topology",
        f"nodes={topology_nodes}, thresholded_positive_edges={topology_edges}, null_draws={null_completed}",
    )

    plot_modules = (
        "clustermaps", "diversity", "indicator_analysis", "network_analysis",
        "non_target_filtering", "outlier_detection", "sankey", "taxonomy", "voc_correlation",
    )
    missing_svgs = [name for name in plot_modules if not any((modules / name / "plots").rglob("*.svg"))]
    audit.check(not missing_svgs, "plot_contract", f"modules without SVG={missing_svgs or 'none'}")

    manifest_path = results / "summary/tables/module_output_manifest.tsv"
    manifest = read_table(manifest_path)
    bad_manifest = []
    for row in manifest.itertuples(index=False):
        path = results / str(row.relative_path)
        if not path.is_file() or sha256(path) != str(row.sha256):
            bad_manifest.append(str(row.relative_path))
    audit.check(not bad_manifest, "manifest_integrity", f"invalid entries={len(bad_manifest)}")

    print(f"\nValidated dataset with {len(supplied)} samples against {results}")
    if audit.failures:
        print("Validation failed:", file=sys.stderr)
        for failure in audit.failures:
            print(f"  - {failure}", file=sys.stderr)
        raise SystemExit(1)
    print("All mock-run checks passed.")


if __name__ == "__main__":
    main()
