#!/usr/bin/env python3
"""Infer auditable genome-species proportionality co-occurrence networks."""

from __future__ import annotations

import argparse
import json
import math
import platform
import re
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t" if path.suffix.lower() in {".tsv", ".tab"} else ",", low_memory=False)


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def multiplicative_replacement(counts: np.ndarray, fraction: float) -> tuple[np.ndarray, pd.DataFrame]:
    """Close samples and replace zeros by a fraction of their smallest positive part."""
    if not 0 < fraction < 1:
        raise ValueError("zero-replacement fraction must be between zero and one")
    x = np.asarray(counts, dtype=float)
    if x.ndim != 2 or np.any(~np.isfinite(x)) or np.any(x < 0):
        raise ValueError("counts must be a finite, nonnegative feature-by-sample matrix")
    out = np.zeros_like(x)
    audits = []
    for column in range(x.shape[1]):
        values = x[:, column]
        total = float(values.sum())
        if total <= 0:
            raise ValueError(f"sample column {column} has zero total abundance")
        composition = values / total
        zero = composition <= 0
        zero_n = int(zero.sum())
        if zero_n == 0:
            out[:, column] = composition
            delta = 0.0
        else:
            positive = composition[~zero]
            if positive.size == 0:
                raise ValueError(f"sample column {column} contains no positive feature")
            delta = float(fraction * positive.min())
            # Guarantee positive remaining mass even for exceptionally sparse samples.
            delta = min(delta, 0.5 / zero_n)
            remaining = 1.0 - zero_n * delta
            out[zero, column] = delta
            out[~zero, column] = composition[~zero] * remaining / positive.sum()
        audits.append({
            "sample_column_index": column,
            "library_total": total,
            "zero_features_n": zero_n,
            "replacement_part": delta,
            "closed_sum": float(out[:, column].sum()),
        })
    return out, pd.DataFrame(audits)


def clr(composition: np.ndarray) -> np.ndarray:
    logged = np.log(composition)
    return logged - logged.mean(axis=0, keepdims=True)


def rho_p_matrix(clr_matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Lovell/Erb proportionality rho and pairwise log-ratio variance."""
    values = np.asarray(clr_matrix, dtype=float)
    variances = np.var(values, axis=1, ddof=1)
    covariance = np.cov(values, ddof=1)
    vlr = variances[:, None] + variances[None, :] - 2.0 * covariance
    denominator = variances[:, None] + variances[None, :]
    rho = np.divide(denominator - vlr, denominator, out=np.zeros_like(vlr), where=denominator > 0)
    rho = np.clip((rho + rho.T) / 2.0, -1.0, 1.0)
    vlr = np.maximum((vlr + vlr.T) / 2.0, 0.0)
    np.fill_diagonal(rho, 1.0)
    np.fill_diagonal(vlr, 0.0)
    return rho, vlr


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.clip(adjusted, 0.0, 1.0)
    return result


def collapse_metadata(metadata: pd.DataFrame, sample_col: str, required: list[str]) -> pd.DataFrame:
    missing = {sample_col, *required} - set(metadata.columns)
    if missing:
        raise ValueError(f"metadata lacks columns: {sorted(missing)}")
    rows = []
    for sample, frame in metadata.groupby(sample_col, dropna=False, sort=False):
        record = {sample_col: str(sample)}
        for column in required:
            observed = frame[column].dropna().astype(str).str.strip()
            observed = observed.loc[observed.ne("")].unique()
            if len(observed) > 1:
                raise ValueError(f"metadata has conflicting {column} values for sample {sample}: {observed.tolist()}")
            record[column] = observed[0] if len(observed) else np.nan
        rows.append(record)
    return pd.DataFrame(rows)


def canonical_number(value: object) -> str:
    number = float(value)
    return str(int(number)) if number.is_integer() else f"{number:g}"


def season_from_month(month: object) -> str:
    value = int(float(month))
    if value in {12, 1, 2}:
        return "Winter"
    if value in {3, 4, 5}:
        return "Spring"
    if value in {6, 7, 8}:
        return "Summer"
    if value in {9, 10, 11}:
        return "Fall"
    raise ValueError(f"invalid calendar month: {month}")


def align_metadata(
    raw: pd.DataFrame,
    sample_ids: list[str],
    *,
    sample_col: str,
    cruise_col: str,
    season_col: str,
    depth_col: str,
    month_col: str,
    sample_id_regex: str,
    required: list[str],
) -> pd.DataFrame:
    """Align direct sample metadata or construct it from an explicit SI sample regex."""
    if sample_col in raw.columns:
        return collapse_metadata(raw, sample_col, required).set_index(sample_col).reindex(sample_ids)
    if cruise_col not in raw.columns or month_col not in raw.columns:
        raise ValueError(
            f"metadata lacks {sample_col}; sample parsing therefore requires {cruise_col} and {month_col}"
        )
    pattern = re.compile(sample_id_regex)
    parsed = []
    for sample in sample_ids:
        match = pattern.fullmatch(str(sample))
        if not match or not {"cruise", "depth"}.issubset(match.groupdict()):
            raise ValueError(f"sample ID does not match configured regex with cruise/depth groups: {sample}")
        parsed.append({
            sample_col: str(sample),
            cruise_col: canonical_number(match.group("cruise")),
            depth_col: canonical_number(match.group("depth")),
        })
    result = pd.DataFrame(parsed)
    cruise_reference = raw.copy()
    cruise_reference[cruise_col] = cruise_reference[cruise_col].map(canonical_number)
    if season_col not in cruise_reference.columns:
        cruise_reference[season_col] = cruise_reference[month_col].map(season_from_month)
    reference_columns = list(dict.fromkeys([cruise_col, season_col, *[x for x in required if x != depth_col]]))
    reference_rows = []
    for cruise, frame in cruise_reference.groupby(cruise_col, sort=False):
        record = {cruise_col: cruise}
        for column in reference_columns:
            if column == cruise_col:
                continue
            values = frame[column].dropna().astype(str).str.strip().unique()
            if len(values) > 1:
                raise ValueError(f"cruise metadata has conflicting {column} values for cruise {cruise}")
            record[column] = values[0] if len(values) else np.nan
        reference_rows.append(record)
    result = result.merge(pd.DataFrame(reference_rows), on=cruise_col, how="left", validate="many_to_one")
    return result.set_index(sample_col).reindex(sample_ids)


def bootstrap_indices(metadata: pd.DataFrame, cruise_col: str, season_col: str, rng: np.random.Generator) -> np.ndarray:
    selected: list[int] = []
    for _, season_frame in metadata.groupby(season_col, sort=True):
        cruises = season_frame[cruise_col].dropna().astype(str).unique()
        if len(cruises) == 0:
            continue
        sampled = rng.choice(cruises, size=len(cruises), replace=True)
        for cruise in sampled:
            selected.extend(season_frame.index[season_frame[cruise_col].astype(str).eq(cruise)].tolist())
    return np.asarray(selected, dtype=int)


def permute_within_strata(values: np.ndarray, strata: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    permuted = np.asarray(values, dtype=float).copy()
    for feature in range(permuted.shape[0]):
        for stratum in np.unique(strata):
            indices = np.flatnonzero(strata == stratum)
            if len(indices) > 1:
                permuted[feature, indices] = permuted[feature, rng.permutation(indices)]
    return permuted


def network_outputs(features: list[str], selected: pd.DataFrame, outdir: Path, prefix: str, seed: int) -> None:
    graph = nx.Graph()
    for feature in features:
        graph.add_node(feature, name=feature, Taxon=feature)
    for row in selected.itertuples(index=False):
        graph.add_edge(
            row.feature_1, row.feature_2,
            weight=float(row.rho_p),
            abs_weight=abs(float(row.rho_p)),
            direction=str(row.direction),
            bootstrap_recovery=float(row.bootstrap_recovery),
            sign_consistency=float(row.sign_consistency),
            q_value=float(row.q_value),
        )

    if graph.number_of_edges():
        distance_graph = graph.copy()
        for u, v, data in distance_graph.edges(data=True):
            data["distance"] = 1.0 / max(abs(float(data["weight"])), 1e-12)
        degree = dict(graph.degree())
        betweenness = nx.betweenness_centrality(distance_graph, weight="distance", normalized=True)
        closeness = nx.closeness_centrality(distance_graph, distance="distance")
        # Eigenvector centrality is not uniquely defined across disconnected
        # components. Calculate it within each connected component so retained
        # isolates and small components remain explicit and deterministic.
        eigenvector = {feature: 0.0 for feature in features}
        for members in nx.connected_components(graph):
            component = graph.subgraph(members)
            if component.number_of_nodes() < 2 or component.number_of_edges() == 0:
                continue
            component_values = nx.eigenvector_centrality(
                component, weight="abs_weight", max_iter=2000, tol=1e-10
            )
            eigenvector.update({node: float(value) for node, value in component_values.items()})
        communities = list(nx.community.louvain_communities(graph, weight="abs_weight", seed=seed))
    else:
        degree = {feature: 0 for feature in features}
        betweenness = {feature: 0.0 for feature in features}
        closeness = {feature: 0.0 for feature in features}
        eigenvector = {feature: 0.0 for feature in features}
        communities = [{feature} for feature in features]
    cluster_lookup = {}
    for number, members in enumerate(sorted(communities, key=lambda group: min(features.index(x) for x in group)), start=1):
        for member in members:
            cluster_lookup[member] = number

    cluster_prefix = {
        "genome_metagenome": "MG-PC",
        "genome_metatranscriptome": "MT-PC",
    }.get(prefix, "GN-PC")

    nodes = pd.DataFrame({
        "GraphML_ID": [f"n{i}" for i in range(len(features))],
        "Taxon": features,
        "Degree": [degree[x] for x in features],
        "Degree_thresholded": [degree[x] for x in features],
        "Betweenness": [betweenness[x] for x in features],
        "Betweenness_raw": [betweenness[x] for x in features],
        "Betweenness_norm": [betweenness[x] for x in features],
        "Closeness": [closeness[x] for x in features],
        "EigenCentral": [eigenvector[x] for x in features],
    })
    nodes.to_csv(outdir / f"{prefix}_node_features.csv", index=False)
    clusters = pd.DataFrame({
        "Taxon": features,
        "proportionality_cluster_id": [cluster_lookup[x] for x in features],
        "proportionality_cluster": [f"{cluster_prefix}{cluster_lookup[x]}" for x in features],
        "node_stability": np.nan,
        "graph_variant": "selected_signed_proportionality",
        "method": "louvain_absolute_rho_p" if graph.number_of_edges() else "singleton",
    })
    clusters.to_csv(outdir / f"{prefix}_proportionality_clusters.tsv", sep="\t", index=False)
    nx.write_graphml(graph, outdir / f"{prefix}_network_all.graphml")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--metadata-sample-col", required=True)
    parser.add_argument("--cruise-col", required=True)
    parser.add_argument("--season-col", required=True)
    parser.add_argument("--depth-col", required=True)
    parser.add_argument("--month-col", default="Month")
    parser.add_argument("--sample-id-regex", default=r"^SI(?P<cruise>[0-9]+)_(?P<depth>[0-9]+(?:[.][0-9]+)?)m$")
    parser.add_argument("--zero-replacement-fraction", type=float, default=0.65)
    parser.add_argument("--min-relative-abundance", type=float, default=0.0)
    parser.add_argument("--min-prevalence", type=float, default=0.05)
    parser.add_argument("--min-abs-rho", type=float, default=0.30)
    parser.add_argument("--bootstrap-iterations", type=int, default=1000)
    parser.add_argument("--min-bootstrap-recovery", type=float, default=0.80)
    parser.add_argument("--min-sign-consistency", type=float, default=0.90)
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--permutation-strata", default="Season,Depth")
    parser.add_argument("--max-q", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    if args.bootstrap_iterations < 1 or args.permutations < 1:
        raise ValueError("bootstrap iterations and permutations must both be positive")
    for label, value in (
        ("minimum relative abundance", args.min_relative_abundance),
        ("minimum prevalence", args.min_prevalence),
        ("minimum bootstrap recovery", args.min_bootstrap_recovery),
        ("minimum sign consistency", args.min_sign_consistency),
        ("maximum q value", args.max_q),
    ):
        if not 0 <= value <= 1:
            raise ValueError(f"{label} must be between zero and one")
    if not 0 <= args.min_abs_rho <= 1:
        raise ValueError("minimum absolute rho-p must be between zero and one")

    matrix = read_table(args.matrix).set_index("Feature_ID")
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    if matrix.isna().any().any() or (matrix < 0).any().any():
        raise ValueError("species matrix contains missing or negative values")
    relative = matrix.divide(matrix.sum(axis=0), axis=1)
    keep = relative.max(axis=1).ge(args.min_relative_abundance) & matrix.gt(0).mean(axis=1).ge(args.min_prevalence)
    filter_audit = pd.DataFrame({
        "feature_id": matrix.index,
        "prevalence": matrix.gt(0).mean(axis=1).to_numpy(),
        "maximum_relative_abundance": relative.max(axis=1).to_numpy(),
        "passes_filter": keep.to_numpy(),
    })
    filter_audit.to_csv(args.outdir / "feature_filter_audit.tsv", sep="\t", index=False)
    matrix = matrix.loc[keep]
    if matrix.shape[0] < 3:
        raise ValueError("fewer than three species passed co-occurrence filtering")

    strata_columns = [x.strip() for x in args.permutation_strata.split(",") if x.strip()]
    required_metadata = list(dict.fromkeys([args.cruise_col, args.season_col, args.depth_col, *strata_columns]))
    metadata = align_metadata(
        read_table(args.metadata), matrix.columns.astype(str).tolist(),
        sample_col=args.metadata_sample_col,
        cruise_col=args.cruise_col,
        season_col=args.season_col,
        depth_col=args.depth_col,
        month_col=args.month_col,
        sample_id_regex=args.sample_id_regex,
        required=required_metadata,
    )
    sample_audit = metadata.reset_index().rename(columns={"index": args.metadata_sample_col})
    sample_audit["matched_metadata"] = metadata[required_metadata].notna().all(axis=1).to_numpy()
    sample_audit.to_csv(args.outdir / "sample_metadata_audit.tsv", sep="\t", index=False)
    if not sample_audit["matched_metadata"].all():
        missing = sample_audit.loc[~sample_audit["matched_metadata"], args.metadata_sample_col].tolist()
        raise ValueError(f"missing complete resampling metadata for {len(missing)} assay samples: {missing[:10]}")
    metadata = metadata.reset_index(drop=True)

    composition, zero_audit = multiplicative_replacement(matrix.to_numpy(float), args.zero_replacement_fraction)
    zero_audit.insert(0, "sample_id", matrix.columns)
    zero_audit.to_csv(args.outdir / "zero_replacement_audit.tsv", sep="\t", index=False)
    clr_values = clr(composition)
    pd.DataFrame(clr_values, index=matrix.index, columns=matrix.columns).rename_axis("Feature_ID").reset_index().to_csv(
        args.outdir / "species_clr_matrix.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    rho, vlr = rho_p_matrix(clr_values)
    iu = np.triu_indices(matrix.shape[0], 1)
    observed = rho[iu]

    rng = np.random.default_rng(args.seed)
    bootstrap = np.empty((args.bootstrap_iterations, len(observed)), dtype=float)
    for iteration in range(args.bootstrap_iterations):
        indices = bootstrap_indices(metadata, args.cruise_col, args.season_col, rng)
        bootstrap[iteration] = rho_p_matrix(clr_values[:, indices])[0][iu]
    recovery = np.mean(np.abs(bootstrap) >= args.min_abs_rho, axis=0)
    sign_consistency = np.mean(np.sign(bootstrap) == np.sign(observed)[None, :], axis=0)

    strata = metadata[strata_columns].astype(str).agg("||".join, axis=1).to_numpy()
    exceed = np.zeros(len(observed), dtype=int)
    for _ in range(args.permutations):
        null_values = permute_within_strata(clr_values, strata, rng)
        null_rho = rho_p_matrix(null_values)[0][iu]
        exceed += np.abs(null_rho) >= np.abs(observed)
    p_values = (exceed + 1.0) / (args.permutations + 1.0)
    q_values = bh_adjust(p_values)

    features = matrix.index.astype(str).tolist()
    results = pd.DataFrame({
        "feature_1": [features[i] for i in iu[0]],
        "feature_2": [features[j] for j in iu[1]],
        "rho_p": observed,
        "direction": np.where(observed >= 0, "positive", "negative"),
        "variation_of_log_ratio": vlr[iu],
        "bootstrap_median_rho_p": np.median(bootstrap, axis=0),
        "bootstrap_ci_lower": np.quantile(bootstrap, 0.025, axis=0),
        "bootstrap_ci_upper": np.quantile(bootstrap, 0.975, axis=0),
        "bootstrap_recovery": recovery,
        "sign_consistency": sign_consistency,
        "permutation_p_value": p_values,
        "q_value": q_values,
    })
    results["passes_effect_size"] = results["rho_p"].abs().ge(args.min_abs_rho)
    results["passes_bootstrap_recovery"] = results["bootstrap_recovery"].ge(args.min_bootstrap_recovery)
    results["passes_sign_consistency"] = results["sign_consistency"].ge(args.min_sign_consistency)
    results["passes_adjusted_significance"] = results["q_value"].le(args.max_q)
    criteria = [
        "passes_effect_size", "passes_bootstrap_recovery",
        "passes_sign_consistency", "passes_adjusted_significance",
    ]
    results["selected_edge"] = results[criteria].all(axis=1)
    results["failed_criteria"] = results.apply(
        lambda row: ";".join(column.replace("passes_", "") for column in criteria if not truthy(row[column])), axis=1
    )
    results.to_csv(args.outdir / f"{args.prefix}_all_pair_statistics.tsv", sep="\t", index=False)
    selected = results.loc[results["selected_edge"]].copy()
    selected.to_csv(args.outdir / f"{args.prefix}_selected_edges.tsv", sep="\t", index=False)
    selected[["feature_1", "feature_2", "rho_p", "direction", "bootstrap_recovery", "sign_consistency", "q_value"]].to_csv(
        args.outdir / f"{args.prefix}_edge_list.csv", index=False
    )
    network_outputs(features, selected, args.outdir, args.prefix, args.seed)

    summary = pd.DataFrame([{
        "assay_samples": matrix.shape[1],
        "species_tested": matrix.shape[0],
        "pairs_tested": len(results),
        "selected_edges": len(selected),
        "positive_edges": int(selected["direction"].eq("positive").sum()),
        "negative_edges": int(selected["direction"].eq("negative").sum()),
        "isolated_species": int(matrix.shape[0] - len(set(selected["feature_1"]) | set(selected["feature_2"]))),
    }])
    summary.to_csv(args.outdir / f"{args.prefix}_inference_summary.tsv", sep="\t", index=False)
    parameters = {
        "association_metric": "rho_p proportionality",
        "rho_formula": "1 - Var(log(x_i/x_j)) / (Var(clr_i) + Var(clr_j))",
        "zero_replacement": "sample-wise multiplicative replacement",
        "zero_replacement_fraction_of_minimum_positive_part": args.zero_replacement_fraction,
        "minimum_relative_abundance": args.min_relative_abundance,
        "minimum_prevalence": args.min_prevalence,
        "minimum_absolute_rho_p": args.min_abs_rho,
        "bootstrap_iterations": args.bootstrap_iterations,
        "bootstrap_unit": f"whole {args.cruise_col} resampled with replacement within {args.season_col}",
        "minimum_bootstrap_recovery": args.min_bootstrap_recovery,
        "minimum_sign_consistency": args.min_sign_consistency,
        "permutations": args.permutations,
        "permutation_scheme": f"independent feature shuffling within {' x '.join(strata_columns)} strata",
        "maximum_bh_q_value": args.max_q,
        "random_seed": args.seed,
    }
    (args.outdir / f"{args.prefix}_parameters.json").write_text(json.dumps(parameters, indent=2) + "\n")
    pd.DataFrame([
        {"software": "Python", "version": platform.python_version()},
        {"software": "NumPy", "version": np.__version__},
        {"software": "pandas", "version": pd.__version__},
        {"software": "NetworkX", "version": nx.__version__},
    ]).to_csv(args.outdir / f"{args.prefix}_software_versions.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
