#!/usr/bin/env python3
"""Summarize a SPIEC-EASI network and compare clustering with a null model.

This is the reusable ASPIRE implementation of the network-topology analysis
deposited with the SPARK supplementary code.  The null model preserves the
observed degree sequence with NetworkX's configuration model, simplifies the
result to an undirected graph, and removes self-loops before calculating
clustering statistics.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import networkx as nx
import numpy as np


SUMMARY_COLUMNS = (
    "n_nodes",
    "n_edges",
    "density",
    "average_degree",
    "transitivity",
    "average_clustering",
    "weighted_average_clustering",
    "n_connected_components",
    "largest_component_size",
    "component_sizes",
    "null_mean_transitivity",
    "null_mean_average_clustering",
    "fold_over_null_transitivity",
    "fold_over_null_average_clustering",
    "empirical_p_transitivity",
    "empirical_p_average_clustering",
    "n_null_requested",
    "n_null_completed",
    "seed",
)


def load_graph(path: Path) -> nx.Graph:
    graph = nx.Graph(nx.read_graphml(path))
    graph.remove_edges_from(nx.selfloop_edges(graph))
    return graph


def basic_stats(graph: nx.Graph) -> dict[str, object]:
    degrees = [degree for _, degree in graph.degree()]
    components = sorted(
        (len(component) for component in nx.connected_components(graph)), reverse=True
    )
    try:
        weighted = nx.average_clustering(graph, weight="weight")
    except Exception:
        weighted = math.nan
    return {
        "n_nodes": graph.number_of_nodes(),
        "n_edges": graph.number_of_edges(),
        "density": nx.density(graph),
        "average_degree": float(np.mean(degrees)) if degrees else 0.0,
        "transitivity": nx.transitivity(graph),
        "average_clustering": nx.average_clustering(graph),
        "weighted_average_clustering": weighted,
        "n_connected_components": len(components),
        "largest_component_size": components[0] if components else 0,
        "component_sizes": ",".join(map(str, components)),
    }


def null_comparison(
    graph: nx.Graph, n_null: int, seed: int
) -> tuple[dict[str, object], list[dict[str, object]]]:
    observed_transitivity = nx.transitivity(graph)
    observed_clustering = nx.average_clustering(graph)
    degrees = [degree for _, degree in graph.degree()]
    rng = np.random.default_rng(seed)
    draws: list[dict[str, object]] = []

    for index in range(n_null):
        draw_seed = int(rng.integers(0, 2**31 - 1))
        randomized = nx.Graph(nx.configuration_model(degrees, seed=draw_seed))
        randomized.remove_edges_from(nx.selfloop_edges(randomized))
        draws.append(
            {
                "draw": index + 1,
                "seed": draw_seed,
                "n_edges": randomized.number_of_edges(),
                "transitivity": nx.transitivity(randomized),
                "average_clustering": nx.average_clustering(randomized),
            }
        )
        if (index + 1) % 100 == 0:
            print(f"[INFO] null draws complete: {index + 1}/{n_null}", file=sys.stderr)

    null_transitivity = np.asarray([row["transitivity"] for row in draws], dtype=float)
    null_clustering = np.asarray([row["average_clustering"] for row in draws], dtype=float)
    mean_transitivity = float(null_transitivity.mean()) if len(draws) else math.nan
    mean_clustering = float(null_clustering.mean()) if len(draws) else math.nan
    summary = {
        "null_mean_transitivity": mean_transitivity,
        "null_mean_average_clustering": mean_clustering,
        "fold_over_null_transitivity": (
            observed_transitivity / mean_transitivity if mean_transitivity > 0 else math.nan
        ),
        "fold_over_null_average_clustering": (
            observed_clustering / mean_clustering if mean_clustering > 0 else math.nan
        ),
        "empirical_p_transitivity": (
            float(np.mean(null_transitivity >= observed_transitivity))
            if len(draws)
            else math.nan
        ),
        "empirical_p_average_clustering": (
            float(np.mean(null_clustering >= observed_clustering))
            if len(draws)
            else math.nan
        ),
        "n_null_completed": len(draws),
    }
    return summary, draws


def write_key_value_table(path: Path, values: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("metric\tvalue\n")
        for key in SUMMARY_COLUMNS:
            value = values.get(key, math.nan)
            handle.write(f"{key}\t{value}\n")


def write_draws(path: Path, draws: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = ("draw", "seed", "n_edges", "transitivity", "average_clustering")
    with path.open("w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in draws:
            handle.write("\t".join(str(row[column]) for column in columns) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graphml", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--null-output", type=Path)
    parser.add_argument("--n-null", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--skip-null", action="store_true")
    args = parser.parse_args()
    if args.n_null < 0:
        parser.error("--n-null must be non-negative")
    if not args.graphml.is_file():
        parser.error(f"GraphML file not found: {args.graphml}")
    return args


def main() -> None:
    args = parse_args()
    graph = load_graph(args.graphml)
    print(
        f"[INFO] loaded {graph.number_of_nodes()} nodes and "
        f"{graph.number_of_edges()} edges from {args.graphml}",
        file=sys.stderr,
    )
    summary = basic_stats(graph)
    draws: list[dict[str, object]] = []
    if args.skip_null or args.n_null == 0:
        summary.update(
            {
                "null_mean_transitivity": math.nan,
                "null_mean_average_clustering": math.nan,
                "fold_over_null_transitivity": math.nan,
                "fold_over_null_average_clustering": math.nan,
                "empirical_p_transitivity": math.nan,
                "empirical_p_average_clustering": math.nan,
                "n_null_completed": 0,
            }
        )
    else:
        null_summary, draws = null_comparison(graph, args.n_null, args.seed)
        summary.update(null_summary)
    summary["n_null_requested"] = args.n_null
    summary["seed"] = args.seed
    write_key_value_table(args.output, summary)
    if args.null_output is not None:
        write_draws(args.null_output, draws)


if __name__ == "__main__":
    main()
