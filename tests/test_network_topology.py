from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

import networkx as nx


SCRIPT = Path(__file__).parents[1] / "processes/network_topology/network_topology_stats.py"
SPEC = importlib.util.spec_from_file_location("network_topology_stats", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class NetworkTopologyTest(unittest.TestCase):
    def test_triangle_basic_statistics(self) -> None:
        stats = MODULE.basic_stats(nx.cycle_graph(3))
        self.assertEqual(stats["n_nodes"], 3)
        self.assertEqual(stats["n_edges"], 3)
        self.assertEqual(stats["n_connected_components"], 1)
        self.assertEqual(stats["largest_component_size"], 3)
        self.assertAlmostEqual(stats["transitivity"], 1.0)
        self.assertAlmostEqual(stats["average_clustering"], 1.0)

    def test_null_model_is_seeded_and_writes_machine_readable_tables(self) -> None:
        graph = nx.barabasi_albert_graph(20, 2, seed=7)
        first, draws_a = MODULE.null_comparison(graph, 8, 42)
        second, draws_b = MODULE.null_comparison(graph, 8, 42)
        self.assertEqual(draws_a, draws_b)
        self.assertEqual(first, second)
        self.assertEqual(first["n_null_completed"], 8)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            MODULE.write_key_value_table(root / "summary.tsv", {**MODULE.basic_stats(graph), **first, "n_null_requested": 8, "seed": 42})
            MODULE.write_draws(root / "draws.tsv", draws_a)
            self.assertIn("metric\tvalue", (root / "summary.tsv").read_text())
            self.assertEqual(len((root / "draws.tsv").read_text().splitlines()), 9)


if __name__ == "__main__":
    unittest.main()
