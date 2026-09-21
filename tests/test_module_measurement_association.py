from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "processes" / "measurement_association" / "module_measurement_association.py"
SPEC = importlib.util.spec_from_file_location("module_measurement_association", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class ModuleMeasurementAssociationTests(unittest.TestCase):
    def test_anchor_taxonomy_parses_prefixed_silva_ranks(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "taxonomy.tsv"
            pd.DataFrame(
                {
                    "Feature ID": ["ASV1;size=12"],
                    "Taxon": [
                        "d__Bacteria; p__Proteobacteria; f__Nitrosomonadaceae; g__Nitrosomonas"
                    ],
                }
            ).to_csv(path, sep="\t", index=False)
            observed = MODULE.load_anchor_taxonomy(path).iloc[0]

        self.assertEqual(observed["ASV_ID"], "ASV1")
        self.assertEqual(observed["asv_phylum"], "Proteobacteria")
        self.assertEqual(observed["asv_family"], "Nitrosomonadaceae")
        self.assertEqual(observed["asv_genus"], "Nitrosomonas")

    def test_anchor_taxonomy_display_cleans_raw_silva_values(self) -> None:
        row = type(
            "TaxonomyRow",
            (),
            {
                "ASV_ID": "ASV12",
                "asv_phylum": "Bacteroidota",
                "asv_family": "Flavobacteriaceae",
                "asv_genus": "NS5_marine_group",
            },
        )()
        self.assertEqual(
            MODULE.anchor_taxonomy_label(row),
            "ASV12 — P: Bacteroidota; F: Flavobacteriaceae; G: NS5 marine group",
        )
        self.assertEqual(MODULE.clean_taxonomy_display("uncultured"), "unclassified")

    def test_network_anchor_is_top_rank_in_any_centrality_metric(self) -> None:
        membership = pd.DataFrame(
            {"ASV_ID": ["A", "B", "C", "D"], "module_label": ["M1"] * 4}
        )
        features = pd.DataFrame(
            {
                "Taxon": ["A", "B", "C", "D"],
                "Degree": [10, 3, 2, 1],
                "EigenCentral": [0.1, 0.9, 0.2, 0.0],
                "Betweenness": [0.1, 0.2, 0.8, 0.0],
            }
        )
        observed = MODULE.network_anchor_table(features, membership, 1).set_index("ASV_ID")

        self.assertTrue(bool(observed.loc["A", "is_ecological_anchor"]))
        self.assertTrue(bool(observed.loc["B", "is_ecological_anchor"]))
        self.assertTrue(bool(observed.loc["C", "is_ecological_anchor"]))
        self.assertFalse(bool(observed.loc["D", "is_ecological_anchor"]))

    def test_global_network_metrics_separate_prevalence_and_connector_roles(self) -> None:
        if MODULE.nx is None:
            self.skipTest("networkx unavailable")
        nodes = pd.DataFrame(
            {
                "ASV_ID": ["A", "B", "C", "D"],
                "module_label": ["M1", "M1", "M2", "M2"],
                "Degree": [2, 1, 2, 1],
                "Degree_thresholded": [2, 1, 1, 1],
                "Betweenness": [0.5, 0.0, 0.5, 0.0],
                "Betweenness_raw": [1.0, 0.0, 1.0, 0.0],
                "Betweenness_norm": [0.5, 0.0, 0.5, 0.0],
                "Closeness": [0.8, 0.5, 0.8, 0.5],
                "EigenCentral": [0.8, 0.2, 0.7, 0.1],
            }
        )
        graph = MODULE.nx.Graph()
        graph.add_edge("A", "B", weight=0.2)
        graph.add_edge("A", "C", weight=0.5)
        graph.add_edge("C", "D", weight=0.3)
        raw = pd.DataFrame(
            {"S1": [2, 1, 0, 0], "S2": [3, 0, 1, 0]},
            index=["A", "B", "C", "D"],
        )
        observed = MODULE.add_global_network_metrics(nodes, graph, raw).set_index("ASV_ID")

        self.assertEqual(observed.loc["A", "prevalence"], 1.0)
        self.assertEqual(observed.loc["B", "prevalence"], 0.5)
        self.assertAlmostEqual(observed.loc["A", "Strength"], 0.7)
        self.assertEqual(int(observed.loc["A", "cross_module_degree"]), 1)
        self.assertAlmostEqual(observed.loc["A", "Participation"], 0.5)
        self.assertEqual(int(observed.loc["A", "global_rank_prevalence"]), 1)
        self.assertIn("global_rank_closeness", observed.columns)
        self.assertIn("global_percentile_within_module_degree_z", observed.columns)

    def test_all_informative_node_metrics_receive_locked_layout_plots(self) -> None:
        expected = {
            "prevalence", "mean_relative_abundance", "max_relative_abundance",
            "observed_samples", "degree", "degree_thresholded", "strength",
            "within_module_degree", "within_module_degree_z",
            "cross_module_degree", "betweenness", "betweenness_raw",
            "closeness", "eigenvector", "participation",
        }
        self.assertEqual(set(MODULE.NETWORK_METRIC_SPECS), expected)
        self.assertNotIn("assessed_samples_n", MODULE.NETWORK_METRIC_SPECS)
        self.assertNotIn("Betweenness_norm", MODULE.NETWORK_METRIC_SPECS)

    def test_metric_network_uses_locked_network_analysis_coordinates(self) -> None:
        nodes = pd.DataFrame(
            {
                "ASV_ID": ["A", "B", "C"],
                "ecological_module": ["M1", "M1", "M2"],
            }
        )
        with tempfile.TemporaryDirectory() as tmp:
            layout_path = Path(tmp) / "layout.tsv"
            pd.DataFrame(
                {
                    "ASV_ID": ["C", "A", "B"],
                    "spieceasi_x": [3.5, -1.0, 0.25],
                    "spieceasi_y": [2.0, 4.25, -0.5],
                }
            ).to_csv(layout_path, sep="\t", index=False)
            observed = MODULE.locked_spieceasi_network_layout(
                nodes, layout_path
            ).set_index("ASV_ID")

        self.assertEqual(observed.loc["A", "spieceasi_x"], -1.0)
        self.assertEqual(observed.loc["B", "spieceasi_y"], -0.5)
        self.assertEqual(observed.loc["C", "spieceasi_x"], 3.5)

    def test_environmental_network_callouts_are_external_and_nonoverlapping(self) -> None:
        original = np.array([
            [-1.0, -0.6],
            [0.9, -0.4],
            [0.8, 0.7],
            [-0.7, 0.8],
        ])
        radii = np.array([0.45, 0.55, 0.50, 0.48])
        background = np.array([
            [-2.0, -1.0],
            [2.0, -1.0],
            [2.0, 1.0],
            [-2.0, 1.0],
            [0.0, 0.0],
        ])
        displayed = MODULE.perimeter_module_centers(original, radii, background)

        center = np.array([0.0, 0.0])
        half_span = np.array([2.0, 1.0])
        for index, point in enumerate(displayed):
            beyond_x = abs(point[0] - center[0]) > half_span[0] + radii[index]
            beyond_y = abs(point[1] - center[1]) > half_span[1] + radii[index]
            self.assertTrue(beyond_x or beyond_y)
        for i in range(len(displayed)):
            for j in range(i + 1, len(displayed)):
                self.assertGreaterEqual(
                    np.linalg.norm(displayed[i] - displayed[j]),
                    radii[i] + radii[j] + 0.14 - 1e-9,
                )

    def test_module_scores_use_total_library_relative_abundance(self) -> None:
        counts = pd.DataFrame(
            {"S1": [10, 20, 70], "S2": [0, 25, 75]},
            index=["ASV1", "ASV2", "ASV3"],
        )
        membership = pd.DataFrame(
            {"ASV_ID": ["ASV1", "ASV2"], "module_label": ["M1", "M1"]}
        )
        observed = MODULE.module_scores(counts, membership).set_index("sampleID")

        self.assertAlmostEqual(observed.loc["S1", "module_relative_abundance"], 0.30)
        self.assertAlmostEqual(observed.loc["S2", "module_relative_abundance"], 0.25)
        self.assertEqual(int(observed.loc["S1", "module_asvs_n"]), 2)

    def test_module_correlations_are_pairwise_complete_and_keep_zero(self) -> None:
        scores = pd.DataFrame(
            {
                "sampleID": ["S1", "S2", "S3", "S4"],
                "ecological_module": ["M1"] * 4,
                "module_relative_abundance": [0.0, 0.1, 0.2, 0.3],
            }
        )
        measurements = pd.DataFrame(
            {"Core": [0.0, 1.0, 2.0, 3.0], "Sparse": [0.0, np.nan, 2.0, 3.0]},
            index=["S1", "S2", "S3", "S4"],
        )
        metadata = pd.DataFrame(
            {
                "sampleID": ["S1", "S2", "S3", "S4"],
                "Cruise": [1, 1, 2, 2],
                "Depth": [10, 20, 10, 20],
            }
        )
        observed = MODULE.module_measurement_tests(
            scores, measurements, metadata, "sampleID", "Cruise", "Depth", {"Sparse"}
        ).set_index("measurement")

        self.assertEqual(int(observed.loc["Core", "samples_measured"]), 4)
        self.assertEqual(int(observed.loc["Sparse", "samples_measured"]), 3)
        self.assertEqual(observed.loc["Core", "measurement_class"], "core")
        self.assertEqual(observed.loc["Sparse", "measurement_class"], "sparse")
        self.assertAlmostEqual(float(observed.loc["Core", "rho"]), 1.0)
        self.assertAlmostEqual(float(observed.loc["Sparse", "rho"]), 1.0)

    def test_coverage_keeps_configured_but_unavailable_measurements(self) -> None:
        coverage = pd.DataFrame(
            {
                "measurement": ["Oxygen"],
                "measurement_class": ["core"],
                "samples_total": [4],
                "samples_measured": [4],
                "samples_missing": [0],
                "measured_fraction": [1.0],
                "cruises_measured": [2],
                "depths_measured": [2],
                "minimum_depth": [10.0],
                "maximum_depth": [20.0],
            }
        )
        audit = pd.DataFrame(
            {
                "measurement": ["Oxygen", "Iron"],
                "overlapping_samples": [4, 4],
                "measured_samples": [4, 0],
                "exclusion_reason": ["retained", "fewer_than_three_measured_samples"],
            }
        )
        observed = MODULE.add_unavailable_measurements_to_coverage(
            coverage, audit, {"Iron"}
        ).set_index("measurement")

        self.assertEqual(observed.loc["Iron", "measurement_class"], "sparse")
        self.assertEqual(int(observed.loc["Iron", "samples_measured"]), 0)
        self.assertEqual(
            observed.loc["Iron", "availability_status"],
            "fewer_than_three_measured_samples",
        )

    def test_workflow_exposes_module_analysis_and_no_imputation_contract(self) -> None:
        pipeline = (ROOT / "asv_pipeline.nf").read_text()
        workflow = (ROOT / "workflow" / "modules" / "indicators_diagnostics.nf").read_text()
        network_workflow = (
            ROOT / "workflow" / "modules" / "networks_reporting.nf"
        ).read_text()
        measurement = (
            ROOT / "processes" / "measurement_association" / "measurement_association.py"
        ).read_text()

        self.assertIn("MODULE_MEASUREMENT_ASSOCIATION(", pipeline)
        self.assertIn(
            "moduleMeasurementNetworkMetricTopN: moduleMeasurementNetworkMetricTopN",
            pipeline,
        )
        self.assertIn("process MODULE_MEASUREMENT_ASSOCIATION", workflow)
        self.assertIn('path(network_graph)', workflow)
        self.assertIn('path(spieceasi_layout)', workflow)
        self.assertIn('path(node_features)', workflow)
        self.assertIn('path(taxonomy_table)', workflow)
        self.assertIn('path(prevalence_counts)', workflow)
        self.assertIn('val(network_metric_top_n)', workflow)
        self.assertIn('--prevalence-counts "${prevalence_counts}"', workflow)
        self.assertIn('--spieceasi-layout "${spieceasi_layout}"', workflow)
        self.assertIn('--taxonomy "${taxonomy_table}"', workflow)
        self.assertIn('--anchor-top-n ${asvMagNetworkAnchorTopN}', workflow)
        self.assertIn('--network-metric-top-n ${network_metric_top_n}', workflow)
        self.assertIn("graph_network_stage.layout_all", pipeline)
        self.assertIn("layout_all = stage.layout_all", pipeline)
        self.assertIn(
            'path("spieceasi_network_layout_all.tsv"), emit: layout_all',
            network_workflow,
        )
        self.assertIn("ordination_measurements.notna().all(axis=1)", measurement)
        self.assertNotIn("measurement_matrix.fillna(measurement_matrix.median", measurement)


if __name__ == "__main__":
    unittest.main()
