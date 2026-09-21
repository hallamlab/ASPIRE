from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

nx = None
try:
    import networkx as nx
except ImportError:  # pragma: no cover - exercised only outside the module env.
    nx = None


SCRIPT = Path(__file__).parents[1] / "processes" / "asv_mag_network" / "asv_mag_network.py"
SPEC = importlib.util.spec_from_file_location("asv_mag_network", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
try:
    SPEC.loader.exec_module(MODULE)
except ImportError as exc:  # pragma: no cover - exercised only outside the module env.
    MODULE = None
    MODULE_IMPORT_ERROR = exc
else:
    MODULE_IMPORT_ERROR = None


class AsvMagNetworkTest(unittest.TestCase):
    @unittest.skipIf(MODULE is None, "network dependencies are supplied by the module env")
    def test_mag_species_use_shapes_and_black_nodes_use_light_rims(self) -> None:
        self.assertEqual(
            MODULE.species_marker_map(["Thioglobus", "SUP05"]),
            {"SUP05": "o", "Thioglobus": "s"},
        )
        self.assertEqual(MODULE.species_marker_map(["SUP05"]), {"SUP05": "o"})
        self.assertEqual(MODULE.node_outline_color("#000000"), "#E6E6E6")
        self.assertEqual(MODULE.node_outline_color("#BDBDBD"), "#000000")

    def test_shared_anchor_top_n_is_exported_to_modular_processes(self) -> None:
        root = Path(__file__).resolve().parents[1]
        pipeline = (root / "asv_pipeline.nf").read_text()
        workflow = (root / "workflow" / "modules" / "networks_reporting.nf").read_text()
        config = (root / "examples" / "si.local.yml").read_text()

        self.assertIn("asvMagNetworkAnchorTopN: asvMagNetworkAnchorTopN", pipeline)
        self.assertIn("--anchor-top-n ${asvMagNetworkAnchorTopN}", workflow)
        self.assertIn("--top-n ${asvMagNetworkAnchorTopN}", workflow)
        self.assertIn("anchor_top_n: 1", config)

    def test_vsearch_size_suffix_is_removed_without_changing_base_asv_id(self) -> None:
        if MODULE is None:
            self.skipTest("networkx/seaborn dependencies are supplied by the module env")
        self.assertEqual(MODULE.standardize_asv_id("ASV1;size=7176726"), "ASV1")
        self.assertEqual(MODULE.standardize_asv_id("ASV_1"), "ASV_1")

    def test_suffix_after_double_underscore_mag_id_mode(self) -> None:
        if MODULE is None:
            self.skipTest("networkx/seaborn dependencies are supplied by the module env")
        self.assertEqual(
            MODULE.standardize_mag_id("100m__xPG_SAGs__AB-746_B09_AB-901.fasta.gz", "suffix_after_double_underscore"),
            "AB-746_B09_AB-901",
        )

    def test_all_ecological_modules_can_be_displayed_despite_best_selection(self) -> None:
        if MODULE is None:
            self.skipTest("networkx/seaborn dependencies are supplied by the module env")
        node_features = pd.DataFrame(
            {
                "ASV_ID": ["ASV1", "ASV2", "ASV3"],
                "Degree": [3, 2, 1],
                "EigenCentral": [0.3, 0.2, 0.1],
                "Betweenness": [0.03, 0.02, 0.01],
            }
        )
        modules = pd.DataFrame(
            {
                "Taxon": ["ASV1", "ASV2", "ASV3"],
                "module_label": ["M1", "M2", "M3"],
            }
        )
        selection = pd.DataFrame(
            {
                "module_label": ["M1", "M2", "M3"],
                "is_best": [True, False, False],
            }
        )

        selected = MODULE.annotate_ecological_modules(
            node_features, modules, 1, selection
        )
        self.assertEqual(
            set(selected["display_module_label"]), {"M1", "Other modules"}
        )

        displayed_all = MODULE.annotate_ecological_modules(
            node_features, modules, 1, selection, display_all_modules=True
        )
        self.assertEqual(
            set(displayed_all["display_module_label"]), {"M1", "M2", "M3"}
        )
        self.assertTrue(displayed_all["is_display_ecological_module"].all())
        self.assertEqual(
            MODULE.standardize_mag_id("xPG_SAGs__100m__AB-746_B09_AB-901", "suffix_after_double_underscore"),
            "AB-746_B09_AB-901",
        )

    def test_ecological_anchor_qualifies_by_any_single_centrality_metric(self) -> None:
        if MODULE is None:
            self.skipTest("networkx/seaborn dependencies are supplied by the module env")
        node_features = pd.DataFrame(
            {
                "ASV_ID": ["ASV_degree", "ASV_eigen", "ASV_between", "ASV_other"],
                "Degree": [10, 3, 2, 1],
                "EigenCentral": [0.1, 0.9, 0.2, 0.0],
                "Betweenness": [0.1, 0.2, 0.8, 0.0],
            }
        )
        modules = pd.DataFrame(
            {
                "Taxon": node_features["ASV_ID"],
                "module_label": ["M1", "M1", "M1", "M1"],
            }
        )

        annotated = MODULE.annotate_ecological_modules(
            node_features, modules, anchor_top_n=1, display_all_modules=True
        ).set_index("ASV_ID")

        self.assertTrue(annotated.loc["ASV_degree", "is_ecological_anchor"])
        self.assertTrue(annotated.loc["ASV_eigen", "is_ecological_anchor"])
        self.assertTrue(annotated.loc["ASV_between", "is_ecological_anchor"])
        self.assertFalse(annotated.loc["ASV_other", "is_ecological_anchor"])

    @unittest.skipIf(nx is None or MODULE is None, "networkx/seaborn dependencies are supplied by the module env")
    def test_smoke_exports_paper_and_heterogeneous_networks(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            graph_path = root / "network.graphml"
            graph = nx.Graph()
            graph.add_edge("ASV1", "ASV2", weight=0.42)
            nx.write_graphml(graph, graph_path)

            node_features = root / "node_features.csv"
            pd.DataFrame(
                [
                    {"Taxon": "ASV1", "degree": 1},
                    {"Taxon": "ASV2", "degree": 1},
                ]
            ).to_csv(node_features, index=False)

            taxonomy = root / "taxonomy.tsv"
            pd.DataFrame(
                [
                    {"Feature ID": "ASV1", "Taxon": "d__Bacteria;p__Proteobacteria;f__Nitrosomonadaceae;g__Nitrosomonas"},
                    {"Feature ID": "ASV2", "Taxon": "d__Bacteria;p__Bacteroidota;g__Flavobacterium"},
                ]
            ).to_csv(taxonomy, sep="\t", index=False)

            pairing = root / "pairing.tsv"
            pd.DataFrame(
                [
                    {
                        "ASV_ID": "ASV1",
                        "genome_id": "MAG_A",
                        "pairing_status": "paired_unique",
                        "link_pident": 100.0,
                        "link_qcov": 100.0,
                        "mag_phylum": "Proteobacteria",
                        "mag_genus": "Nitrosomonas",
                    },
                    {
                        "ASV_ID": "ASV2",
                        "genome_id": "MAG_B",
                        "pairing_status": "paired_unique",
                        "link_pident": 100.0,
                        "link_qcov": 100.0,
                        "mag_phylum": "Firmicutes",
                        "mag_genus": "Bacillus",
                    },
                ]
            ).to_csv(pairing, sep="\t", index=False)

            outdir = root / "out"
            old_argv = sys.argv
            try:
                sys.argv = [
                    "asv_mag_network.py",
                    "--graph",
                    str(graph_path),
                    "--node-features",
                    str(node_features),
                    "--asv-mag-pairing",
                    str(pairing),
                    "--taxonomy",
                    str(taxonomy),
                    "--outdir",
                    str(outdir),
                    "--prefix",
                    "smoke",
                    "--asv-taxonomy-source",
                    "ncbi",
                    "--mag-taxonomy-source",
                    "ncbi",
                ]
                MODULE.main()
            finally:
                sys.argv = old_argv

            validation = pd.read_csv(outdir / "validation" / "smoke_taxonomy_validation.tsv", sep="\t")
            accepted = validation.set_index("ASV_ID")["accepted_paper_pair"].to_dict()
            self.assertTrue(bool(accepted["ASV1"]))
            self.assertFalse(bool(accepted["ASV2"]))
            self.assertEqual(
                validation.set_index("ASV_ID").loc["ASV2", "taxonomy_validation_status"],
                "taxonomy_rejected",
            )
            self.assertTrue((outdir / "network" / "smoke_paper.graphml").exists())
            self.assertTrue((outdir / "network" / "smoke_metagenome_heterogeneous.graphml").exists())
            hetero_nodes = pd.read_csv(
                outdir / "network" / "smoke_metagenome_heterogeneous_nodes.tsv",
                sep="\t",
            )
            self.assertIn("MAG_A", set(hetero_nodes["id"].astype(str)))
            self.assertNotIn("MAG_B", set(hetero_nodes["id"].astype(str)))
            asv1 = hetero_nodes.set_index("id").loc["ASV1"]
            self.assertEqual(asv1["asv_phylum"], "Proteobacteria")
            self.assertEqual(asv1["asv_family"], "Nitrosomonadaceae")
            self.assertEqual(asv1["asv_genus"], "Nitrosomonas")


if __name__ == "__main__":
    unittest.main()
