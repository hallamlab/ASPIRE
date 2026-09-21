from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

import pandas as pd


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "processes"
    / "measurement_association"
    / "measurement_association.py"
)
SPEC = importlib.util.spec_from_file_location("measurement_association", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class MeasurementAssociationTaxonomyTests(unittest.TestCase):
    def test_explicit_measurement_table_wins_column_collision(self) -> None:
        metadata = pd.DataFrame(
            {"sampleID": ["S1", "S2"], "Cruise": [1, 1], "Depth": [10, 20], "Oxygen": [999, 999]}
        )
        measurements = pd.DataFrame(
            {"Cruise": [1, 1], "Depth": [10, 20], "Oxygen": [0.0, 4.5]}
        )
        observed = MODULE.merge_measurements(
            metadata,
            measurements,
            "sampleID",
            "sampleID",
            ["Cruise", "Depth"],
            ["Cruise", "Depth"],
        )

        self.assertEqual(observed["Oxygen"].tolist(), [0.0, 4.5])
        self.assertEqual(observed["Oxygen_metadata"].tolist(), [999, 999])

    def test_correlations_are_pairwise_complete_without_imputation(self) -> None:
        asvs = pd.DataFrame({"ASV1": [0.0, 1.0, 2.0, 3.0]}, index=["S1", "S2", "S3", "S4"])
        measurements = pd.DataFrame(
            {"Oxygen": [0.0, 1.0, 2.0, 3.0], "Sparse": [0.0, float("nan"), 2.0, 3.0]},
            index=asvs.index,
        )
        _, observed = MODULE.correlation_tables(asvs, measurements, "both")
        observed = observed.set_index("measurement")

        self.assertEqual(int(observed.loc["Oxygen", "n"]), 4)
        self.assertEqual(int(observed.loc["Sparse", "n"]), 3)
        self.assertAlmostEqual(float(observed.loc["Sparse", "rho"]), 1.0)

    def test_spieceasi_subset_uses_standard_filter_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "spieceasi_filtering_audit.csv"
            pd.DataFrame(
                {
                    "ASV_ID": ["ASV1", "ASV2", "ASV3"],
                    "passes_standard_filter": [True, False, False],
                    "retained_for_asv_mag_link": [False, True, False],
                    "retained_final": [True, True, False],
                }
            ).to_csv(path, index=False)
            observed = MODULE.load_spieceasi_standard_filter(path)

        self.assertEqual(observed, {"ASV1"})

    def test_taxonomy_is_deduplicated_and_missing_values_are_explicit(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "asv_meta.tsv"
            pd.DataFrame(
                {
                    "ASV_ID": ["ASV1", "ASV1", "ASV2"],
                    "Phylum": ["Proteobacteria", "Proteobacteria", ""],
                    "sampleID": ["S1", "S2", "S1"],
                }
            ).to_csv(path, sep="\t", index=False)
            observed = MODULE.load_asv_taxonomy(path, "ASV_ID")

        self.assertEqual(observed["ASV_ID"].tolist(), ["ASV1", "ASV2"])
        self.assertEqual(
            observed["Phylum"].tolist(),
            ["Proteobacteria", "Unclassified"],
        )

    def test_phylum_summary_preserves_asv_level_associations(self) -> None:
        correlations = pd.DataFrame(
            {
                "ASV_ID": ["ASV1", "ASV2", "ASV3"],
                "measurement": ["Oxygen", "Oxygen", "Oxygen"],
                "rho": [0.8, 0.4, -0.5],
                "q_value": [0.01, 0.20, 0.02],
            }
        )
        taxonomy = pd.DataFrame(
            {
                "ASV_ID": ["ASV1", "ASV2", "ASV3"],
                "Phylum": ["A", "A", "B"],
            }
        )
        observed = MODULE.summarize_taxonomic_correlations(correlations, taxonomy)
        phylum_a = observed.set_index("Phylum").loc["A"]
        phylum_b = observed.set_index("Phylum").loc["B"]

        self.assertEqual(int(phylum_a["n_asvs"]), 2)
        self.assertAlmostEqual(float(phylum_a["median_rho"]), 0.6)
        self.assertAlmostEqual(float(phylum_a["fraction_moderate_significant"]), 0.5)
        self.assertAlmostEqual(float(phylum_b["fraction_moderate_significant"]), 1.0)


if __name__ == "__main__":
    unittest.main()
