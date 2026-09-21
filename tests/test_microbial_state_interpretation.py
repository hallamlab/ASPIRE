from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT = (
    Path(__file__).parents[1]
    / "processes"
    / "microbial_state_interpretation"
    / "microbial_state_interpretation.py"
)
SPEC = importlib.util.spec_from_file_location("microbial_state_interpretation", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class MicrobialStateInterpretationTest(unittest.TestCase):
    def test_assignment_merge_reconciles_identical_staged_columns(self) -> None:
        assignments = pd.DataFrame({
            "sampleID": ["S1", "S2"],
            "microbial_compartment": ["MC1", "MC2"],
            "selected_K": [2, 2],
        })
        metadata = pd.DataFrame({
            "sampleID": ["S1", "S2"],
            "microbial_compartment": ["MC1", "MC2"],
            "selected_K": [2.0, 2.0],
            "Depth": [10, 100],
        })
        merged = MODULE.merge_assignments_with_metadata(assignments, metadata, "sampleID")
        self.assertEqual(merged.columns.tolist().count("microbial_compartment"), 1)
        self.assertEqual(merged.microbial_compartment.tolist(), ["MC1", "MC2"])
        self.assertEqual(merged.Depth.tolist(), [10, 100])

    def test_assignment_merge_rejects_conflicting_staged_columns(self) -> None:
        assignments = pd.DataFrame({
            "sampleID": ["S1", "S2"],
            "microbial_compartment": ["MC1", "MC2"],
        })
        metadata = pd.DataFrame({
            "sampleID": ["S1", "S2"],
            "microbial_compartment": ["MC1", "MC3"],
        })
        with self.assertRaisesRegex(ValueError, "Conflicting staged values"):
            MODULE.merge_assignments_with_metadata(assignments, metadata, "sampleID")

    def test_renewal_onsets_are_loaded_without_duplicate_dates(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "renewals.tsv"
            pd.DataFrame({"start_date": ["2012-01-01", "2012-01-01", "2013-01-01"]}).to_csv(
                path, sep="\t", index=False
            )
            onsets = MODULE.load_renewal_onsets(path, "start_date")
        self.assertEqual(onsets.tolist(), [pd.Timestamp("2012-01-01"), pd.Timestamp("2013-01-01")])

    def test_balanced_accuracy_uses_only_classes_present_in_fold(self) -> None:
        truth = np.array(["MC1", "MC1", "MC2", "MC2"])
        predicted = np.array(["MC1", "MC1", "MC1", "MC2"])
        self.assertAlmostEqual(
            MODULE.observed_class_balanced_accuracy(truth, predicted), 0.75
        )

    def test_identical_labels_have_perfect_correspondence(self) -> None:
        values = pd.Series(["A", "A", "B", "B", "C", "C"])
        metrics = MODULE.association_metrics(values, values)
        self.assertAlmostEqual(metrics["adjusted_rand_index"], 1.0)
        self.assertAlmostEqual(metrics["normalized_mutual_information"], 1.0)
        self.assertAlmostEqual(metrics["majority_mapping_accuracy"], 1.0)
        self.assertAlmostEqual(
            metrics["normalized_conditional_entropy_mc_given_environment"], 0.0
        )

    def test_cruise_restricted_null_preserves_block_label_counts(self) -> None:
        frame = pd.DataFrame({
            "mc": ["MC1", "MC2", "MC1", "MC2", "MC1", "MC2"],
            "environment": ["A", "A", "B", "A", "B", "B"],
            "cruise": ["C1", "C1", "C1", "C2", "C2", "C2"],
        })
        observed, null = MODULE.restricted_permutations(
            frame, "mc", "environment", "cruise", 9, np.random.default_rng(42)
        )
        self.assertEqual(len(null), 9)
        self.assertTrue(0 < observed["adjusted_rand_index_permutation_pvalue"] <= 1)

    def test_renewal_test_uses_whole_cruise_labels_within_year(self) -> None:
        frame = pd.DataFrame({
            "phase": ["baseline"] * 4 + ["renewal"] * 4 + ["baseline"] * 4 + ["renewal"] * 4,
            "Depth": [10, 100, 10, 100] * 4,
            "Cruise": np.repeat(["C1", "C2", "C3", "C4"], 4),
            "Year": np.repeat([2020, 2020, 2021, 2021], 4),
        })
        response = np.column_stack([
            frame.phase.eq("baseline").to_numpy(float),
            frame.phase.eq("renewal").to_numpy(float),
        ])
        result = MODULE.cruise_label_depth_adjusted_test(
            response, frame, "phase", "Depth", "Cruise", "Year", 19,
            np.random.default_rng(42),
        )
        self.assertEqual(result["permutation_scheme"], "cruise_labels_within_calendar_year")
        self.assertGreater(result["nonidentity_permutations_n"], 0)
        self.assertEqual(result["cruises_n"], 4)
        self.assertGreater(result["incremental_r_squared"], 0.9)

    def test_basin_comparable_surface_retains_observed_mc_anchors(self) -> None:
        frame = pd.DataFrame({
            "sampleID": ["S1", "S2", "S3", "S4"],
            "Date": pd.to_datetime(["2020-01-01", "2020-01-01", "2020-02-01", "2020-02-01"]),
            "Depth": [10.0, 100.0, 10.0, 100.0],
            "microbial_compartment": ["MC1", "MC2", "MC1", "MC2"],
        })
        surface = MODULE.build_basin_comparable_mc_surface(
            frame,
            ["MC1", "MC2"],
            "Date",
            "Depth",
            110.0,
            1.0,
            7,
            90.0,
            30.0,
        )
        dates = pd.DatetimeIndex(surface["dates"])
        depths = np.asarray(surface["depths"])
        winner = surface["winner"]
        for _, row in frame.iterrows():
            date_index = dates.get_loc(row.Date)
            depth_index = int(np.where(np.isclose(depths, row.Depth))[0][0])
            expected = ["MC1", "MC2"].index(row.microbial_compartment)
            self.assertFalse(bool(winner.mask[depth_index, date_index]))
            self.assertEqual(int(winner[depth_index, date_index]), expected)

    def test_basin_reference_grid_fixes_complete_date_depth_domain(self) -> None:
        frame = pd.DataFrame({
            "sampleID": ["S1", "S2", "S3", "S4"],
            "Date": pd.to_datetime(["2020-01-01", "2020-01-01", "2020-02-01", "2020-02-01"]),
            "Depth": [10.0, 100.0, 10.0, 100.0],
            "microbial_compartment": ["MC1", "MC2", "MC1", "MC2"],
        })
        with tempfile.TemporaryDirectory() as tmpdir:
            reference_path = Path(tmpdir) / "basin_grid.tsv"
            pd.DataFrame({
                "date": np.repeat(pd.to_datetime(["2020-01-01", "2021-01-01"]), 3),
                "depth_m": np.tile([0.0, 100.0, 210.0], 2),
            }).to_csv(reference_path, sep="\t", index=False)
            surface = MODULE.build_basin_comparable_mc_surface(
                frame, ["MC1", "MC2"], "Date", "Depth", 210.0, 1.0,
                7, 90.0, 30.0, reference_path,
            )
        self.assertEqual(
            pd.DatetimeIndex(surface["dates"]).tolist(),
            [pd.Timestamp("2020-01-01"), pd.Timestamp("2021-01-01")],
        )
        np.testing.assert_array_equal(surface["depths"], [0.0, 100.0, 210.0])


if __name__ == "__main__":
    unittest.main()
