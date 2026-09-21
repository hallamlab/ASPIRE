from __future__ import annotations

import importlib.util
import re
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "processes" / "asv_time_depth_curtain" / "asv_time_depth_curtain.py"
SPEC = importlib.util.spec_from_file_location("asv_time_depth_curtain", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class AsvTimeDepthCurtainTests(unittest.TestCase):
    def test_modal_profile_is_deterministic_for_ties(self) -> None:
        data = pd.DataFrame({
            "cruise_index": [0, 0, 0, 0],
            "Depth": [10, 10, 100, 100],
            "hybrid": ["state_b", "state_a", "state_b", "state_b"],
        })
        observed = MODULE.modal_profile(
            data, "Cruise", "Depth", "hybrid", ["state_a", "state_b"]
        ).set_index("depth_m")

        self.assertEqual(observed.loc[10.0, "hybrid_compartment"], "state_a")
        self.assertTrue(bool(observed.loc[10.0, "selection_tied"]))
        self.assertEqual(observed.loc[100.0, "hybrid_compartment"], "state_b")

    def test_smoothed_surface_uses_only_supplied_states(self) -> None:
        profile = pd.DataFrame({
            "cruise_index": [0, 0, 1, 1],
            "depth_m": [10.0, 200.0, 10.0, 200.0],
            "hybrid_compartment": ["oxic__gmm0", "suboxic__gmm4"] * 2,
        })
        labels = ["oxic__gmm0", "suboxic__gmm4"]
        cruise_dates = pd.Series(pd.to_datetime(["2020-02-01", "2020-11-01"]))
        grid, fine_x, depth, codes, support, fine_dates, start, end = MODULE.build_surface(
            profile, cruise_dates, labels, 200.0, 5.0, 2, 0.5
        )

        self.assertEqual(codes.shape, (len(fine_x), len(depth)))
        self.assertEqual(support.shape, (len(labels), len(fine_x), len(depth)))
        self.assertEqual(set(grid["hybrid_compartment_display"]), set(labels))
        self.assertTrue(np.isfinite(grid["maximum_smoothed_support_fraction"]).all())
        self.assertEqual(start, pd.Timestamp("2020-01-01"))
        self.assertEqual(end, pd.Timestamp("2021-01-01"))
        self.assertEqual(grid["calendar_month"].nunique(), 12)
        self.assertEqual(len(fine_dates), 24)

        anchored, anchor_mask, anchor_owner = MODULE.enforce_observed_anchors(
            codes,
            fine_x,
            depth,
            profile,
            cruise_dates,
            labels,
            start,
        )
        self.assertGreater(int(anchor_mask.sum()), 0)
        self.assertTrue((anchor_owner[anchor_mask] >= 0).all())
        for row in profile.itertuples(index=False):
            sample_x = float((cruise_dates.iloc[row.cruise_index] - start).days)
            time_index = int(np.argmin(np.abs(fine_x - sample_x)))
            depth_index = int(np.argmin(np.abs(depth - row.depth_m)))
            self.assertEqual(labels[int(anchored[time_index, depth_index])], row.hybrid_compartment)

    def test_pipeline_wires_basin_backed_asv_curtain_process(self) -> None:
        pipeline = (ROOT / "asv_pipeline.nf").read_text()
        workflow = (ROOT / "workflow" / "modules" / "metadata_ecology.nf").read_text()
        launcher = (ROOT / "run_asv_pipeline.sh").read_text()
        config = (ROOT / "examples" / "si.local.yml").read_text()

        self.assertIn("ASV_TIME_DEPTH_CURTAIN(baseMetadata.map { it })", pipeline)
        self.assertIn("process ASV_TIME_DEPTH_CURTAIN", workflow)
        curtain_process = workflow.split("process ASV_TIME_DEPTH_CURTAIN", 1)[1].split(
            "\nprocess PLOT_UPSET", 1
        )[0]
        self.assertIn(
            "def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0",
            curtain_process,
        )
        self.assertIn('echo "aspire cache generation: ${aspireCacheGeneration}"', curtain_process)
        self.assertIn('"\\${CONDA_PREFIX}/bin/python" "${asvTimeDepthCurtainScriptPath}"', workflow)
        self.assertIn("ASV_TIME_DEPTH_CURTAIN", launcher)
        self.assertIn("asv_time_depth_curtain:\n  enabled: true", config)
        curtain_config = config.split("\nasv_time_depth_curtain:\n", 1)[1].split("\nasv_mag_network:\n", 1)[0]
        self.assertIn("source_mode: basin", curtain_config)
        self.assertIn("basin_grid:", curtain_config)
        self.assertIn("basin_cells:", curtain_config)
        self.assertIn("renewal_events:", curtain_config)
        self.assertIn("asv_point_size: 42.0", curtain_config)
        self.assertIn("asv_point_edge_width: 1.2", curtain_config)
        self.assertIn("maximum_depth: 210.0", curtain_config)
        self.assertIn("--asv-point-edge-width", curtain_process)
        self.assertIn('palette: "oxic__gmm0=#660000,', curtain_config)
        self.assertNotIn("palette: *o2_subcompartment_palette", curtain_config)

        palette_specs = re.findall(
            r'"(oxic__gmm0=#[0-9A-Fa-f]{6},[^"\n]+outlier=#FFFFFF)"',
            config,
        )
        self.assertEqual(len(palette_specs), 13)
        self.assertEqual(len(set(palette_specs)), 1)
        self.assertEqual(
            MODULE.parse_palette(palette_specs[0]),
            {
                "oxic__gmm0": "#660000", "oxic__gmm1": "#D40000",
                "oxic__gmm2": "#FF4444", "oxic__gmm3": "#FFB2B2",
                "dysoxic__gmm0": "#003300", "dysoxic__gmm1": "#8CFF8C",
                "suboxic__gmm0": "#26667C", "suboxic__gmm1": "#3897B7",
                "suboxic__gmm2": "#6AB8D1", "suboxic__gmm3": "#A5D4E3",
                "suboxic__gmm4": "#E0F0F5", "anoxic__gmm0": "#330033",
                "anoxic__gmm1": "#890089", "anoxic__gmm2": "#DF00DF",
                "anoxic__gmm3": "#FF36FF", "anoxic__gmm4": "#FF8CFF",
                "outlier": "#FFFFFF",
            },
        )

    def test_invalid_palette_alias_fails_instead_of_rendering_gray(self) -> None:
        with self.assertRaisesRegex(ValueError, "no label=color entries"):
            MODULE.validate_palette(
                MODULE.parse_palette("o2_subcompartment_palette"),
                ["oxic__gmm0"],
                "o2_subcompartment_palette",
            )


if __name__ == "__main__":
    unittest.main()
