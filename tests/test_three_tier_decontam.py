from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path

import pandas as pd


SCRIPT = (
    Path(__file__).parents[1]
    / "processes/three_tier_decontam/pipeline/build_decontam_metadata.py"
)


class ThreeTierMetadataTest(unittest.TestCase):
    def test_custom_output_sample_column_is_preserved(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            counts = root / "counts.tsv"
            metadata = root / "metadata.tsv"
            output = root / "decontam_metadata.tsv"
            counts.write_text("ASV\tS-1\tNegative-96\nASV1\t10\t1\n")
            metadata.write_text(
                "sample_id\tType_Group\tDNA_conc\n"
                "S_1\tBAL\t2.5\n"
                "Negative_96\tControl\t1.0\n"
            )

            subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "--counts",
                    str(counts),
                    "--metadata-in",
                    str(metadata),
                    "--metadata-out",
                    str(output),
                    "--sample-col-in",
                    "sample_id",
                    "--sample-col-out",
                    "sample_id",
                    "--type-col",
                    "Type_Group",
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            generated = pd.read_csv(output, sep="\t")
            self.assertEqual(generated.columns[0], "sample_id")
            self.assertNotIn("Sample", generated.columns)
            self.assertEqual(generated["sample_id"].tolist(), ["S-1", "Negative-96"])
            self.assertEqual(generated["is_negative_control"].tolist(), [False, True])


if __name__ == "__main__":
    unittest.main()
