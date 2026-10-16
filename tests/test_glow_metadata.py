"""Preserve valid trait names when prepared files reach the driver."""
import tempfile
import unittest
from pathlib import Path

class MetadataTests(unittest.TestCase):
    def test_batch_trait_names_round_trip_without_missing_value_coercion(self):
        from scripts.glow.spark_main import read_phenotype_batches
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "batches.tsv"
            names = ["NA", "NULL", "NaN", "disease"]
            path.write_text("batch_index\ttrait_index\ttrait_name\n" +
                            "".join(f"0\t{i}\t{name}\n" for i, name in enumerate(names)))
            result = read_phenotype_batches(path)
            self.assertEqual(result.trait_name.tolist(), names)
