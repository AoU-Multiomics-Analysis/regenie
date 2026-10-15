import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


class PhenotypeBatchTests(unittest.TestCase):
    def setUp(self):
        script = ROOT / "scripts/prepare_batches.py"
        self.assertTrue(script.is_file(), "Phenotype batch planner is missing")
        spec = importlib.util.spec_from_file_location("prepare_batches", script)
        self.batch = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.batch)
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.phenotypes = self.work / "phenotypes ' $(touch INJECTED).txt"
        self.phenotypes.write_text("FID IID a b c d e\n0 sample1 1 2 3 NA 5\n")

    def test_all_traits_are_partitioned_once_in_header_order(self):
        for selection in (None, "", "   "):
            with self.subTest(selection=selection):
                selected, batches = self.batch.phenotype_batches(self.phenotypes, selection, 2)
                self.assertEqual(selected, "a,b,c,d,e")
                self.assertEqual([entry["phenotype_columns"] for entry in batches], ["a,b", "c,d", "e"])
                self.assertEqual([entry["batch_index"] for entry in batches], [0, 1, 2])

    def test_subset_preserves_explicit_selection_order(self):
        selected, batches = self.batch.phenotype_batches(self.phenotypes, "e,a,c", 2)
        self.assertEqual(selected, "e,a,c")
        self.assertEqual([entry["phenotype_columns"] for entry in batches], ["e,a", "c"])

    def test_zero_or_large_batch_size_produces_one_job(self):
        for size in (0, 5, 100):
            with self.subTest(size=size):
                _, batches = self.batch.phenotype_batches(self.phenotypes, None, size)
                self.assertEqual(batches, [{"batch_index": 0, "phenotype_columns": "a,b,c,d,e"}])

    def test_invalid_selection_and_batch_size_fail_before_computation(self):
        for selection, size in (("missing", 2), ("a,a", 2), ("a;touch", 2), (None, -1)):
            with self.subTest(selection=selection, size=size), self.assertRaises(ValueError):
                self.batch.phenotype_batches(self.phenotypes, selection, size)

    def test_bad_header_is_rejected(self):
        for header in ("IID a b", "FID IID a a", "FID IID"):
            self.phenotypes.write_text(header + "\n")
            with self.subTest(header=header), self.assertRaises(ValueError):
                self.batch.phenotype_batches(self.phenotypes, None, 2)

    def test_cli_writes_batch_metadata_and_preserves_source(self):
        before = self.phenotypes.read_bytes()
        result = subprocess.run([sys.executable, str(ROOT / "scripts/prepare_batches.py"),
                                 "--phenotypes", str(self.phenotypes), "--batch-size", "2"],
                                cwd=self.work, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(len(json.loads((self.work / "phenotype_batches.json").read_text())), 3)
        self.assertEqual((self.work / "selected_phenotypes.txt").read_text().strip(), "a,b,c,d,e")
        self.assertEqual(self.phenotypes.read_bytes(), before)
        self.assertFalse((self.work / "INJECTED").exists())

    def test_cli_rejects_unlocalized_phenotype_file(self):
        result = subprocess.run([sys.executable, str(ROOT / "scripts/prepare_batches.py"),
                                 "--phenotypes", "gs://input-bucket/phenotypes.txt", "--batch-size", "2"],
                                cwd=self.work, text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("unresolved cloud URI", result.stderr)


if __name__ == "__main__":
    unittest.main()
