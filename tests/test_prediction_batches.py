import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


class PredictionBatchTests(unittest.TestCase):
    def setUp(self):
        script = ROOT / "scripts/collect_predictions.py"
        self.assertTrue(script.is_file(), "Prediction collector is missing")
        spec = importlib.util.spec_from_file_location("collect_predictions", script)
        self.collector = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.collector)
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)

    def model(self, index, traits):
        source = self.work / f"model {index} ' input"
        source.mkdir()
        paths = []
        rows = []
        for i, trait in enumerate(traits, 1):
            path = source / f"fit_{i}.loco"
            path.write_text(trait + " prediction\n")
            paths.append(path)
            rows.append(trait + " gs://old-bucket/fit_" + str(i) + ".loco\n")
        listing = source / "fit_pred.list"
        listing.write_text("".join(rows))
        return listing, paths

    def test_shared_fit_is_partitioned_without_renaming_traits_or_losing_predictions(self):
        listing, files = self.model(0, ["a", "b", "c"])
        result = self.collector.collect_predictions([listing], files[::-1], [3], "a,b,c", 2, self.work / "results")
        self.assertEqual([row["phenotype_columns"] for row in result], ["a,b", "c"])
        self.assertEqual([len(row["prediction_files"]) for row in result], [2, 1])
        for row in result:
            for trait, filename in zip(row["phenotype_columns"].split(","), row["prediction_files"]):
                self.assertEqual(Path(filename).read_text(), trait + " prediction\n")
            self.assertNotIn("gs://", Path(row["prediction_list"]).read_text())

    def test_separate_fits_with_identical_basenames_keep_correct_trait_predictions(self):
        first, first_files = self.model(0, ["a", "b"])
        second, second_files = self.model(1, ["c"])
        result = self.collector.collect_predictions([first, second], first_files + second_files,
                                                    [2, 1], "c,a,b", 1, self.work / "results")
        self.assertEqual([row["phenotype_columns"] for row in result], ["c", "a", "b"])
        names = []
        for row in result:
            names.extend(Path(value).name for value in row["prediction_files"])
            self.assertEqual(Path(row["prediction_files"][0]).read_text(), row["phenotype_columns"] + " prediction\n")
        self.assertEqual(len(set(names)), 3)
        self.assertEqual(len((self.work / "results/all_pred.list").read_text().splitlines()), 3)
        self.assertEqual(json.loads((self.work / "results/prediction_batches.json").read_text()), result)

    def test_bad_counts_and_missing_or_extra_traits_fail(self):
        listing, files = self.model(0, ["a", "b"])
        for counts, selected in (([1], "a,b"), ([2, 1], "a,b"), ([2], "a"), ([2], "a,b,c")):
            with self.subTest(counts=counts, selected=selected), self.assertRaises(ValueError):
                self.collector.collect_predictions([listing], files, counts, selected, 1, self.work / "results")

    def test_duplicate_traits_across_models_fail(self):
        first, first_files = self.model(0, ["a"])
        second, second_files = self.model(1, ["a"])
        with self.assertRaisesRegex(ValueError, "Duplicate.*phenotype"):
            self.collector.collect_predictions([first, second], first_files + second_files, [1, 1], "a", 1, self.work / "results")

    def test_unlocalized_files_are_rejected(self):
        listing, files = self.model(0, ["a"])
        for lists, predictions in ((["gs://bucket/fit_pred.list"], files), ([listing], ["gs://bucket/fit_1.loco"])):
            with self.subTest(lists=lists), self.assertRaisesRegex(ValueError, "unresolved cloud URI"):
                self.collector.collect_predictions(lists, predictions, [1], "a", 1, self.work / "results")


if __name__ == "__main__":
    unittest.main()
