"""Check declared file collection with separate producer and cloud directories.

This models the Google backend's collection boundary, without a cloud job:
literal File outputs and glob results are copied, but paths inside JSON are not.
The Cromwell tests separately exercise the complete workflow graph.
"""

import glob
import hashlib
from pathlib import Path
import shutil
import unittest

import WDL

import test_wdl


class CollectedOutputs(WDL.StdLib.Base):
    def __init__(self, execution, storage):
        super().__init__("1.0")
        self.execution = execution.resolve()
        self.storage = storage
        self.glob = WDL.StdLib.StaticFunction(
            "glob", [WDL.Type.String()], WDL.Type.Array(WDL.Type.File()), self.collect_glob)

    def _devirtualize_filename(self, filename):
        return str(self.execution / filename)

    def collect(self, path, relative):
        target = self.storage / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        return "gs://test-output/" + relative

    def collect_glob(self, pattern):
        directory = "glob-" + hashlib.md5(pattern.value.encode()).hexdigest()
        paths = sorted(glob.glob(str(self.execution / pattern.value)))
        return WDL.Value.Array(WDL.Type.File(), [
            WDL.Value.File(self.collect(path, directory + "/" + Path(path).name)) for path in paths])

    def evaluate(self, task):
        env = WDL.Env.Bindings()
        for output in task.outputs:
            value = output.expr.eval(env, self).coerce(output.type)
            # The cloud backend can collect a literal File output before it
            # evaluates output manifests. It cannot discover manifest contents.
            if isinstance(output.type, WDL.Type.File) and isinstance(output.expr, WDL.Expr.String):
                relative = value.value
                value = WDL.Value.File(self.collect(self.execution / relative, relative))
            def map_output(file):
                if file.value.startswith("gs://"):
                    return file.value
                return "gs://test-output/" + str(Path(file.value).relative_to(self.execution))
            value = WDL.Value.rewrite_paths(value, map_output)
            env = env.bind(output.name, value)
        # Producer storage is gone before the next task can read its inputs.
        shutil.rmtree(self.execution)
        return env


class CloudOutputTests(unittest.TestCase):
    def setUp(self):
        self.commands = test_wdl.WdlLocalizationTests()
        self.commands.setUp()
        self.addCleanup(self.commands.doCleanups)

    def collected(self, task_name, inputs):
        commands = self.commands
        work = commands.run_task(commands.render(task_name, inputs), "producer")
        task = next(t for t in commands.doc.tasks if t.name == task_name)
        storage = commands.work / "cloud"
        outputs = CollectedOutputs(work, storage).evaluate(task)
        paths = []
        for binding in outputs:
            def readable(file):
                self.assertTrue(file.value.startswith("gs://test-output/"))
                path = storage / file.value.removeprefix("gs://test-output/")
                self.assertTrue(path.is_file(), "File output was not collected: " + file.value)
                self.assertGreater(path.stat().st_size, 0)
                paths.append(path)
                return file.value
            WDL.Value.rewrite_paths(binding.value, readable)
        return outputs, paths

    def test_split_genotypes_survive_producer_removal(self):
        outputs, paths = self.collected("SplitGenotypes", {})
        self.assertEqual(sum(path.suffix == ".pgen" for path in paths), 1)
        self.assertEqual(sum(path.suffix == ".pvar" for path in paths), 1)
        self.assertEqual(sum(path.suffix == ".psam" for path in paths), 1)

    def test_batch_prediction_lists_and_loco_files_survive_producer_removal(self):
        commands = self.commands
        _, paths = self.collected("CollectPredictions", {
            "prediction_lists": [str(commands.localized / "fit_pred.list")],
            "prediction_files": [str(commands.localized / "fit_1.loco")],
            "prediction_file_counts": [1], "batch_size": 1,
        })
        self.assertIn("batch_0_pred.list", [path.name for path in paths])
        self.assertIn("trait_1.loco", [path.name for path in paths])

    def test_chromosome_file_arrays_keep_numeric_order_and_matched_triplets(self):
        (self.commands.localized / "data.pvar").write_text(
            "#CHROM POS ID REF ALT\n1 1 a A G\n2 2 b A G\n10 3 c A G\nX 4 d A G\nMT 5 e A G\n")
        outputs, _ = self.collected("SplitGenotypes", {})
        for suffix in ("pgen", "pvar", "psam"):
            names = [Path(file.value).name for file in outputs.resolve(suffix + "_files").value]
            self.assertEqual(names, ["chr1." + suffix, "chr2." + suffix, "chr10." + suffix,
                                     "chrX." + suffix, "chrMT." + suffix])

    def test_batch_indices_select_correct_collected_files_after_producer_removal(self):
        commands = self.commands
        traits = [f"p{i}" for i in range(1, 13)]
        sources = []
        rows = []
        for index, trait in enumerate(traits, 1):
            source = commands.localized / f"fit_{index}.loco"
            source.write_text(trait + " prediction\n")
            sources.append(str(source))
            rows.append(trait + " gs://old-model/" + source.name + "\n")
        listing = commands.localized / "fit_pred.list"
        listing.write_text("".join(rows))
        outputs, _ = self.collected("CollectPredictions", {
            "prediction_lists": [str(listing)], "prediction_files": sources,
            "prediction_file_counts": [12], "batch_size": 1,
            "phenotype_columns": ",".join(traits),
        })
        files = outputs.resolve("collected_prediction_files").json
        lists = outputs.resolve("batch_prediction_lists").json
        indices = outputs.resolve("batch_indices").json
        storage = commands.work / "cloud"
        for trait, batch in zip(traits, indices):
            self.assertEqual(batch["phenotype_columns"], trait)
            remote_prediction = files[batch["prediction_file_indices"][0]]
            remote_list = lists[batch["prediction_list_index"]]
            self.assertEqual((storage / remote_prediction.removeprefix("gs://test-output/")).read_text(), trait + " prediction\n")
            self.assertEqual((storage / remote_list.removeprefix("gs://test-output/")).read_text().split()[0], trait)


if __name__ == "__main__":
    unittest.main()
