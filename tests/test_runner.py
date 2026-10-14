import argparse
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class RunnerTests(unittest.TestCase):
    def setUp(self):
        script = ROOT / "scripts/run_regenie.py"
        self.assertTrue(script.is_file(), "The localized REGENIE runner is missing")
        spec = importlib.util.spec_from_file_location("run_regenie", script)
        self.runner = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = self.runner
        spec.loader.exec_module(self.runner)
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)

    def file(self, name, content="data\n"):
        path = self.work / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
        return path

    def args(self, step="step1", *extra):
        return self.runner.parser().parse_args([
            step, "--pgen", str(self.file("a.pgen")),
            "--pvar", str(self.file("b.pvar")),
            "--psam", str(self.file("c.psam")),
            "--phenotypes", str(self.file("phenotypes.txt", "FID IID bmi\n0 s1 2\n")),
            "--phenotype-columns", "bmi", "--trait-type", "quantitative", *extra,
        ])

    def test_unresolved_cloud_uri_reports_localization_error(self):
        for uri in ("gs://bucket/x.pgen", "s3://bucket/x.pgen", "https://host/x.pgen"):
            with self.subTest(uri=uri), self.assertRaisesRegex(ValueError, "localization"):
                self.runner.local_file(uri, "PGEN")

    def test_missing_input_fails_before_computation(self):
        with self.assertRaisesRegex(ValueError, "readable"):
            self.runner.local_file(self.work / "missing.pgen", "PGEN")

    def test_stage_triplet_with_different_names_and_special_characters(self):
        files = [self.file("a space/'$(touch NEVER); " + suffix) for suffix in ("pgen", "pvar", "psam")]
        files[2].write_text("#FID IID SEX\n0 s1 0\n")
        prefix = self.runner.stage_genotypes(*files, self.work / "staged")
        for source, suffix in zip(files, ("pgen", "pvar", "psam")):
            self.assertEqual(Path(str(prefix) + "." + suffix).read_bytes(), source.read_bytes())
        self.assertFalse((self.work / "NEVER").exists())

    def test_iid_only_psam_gets_zero_fid_without_changing_source(self):
        original = "##source=PLINK2\n#IID\tSEX\ns1\t1\ns2\t2\n"
        psam = self.file("samples.psam", original)
        prefix = self.runner.stage_genotypes(self.file("g.pgen"), self.file("v.pvar"), psam, self.work / "staged")
        self.assertEqual(Path(str(prefix) + ".psam").read_text(), "##source=PLINK2\n#FID\tIID\tSEX\n0\ts1\t1\n0\ts2\t2\n")
        self.assertEqual(psam.read_text(), original)

    def test_existing_fid_values_are_preserved(self):
        original = "#FID IID SEX\nfamily1 s1 1\nfamily2 s2 2\n"
        psam = self.file("samples.psam", original)
        prefix = self.runner.stage_genotypes(self.file("g.pgen"), self.file("v.pvar"), psam, self.work / "staged")
        self.assertEqual(Path(str(prefix) + ".psam").read_text(), original)

    def test_minimal_psam_gets_unknown_sex_for_regenie_compatibility(self):
        for index, (content, expected) in enumerate((
            ("#IID\ns1\n", "#FID\tIID\tSEX\n0\ts1\t0\n"),
            ("#FID IID\nfamily1 s1\n", "#FID\tIID\tSEX\nfamily1\ts1\t0\n"),
        )):
            with self.subTest(source=content):
                psam = self.file("samples.psam", content)
                prefix = self.runner.stage_genotypes(self.file("g.pgen"), self.file("v.pvar"), psam, self.work / f"staged{index}")
                self.assertEqual(Path(str(prefix) + ".psam").read_text(), expected)
                self.assertEqual(psam.read_text(), content)

    def test_malformed_psam_is_rejected(self):
        for index, content in enumerate(("#IID SEX\ns1\n", "IID SEX\ns1 2\n", "#FID SEX IID\n0 2 s1\n")):
            with self.subTest(content=content), self.assertRaisesRegex(ValueError, "PSAM"):
                self.runner.stage_genotypes(self.file("g.pgen"), self.file("v.pvar"), self.file("samples.psam", content), self.work / f"staged{index}")

    def test_duplicate_sample_ids_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            self.runner.stage_genotypes(self.file("g.pgen"), self.file("v.pvar"), self.file("samples.psam", "#IID\ns1\ns1\n"), self.work / "staged")

    def test_prediction_paths_are_rebound_to_explicit_localized_files(self):
        pred = self.file("old.list", "bmi gs://old-bucket/fit_1.loco\nheight /old/task/fit_2.loco\n")
        files = [self.file("new local/fit_2.loco"), self.file("new local/fit_1.loco")]
        output = self.runner.rebind_predictions(pred, files, self.work / "predictions")
        rows = [line.split(maxsplit=1) for line in output.read_text().splitlines()]
        self.assertEqual([row[0] for row in rows], ["bmi", "height"])
        self.assertEqual([Path(row[1]).read_text() for row in rows], ["data\n", "data\n"])
        self.assertNotIn("gs://", output.read_text())
        self.assertNotIn("/old/task", output.read_text())

    def test_prediction_path_list_does_not_localize_files_by_itself(self):
        pred = self.file("old.list", "bmi gs://bucket/fit_1.loco\n")
        with self.assertRaisesRegex(ValueError, "explicit|localized"):
            self.runner.rebind_predictions(pred, [], self.work / "predictions")

    def test_duplicate_prediction_basenames_are_rejected(self):
        pred = self.file("old.list", "bmi /old/fit_1.loco\n")
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            self.runner.rebind_predictions(pred, [self.file("a/fit_1.loco"), self.file("b/fit_1.loco")], self.work / "predictions")

    def test_quantitative_command_has_no_binary_correction(self):
        args = self.args("step2", "--prediction-list", "pred.list", "--prediction-files", "files.txt")
        command = self.runner.regenie_command(args, "inputs/genotypes", "inputs/pred.list")
        self.assertIn("--qt", command)
        self.assertNotIn("--firth", command)
        self.assertNotIn("--cc12", command)
        self.assertEqual(command[command.index("--pred") + 1], "inputs/pred.list")

    def test_binary_command_uses_approximate_firth_and_requested_encoding(self):
        args = self.args("step2", "--trait-type", "binary", "--binary-encoding", "12")
        command = self.runner.regenie_command(args, "inputs/genotypes", "pred.list")
        for flag in ("--bt", "--firth", "--approx", "--cc12"):
            self.assertIn(flag, command)

    def test_absent_optional_arguments_are_omitted(self):
        command = self.runner.regenie_command(self.args(), "inputs/genotypes")
        for flag in ("--covarFile", "--keep", "--remove", "--catCovarList"):
            self.assertNotIn(flag, command)
        self.assertNotIn("", command)

    def test_binary_rint_is_rejected(self):
        args = self.args("step1", "--trait-type", "binary", "--rint")
        with self.assertRaisesRegex(ValueError, "quantitative"):
            self.runner.validate_options(args)

    def test_invalid_column_names_and_duplicate_columns_are_rejected(self):
        for value in ("bmi,bmi", "bmi;touch x", "bmi\nheight", "bmi,"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.runner.columns(value)

    def test_unknown_phenotype_column_is_rejected(self):
        args = self.args("step1", "--phenotype-columns", "missing")
        with self.assertRaisesRegex(ValueError, "missing"):
            self.runner.validate_tables(args)

    def test_sample_keep_and_remove_are_combined_before_regenie(self):
        keep = self.file("keep, samples.txt", "0 s1\n0 s2\n0 s3\n")
        remove = self.file("remove, samples.txt", "0 s2\n")
        args = self.args("step1", "--keep", str(keep), "--remove", str(remove))
        self.assertTrue(hasattr(self.runner, "stage_auxiliary"), "Auxiliary input staging is missing")
        self.runner.stage_auxiliary(args, self.work / "staged")
        self.assertIsNone(args.remove)
        self.assertEqual(Path(args.keep).read_text(), "0 s1\n0 s3\n")

    def test_filter_paths_with_commas_are_staged_as_single_files(self):
        variant_ids = self.file("variants, common.txt", "v1\nv2\n")
        args = self.args("step1", "--extract", str(variant_ids))
        self.assertTrue(hasattr(self.runner, "stage_auxiliary"), "Auxiliary input staging is missing")
        self.runner.stage_auxiliary(args, self.work / "staged")
        command = self.runner.regenie_command(args, "inputs/genotypes")
        value = command[command.index("--extract") + 1]
        self.assertNotIn(",", value)
        self.assertEqual(Path(value).read_text(), "v1\nv2\n")

    def test_missing_requested_association_trait_is_an_error(self):
        self.file("results/assoc_bmi.regenie.gz", "compressed result")
        self.assertTrue(hasattr(self.runner, "verify_association_outputs"), "Per-trait output check is missing")
        with self.assertRaisesRegex(ValueError, "height"):
            self.runner.verify_association_outputs("bmi,height", self.work / "results")

    def test_firth_threshold_must_be_strictly_between_zero_and_one(self):
        for value in ("0", "1", "nan", "inf"):
            args = self.args("step2", "--firth-threshold", value, "--prediction-list", "pred.list", "--prediction-files", "files.txt")
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "Firth threshold"):
                self.runner.validate_options(args)


if __name__ == "__main__":
    unittest.main()
