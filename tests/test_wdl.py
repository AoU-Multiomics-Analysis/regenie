import os
import shlex
import subprocess
import tempfile
import unittest
from pathlib import Path

import WDL

ROOT = Path(__file__).resolve().parents[1]


class LocalStdLib(WDL.StdLib.Base):
    def _devirtualize_filename(self, filename):
        return filename


class WdlLocalizationTests(unittest.TestCase):
    def setUp(self):
        wdl = ROOT / "workflows/regenie.wdl"
        self.assertTrue(wdl.is_file(), "REGENIE WDL is missing")
        self.doc = WDL.load(str(wdl))
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.localized = self.work / "localized ' $(touch INJECTED)"
        self.localized.mkdir()
        for name, text in {
            "data.pgen": "pgen", "data.pvar": "#CHROM POS ID REF ALT\n1 1 v1 A G\n",
            "data.psam": "#FID IID\n0 s1\n", "phenotypes.txt": "FID IID bmi\n0 s1 2\n",
            "variants.txt": "v1\n", "fit_1.loco": "prediction\n",
            "fit_pred.list": "bmi gs://old-bucket/fit_1.loco\n",
        }.items():
            (self.localized / name).write_text(text)
        # REGENIE itself is external. This executable records the runner's real
        # arguments and supplies its documented output files. Real REGENIE is
        # tested in the separate GitHub Actions image smoke test.
        bin_dir = self.work / "bin"
        bin_dir.mkdir()
        fake = bin_dir / "regenie"
        fake.write_text('''#!/usr/bin/env python3
import json, sys
from pathlib import Path
args = sys.argv[1:]
Path("argv.json").write_text(json.dumps(args))
out = args[args.index("--out") + 1]
Path(out + ".log").write_text("test executable\\n")
if args[args.index("--step") + 1] == "1":
    Path(out + "_1.loco").write_text("prediction\\n")
    Path(out + "_pred.list").write_text("bmi " + out + "_1.loco\\n")
else:
    pred = Path(args[args.index("--pred") + 1]).read_text().split(maxsplit=1)[1].strip()
    assert Path(pred).read_text() == "prediction\\n"
    Path(out + "_bmi.regenie.gz").write_text("association\\n")
''')
        fake.chmod(0o755)
        plink = bin_dir / "plink2"
        plink.write_text('''#!/usr/bin/env python3
import json, shutil, sys
from pathlib import Path
args = sys.argv[1:]
Path("plink_argv.json").write_text(json.dumps(args))
prefix = args[args.index("--pfile") + 1]
out = args[args.index("--out") + 1]
assert Path(prefix + ".pvar").is_file()
assert Path(args[args.index("--extract") + 1]).read_text() == "v1\\n"
assert Path(args[args.index("--keep") + 1]).read_text() == "0 s1\\n"
assert "--remove" not in args
for suffix in ("pgen", "pvar", "psam"):
    shutil.copyfile(prefix + "." + suffix, out + "." + suffix)
Path(out + ".log").write_text("PLINK test executable\\n")
''')
        plink.chmod(0o755)
        self.env = dict(os.environ, PATH=str(bin_dir) + os.pathsep + os.environ["PATH"])

    def render(self, task_name, extra=None):
        task = next(t for t in self.doc.tasks if t.name == task_name)
        inputs = {
            "genotypes": {"pgen": str(self.localized / "data.pgen"), "pvar": str(self.localized / "data.pvar"), "psam": str(self.localized / "data.psam")},
            "docker_image": "test:image",
        }
        if task_name != "SubsetPgen":
            inputs.update(phenotypes=str(self.localized / "phenotypes.txt"), phenotype_columns="bmi",
                          trait_type="quantitative")
        inputs.update(extra or {})
        # Start with cloud File values, including files inside the struct.
        # Simulate the engine's localization before declarations and commands.
        def cloud_values(value):
            if isinstance(value, dict):
                return {key: cloud_values(item) for key, item in value.items()}
            if isinstance(value, list):
                return [cloud_values(item) for item in value]
            if isinstance(value, str) and value.startswith(str(self.localized) + "/"):
                return "gs://input-bucket/" + Path(value).name
            return value

        env = WDL.values_from_json(cloud_values(inputs), task.available_inputs)
        env = WDL.Value.rewrite_env_paths(env, lambda file: str(self.localized / file.value.rsplit("/", 1)[-1]))
        stdlib = LocalStdLib("1.0")
        for decl in [*task.inputs, *task.postinputs]:
            if not env.has_binding(decl.name):
                value = decl.expr.eval(env, stdlib) if decl.expr else WDL.Value.Null()
                env = env.bind(decl.name, value)
        return task.command.eval(env, stdlib).value

    def run_task(self, command, name):
        work = self.work / name
        work.mkdir()
        # Use the production runner while retaining the WDL command's arguments.
        command = command.replace("/opt/regenie/scripts/run_regenie.py", shlex.quote(str(ROOT / "scripts/run_regenie.py")))
        command = command.replace("/opt/regenie/scripts/subset_genotypes.py", shlex.quote(str(ROOT / "scripts/subset_genotypes.py")))
        result = subprocess.run(["bash", "-c", command], cwd=work, env=self.env, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertFalse((work / "INJECTED").exists())
        return work

    def test_subset_localizes_struct_and_optional_filters_before_plink(self):
        self.assertIn("SubsetPgen", [task.name for task in self.doc.tasks], "PLINK subset task is missing")
        (self.localized / "data.psam").write_text("#IID\ns1\n")
        (self.localized / "keep.txt").write_text("0 s1\n0 s2\n")
        (self.localized / "remove.txt").write_text("0 s2\n")
        command = self.render("SubsetPgen", {
            "variant_ids": str(self.localized / "variants.txt"),
            "keep_samples": str(self.localized / "keep.txt"),
            "remove_samples": str(self.localized / "remove.txt"),
        })
        work = self.run_task(command, "subset")
        self.assertEqual((work / "results/subset.pgen").read_text(), "pgen")
        self.assertEqual((work / "results/subset.psam").read_text(), "#FID\tIID\tSEX\n0\ts1\t0\n")
        self.assertEqual((self.localized / "data.psam").read_text(), "#IID\ns1\n")

    def test_subset_rejects_unlocalized_cloud_files_before_plink(self):
        script = ROOT / "scripts/subset_genotypes.py"
        self.assertTrue(script.is_file(), "PLINK subset runner is missing")
        result = subprocess.run(["python3", str(script), "--pgen", "gs://bucket/data.pgen",
                                 "--pvar", str(self.localized / "data.pvar"),
                                 "--psam", str(self.localized / "data.psam")],
                                cwd=self.work, env=self.env, text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("localization error", result.stderr)
        self.assertFalse((self.work / "plink_argv.json").exists())

    def test_workflow_routes_subset_file_structs_to_regenie(self):
        genotype_type = next(t for t in self.doc.tasks if t.name == "FitModel").available_inputs.resolve("genotypes").type
        raw = WDL.Value.from_json(genotype_type, {suffix: f"gs://raw/data.{suffix}" for suffix in ("pgen", "pvar", "psam")})
        subset = WDL.Value.from_json(genotype_type, {suffix: f"gs://subset/data.{suffix}" for suffix in ("pgen", "pvar", "psam")})
        stdlib = LocalStdLib("1.0")
        fit = next(n for n in self.doc.workflow.body if isinstance(n, WDL.Tree.Call) and n.name == "FitModel")
        env = WDL.Env.Bindings().bind("step1_genotypes", raw).bind("SubsetStep1.subset_genotypes", subset)
        result = fit.inputs["genotypes"].eval(env, stdlib)
        self.assertEqual(result.json, {suffix: f"gs://subset/data.{suffix}" for suffix in ("pgen", "pvar", "psam")})
        self.assertTrue(all(isinstance(value, WDL.Value.File) for value in result.value.values()))
        scatter = next(n for n in self.doc.workflow.body if isinstance(n, WDL.Tree.Scatter))
        condition = next(n for n in scatter.body if isinstance(n, WDL.Tree.Conditional))
        selected = next(n for n in scatter.body if isinstance(n, WDL.Tree.Decl))
        test = next(n for n in scatter.body if isinstance(n, WDL.Tree.Call))
        for filters, expected_subset in (({}, False), ({"association_variant_ids": "variants.txt"}, True),
                                         ({"keep_samples": "keep.txt"}, True), ({"remove_samples": "remove.txt"}, True)):
            with self.subTest(filters=filters):
                env = WDL.Env.Bindings().bind("genotype_set", raw)
                for name in ("association_variant_ids", "keep_samples", "remove_samples"):
                    env = env.bind(name, WDL.Value.File(filters[name]) if name in filters else WDL.Value.Null())
                self.assertEqual(condition.expr.eval(env, stdlib).value, expected_subset)
                env = env.bind("SubsetAssociation.subset_genotypes", subset if expected_subset else WDL.Value.Null())
                env = env.bind("selected_genotypes", selected.expr.eval(env, stdlib))
                result = test.inputs["genotypes"].eval(env, stdlib)
                self.assertEqual(result.json, (subset if expected_subset else raw).json)
                self.assertTrue(all(isinstance(value, WDL.Value.File) for value in result.value.values()))

    def test_step1_renders_localized_files_without_shell_expansion(self):
        command = self.render("FitModel", {"variant_ids": str(self.localized / "variants.txt")})
        work = self.run_task(command, "step1")
        self.assertEqual((work / "inputs/genotypes.pgen").read_text(), "pgen")
        self.assertTrue((work / "results/fit_pred.list").is_file())

    def test_step2_lists_localized_files_at_command_time_and_rebinds_cloud_metadata(self):
        command = self.render("TestVariants", {
            "prediction_list": str(self.localized / "fit_pred.list"),
            "prediction_files": [str(self.localized / "fit_1.loco")],
        })
        work = self.run_task(command, "step2")
        self.assertEqual((work / "prediction_files.txt").read_text().strip(), str(self.localized / "fit_1.loco"))
        self.assertNotIn("gs://", (work / "inputs/predictions/predictions.list").read_text())

    def test_both_tasks_normalize_localized_psam_without_fid(self):
        (self.localized / "data.psam").write_text("#IID SEX\ns1 2\n")
        inputs = {
            "FitModel": {"variant_ids": str(self.localized / "variants.txt")},
            "TestVariants": {"prediction_list": str(self.localized / "fit_pred.list"), "prediction_files": [str(self.localized / "fit_1.loco")]},
        }
        for task_name, extra in inputs.items():
            with self.subTest(task=task_name):
                work = self.run_task(self.render(task_name, extra), task_name)
                self.assertEqual((work / "inputs/genotypes.psam").read_text(), "#FID\tIID\tSEX\n0\ts1\t2\n")
        self.assertEqual((self.localized / "data.psam").read_text(), "#IID SEX\ns1 2\n")

    def test_present_optional_files_and_binary_settings_are_rendered_safely(self):
        (self.localized / "covariates.txt").write_text("FID IID age batch\n0 s1 40 A\n")
        (self.localized / "keep.txt").write_text("0 s1\n")
        (self.localized / "remove.txt").write_text("0 s2\n")
        command = self.render("TestVariants", {
            "prediction_list": str(self.localized / "fit_pred.list"),
            "prediction_files": [str(self.localized / "fit_1.loco")],
            "covariates": str(self.localized / "covariates.txt"),
            "covariate_columns": "age,batch", "categorical_covariates": "batch",
            "keep_samples": str(self.localized / "keep.txt"), "remove_samples": str(self.localized / "remove.txt"),
            "variant_ids": str(self.localized / "variants.txt"), "trait_type": "binary", "binary_encoding": "12",
        })
        work = self.run_task(command, "optional")
        import json
        args = json.loads((work / "argv.json").read_text())
        for flag in ("--bt", "--cc12", "--firth", "--approx", "--covarFile", "--keep", "--extract"):
            self.assertIn(flag, args)
        self.assertNotIn("--remove", args)
        self.assertEqual((work / args[args.index("--covarFile") + 1]).read_text(), "FID IID age batch\n0 s1 40 A\n")
        self.assertEqual((work / args[args.index("--keep") + 1]).read_text(), "0 s1\n")


if __name__ == "__main__":
    unittest.main()
