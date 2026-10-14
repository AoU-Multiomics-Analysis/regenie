"""Exercise the production command templates with Cromwell's local File mapper."""

import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(os.environ.get("CROMWELL_JAR"), "Set CROMWELL_JAR to run the engine localization test")
class CromwellLocalizationTests(unittest.TestCase):
    def test_production_workflow_localizes_every_file_before_command_quoting(self):
        with tempfile.TemporaryDirectory(prefix="regenie-cromwell-") as directory:
            work = Path(directory)
            data = work / "source"
            data.mkdir()
            for name, text in {
                "data.pgen": "synthetic", "data.pvar": "#CHROM POS ID REF ALT\n1 1 v1 A G\n2 2 v2 T C\n",
                "data.psam": "#IID\ns1\n", "phenotypes.txt": "FID IID bmi\n0 s1 2\n",
                "covariates.txt": "FID IID age\n0 s1 40\n", "keep.txt": "0 s1\n",
                "remove.txt": "0 s2\n", "variants.txt": "v1\nv2\n",
            }.items():
                (data / name).write_text(text)
            # This probe checks the actual command arguments and supplies outputs.
            # The image smoke test separately exercises real REGENIE computation.
            probe = work / "check_paths.py"
            probe.write_text('''import json, sys
from pathlib import Path
args = sys.argv[1:]
for flag in ("--pgen", "--pvar", "--psam", "--phenotypes", "--covariates", "--keep", "--remove", "--extract", "--prediction-list"):
    if flag in args:
        value = args[args.index(flag) + 1]
        assert "/inputs/" in value, "File was not mapped to a task input path: " + flag + " " + value
        assert Path(value).is_file(), "Localized file is missing: " + value
if "--prediction-files" in args:
    paths = Path(args[args.index("--prediction-files") + 1]).read_text().splitlines()
    assert paths and all("/inputs/" in value and Path(value).is_file() for value in paths)
result = Path("results")
result.mkdir()
if args[0] == "step1":
    (result / "fit.log").write_text("GENOTYPE_LOCALIZATION_OK\\n")
    (result / "fit_1.loco").write_text("synthetic predictions\\n")
    (result / "fit_pred.list").write_text("bmi results/fit_1.loco\\n")
else:
    (result / "assoc.log").write_text("GENOTYPE_LOCALIZATION_OK\\n")
    (result / "assoc_bmi.regenie.gz").write_bytes(b"synthetic association output")
print("COMMAND_FILE_LOCALIZATION_OK", args[0])
''')
            source = (ROOT / "workflows/regenie.wdl").read_text()
            # Remove Docker only for this local engine test. Keep the production
            # workflow graph, task inputs, declarations, and command expressions.
            source = re.sub(r"(?m)^\s+docker: docker_image\n", "\n", source)
            source = source.replace("/opt/regenie/scripts/run_regenie.py", shlex.quote(str(probe)))
            # Run the production splitter and output manifest reader. Substitute
            # only PLINK computation, so Cromwell must collect the new File
            # outputs and localize them for the second scatter.
            source = source.replace("/opt/regenie/scripts/split_genotypes.py", shlex.quote(str(ROOT / "scripts/split_genotypes.py")))
            bin_dir = work / "bin"
            bin_dir.mkdir()
            plink = bin_dir / "plink2"
            plink.write_bytes((ROOT / "tests/fake_plink2.py").read_bytes())
            plink.chmod(0o755)
            source = source.replace("set -euo pipefail", "set -euo pipefail\n    export PATH=" + shlex.quote(str(bin_dir)) + ":\"$PATH\"")
            workflow = work / "localization.wdl"
            workflow.write_text(source)
            genotypes = {suffix: str(data / f"data.{suffix}") for suffix in ("pgen", "pvar", "psam")}
            inputs = {
                "Regenie.step1_genotypes": genotypes,
                "Regenie.association_genotypes": [genotypes],
                "Regenie.step1_variant_ids": str(data / "variants.txt"),
                "Regenie.association_variant_ids": str(data / "variants.txt"),
                "Regenie.phenotypes": str(data / "phenotypes.txt"),
                "Regenie.covariates": str(data / "covariates.txt"),
                "Regenie.keep_samples": str(data / "keep.txt"),
                "Regenie.remove_samples": str(data / "remove.txt"),
                "Regenie.docker_image": "unused-local-test",
                "Regenie.step1_cpu": 1, "Regenie.step2_cpu": 1,
                "Regenie.step1_memory_gb": 1, "Regenie.step2_memory_gb": 1,
                "Regenie.split_cpu": 1, "Regenie.split_memory_gb": 1,
                "Regenie.max_retries": 0,
            }
            input_file = work / "inputs.json"
            input_file.write_text(json.dumps(inputs))
            config = work / "cromwell.conf"
            config.write_text('include required(classpath("application"))\n'
                              'backend.default = "Local"\ncall-caching.enabled = false\n')
            metadata = work / "metadata.json"
            result = subprocess.run([
                "java", "-Dconfig.file=" + str(config), "-jar", os.environ["CROMWELL_JAR"],
                "run", str(workflow), "-i", str(input_file), "-m", str(metadata),
            ], cwd=work, text=True, capture_output=True, timeout=300)
            self.assertEqual(result.returncode, 0, (result.stdout + result.stderr)[-16000:])
            observed = json.loads(metadata.read_text())
            self.assertEqual(observed["status"], "Succeeded")
            outputs = observed["outputs"]
            self.assertEqual(len(outputs["Regenie.associations_by_shard"]), 2)
            chromosomes = outputs["Regenie.association_chromosomes"]
            self.assertEqual([entry["chromosome"] for entry in chromosomes], ["1", "2"])
            self.assertEqual(len(outputs["Regenie.split_manifests"]), 1)
            for entry in chromosomes:
                for suffix in ("pgen", "pvar", "psam"):
                    self.assertTrue(Path(entry[suffix]).is_file())
                self.assertIn(entry["chromosome"] + " ", Path(entry["pvar"]).read_text())
            for log in [outputs["Regenie.step1_log"], *outputs["Regenie.step2_logs"]]:
                self.assertEqual(Path(log).read_text(), "GENOTYPE_LOCALIZATION_OK\n")


if __name__ == "__main__":
    unittest.main()
