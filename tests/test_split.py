import importlib.util
import sys
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


class SplitTests(unittest.TestCase):
    def setUp(self):
        script = ROOT / "scripts/split_genotypes.py"
        self.assertTrue(script.is_file(), "Chromosome split script is missing")
        spec = importlib.util.spec_from_file_location("split_genotypes", script)
        self.split = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.split)
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)

    def pvar(self, text):
        path = self.work / "data.pvar"
        path.write_text(text)
        return path

    def test_detects_chromosomes_in_numeric_order_and_preserves_selected_counts(self):
        path = self.pvar("##source=test\n#CHROM POS ID REF ALT\nchr2 20 v2 A G\n1 10 v1 T C\nchrX 30 vx A C\n2 40 v3 A T\n")
        self.assertEqual(self.split.chromosome_counts(path), [("1", 1), ("2", 2), ("X", 1)])
        self.assertEqual(self.split.chromosome_counts(path, {"v2", "vx", "absent"}), [("2", 1), ("X", 1)])

    def test_empty_filter_match_does_not_create_empty_chromosome_jobs(self):
        path = self.pvar("#CHROM POS ID REF ALT\n1 10 v1 A G\n")
        self.assertEqual(self.split.chromosome_counts(path, {"absent"}), [])

    def test_rejects_invalid_pvar_before_plink(self):
        for content in ("CHROM POS ID REF ALT\n1 1 v1 A G\n", "#CHROM POS ID REF ALT\n1 1 v1\n",
                        "#CHROM POS ID REF ALT\nunknown 1 v1 A G\n"):
            with self.subTest(content=content), self.assertRaises(ValueError):
                self.split.chromosome_counts(self.pvar(content))

    def test_cli_rejects_unlocalized_files_before_plink(self):
        for flag in ("pgen", "pvar", "psam", "extract"):
            with self.subTest(flag=flag), tempfile.TemporaryDirectory() as directory:
                work = Path(directory)
                files = {"pgen": "synthetic", "pvar": "#CHROM POS ID REF ALT\n1 1 v1 A G\n",
                         "psam": "#IID\ns1\n", "extract": "v1\n"}
                command = [sys.executable, str(ROOT / "scripts/split_genotypes.py")]
                for suffix, content in files.items():
                    path = work / suffix
                    path.write_text(content)
                    command.extend(["--" + suffix, "gs://input-bucket/" + suffix if flag == suffix else str(path)])
                result = subprocess.run(command, cwd=work, capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("unresolved cloud URI", result.stderr)
                self.assertFalse((work / "chromosome_genotypes.json").exists())

    def test_cli_fails_if_variant_filter_retains_no_variants(self):
        for suffix, content in {"pgen": "synthetic", "pvar": "#CHROM POS ID REF ALT\n1 1 v1 A G\n",
                                "psam": "#IID\ns1\n", "extract": "absent\n"}.items():
            (self.work / suffix).write_text(content)
        command = [sys.executable, str(ROOT / "scripts/split_genotypes.py")]
        for suffix in ("pgen", "pvar", "psam", "extract"):
            command.extend(["--" + suffix, str(self.work / suffix)])
        result = subprocess.run(command, cwd=self.work, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("No association variants remain", result.stderr)
        self.assertFalse((self.work / "chromosome_genotypes.json").exists())

    def test_command_uses_requested_chromosome_and_preserves_genotype_data(self):
        args = self.split.parser().parse_args(["--pgen", "a", "--pvar", "b", "--psam", "c",
                                               "--extract", "inputs/extract.txt", "--threads", "3", "--memory-mb", "2048"])
        self.assertEqual(self.split.plink_command(args, Path("inputs/genotypes"), "2"), [
            "plink2", "--pfile", "inputs/genotypes", "--chr", "2", "--make-pgen", "--out", "chromosomes/chr2",
            "--threads", "3", "--memory", "2048", "--extract", "inputs/extract.txt",
        ])


if __name__ == "__main__":
    unittest.main()
