import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


class SubsetTests(unittest.TestCase):
    def setUp(self):
        script = ROOT / "scripts/subset_genotypes.py"
        self.assertTrue(script.is_file(), "PLINK subset runner is missing")
        spec = importlib.util.spec_from_file_location("subset_genotypes", script)
        self.subset = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.subset)
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)

    def args(self, *extra):
        return self.subset.parser().parse_args([
            "--pgen", "data.pgen", "--pvar", "data.pvar", "--psam", "data.psam", *extra,
        ])

    def test_command_selects_variants_and_samples_and_writes_pgen(self):
        command = self.subset.plink_command(self.args(
            "--extract", "variants.txt", "--keep", "keep.txt", "--remove", "remove.txt",
            "--threads", "3", "--memory-mb", "2048",
        ), Path("inputs/genotypes"))
        self.assertEqual(command, [
            "plink2", "--pfile", "inputs/genotypes", "--make-pgen", "--out", "results/subset",
            "--threads", "3", "--memory", "2048", "--extract", "variants.txt",
            "--keep", "keep.txt", "--remove", "remove.txt",
        ])

    def test_absent_filters_are_omitted(self):
        command = self.subset.plink_command(self.args(), Path("inputs/genotypes"))
        for flag in ("--extract", "--keep", "--remove"):
            self.assertNotIn(flag, command)

    def test_output_check_counts_rows_after_pvar_metadata(self):
        prefix = self.work / "subset"
        prefix.with_suffix(".pgen").write_bytes(b"synthetic")
        prefix.with_suffix(".pvar").write_text("##source=test\n#CHROM POS ID REF ALT\n1 10 v1 A G\n2 20 v2 T C\n")
        prefix.with_suffix(".psam").write_text("#IID SEX\ns1 0\ns2 0\ns3 0\n")
        self.assertEqual(self.subset.check_outputs(prefix), (2, 3))

    def test_output_check_rejects_empty_variant_or_sample_sets(self):
        prefix = self.work / "subset"
        prefix.with_suffix(".pgen").write_bytes(b"synthetic")
        for pvar, psam in (("#CHROM POS ID REF ALT\n", "#IID\ns1\n"),
                           ("#CHROM POS ID REF ALT\n1 1 v1 A G\n", "#IID\n")):
            prefix.with_suffix(".pvar").write_text(pvar)
            prefix.with_suffix(".psam").write_text(psam)
            with self.subTest(pvar=pvar, psam=psam), self.assertRaisesRegex(ValueError, "no variants or samples"):
                self.subset.check_outputs(prefix)


if __name__ == "__main__":
    unittest.main()
