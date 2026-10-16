"""Cohort tests detect ID changes, silent sample loss, and invalid designs."""
import json
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

from scripts.glow.contracts import ScientificOptions, require_local_file
from scripts.glow.prepare import prepare_cohort, sample_token


class CohortTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.psam = self.write("data.psam", "#IID\n001\n002\n003\n004\n")
        self.pheno = self.write("pheno.tsv", "FID IID a b\n0 001 1 4\n0 002 2 3\n0 003 3 2\n0 004 4 1\n")
        self.options = ScientificOptions(trait_type="quantitative", sample_block_count=2)

    def write(self, name, text):
        p = self.root / name
        p.write_text(text)
        return p

    def run_cohort(self, **kwargs):
        values = dict(psam=self.psam, phenotypes=self.pheno, covariates=None, keep=None,
                      remove=None, options=self.options, out_dir=self.root / "out")
        values.update(kwargs)
        return prepare_cohort(**values)

    def test_tokens_preserve_pairs_and_leading_zeroes(self):
        self.assertNotEqual(sample_token("a_b", "c"), sample_token("a", "b_c"))
        self.assertNotEqual(sample_token("0", "001"), sample_token("0", "1"))

    def test_missing_fid_preserves_ids_and_source_files(self):
        before = self.psam.read_bytes()
        out = self.run_cohort()
        import pandas as pd
        table = pd.read_csv(out.sample_map, sep="\t", dtype=str)
        self.assertEqual(set(table.IID), {"001", "002", "003", "004"})
        self.assertEqual(set(table.FID), {"0"})
        self.assertEqual(self.psam.read_bytes(), before)
        self.assertEqual(len(json.loads(out.report.read_text())["traits"]), 2)

    def test_missing_values_require_explicit_policy(self):
        self.pheno.write_text("FID IID a b\n0 001 1 NA\n0 002 2 3\n0 003 3 2\n0 004 4 1\n")
        with self.assertRaisesRegex(ValueError, "missing"):
            self.run_cohort()
        out = self.run_cohort(options=replace(self.options, missing_data_policy="complete_case", batch_size=1))
        import pandas as pd
        table = pd.read_csv(out.phenotypes, sep="\t")
        self.assertEqual(len(table), 3)
        self.assertEqual(len(pd.read_csv(out.batches, sep="\t")), 2)
        self.assertEqual(json.loads(out.report.read_text())["removed_incomplete"], 1)

    def test_duplicate_pairs_fail(self):
        self.psam.write_text("#IID\n001\n001\n")
        with self.assertRaisesRegex(ValueError, "duplicate"):
            self.run_cohort()

    def test_remove_wins_over_keep(self):
        keep = self.write("keep", "0 001\n0 002\n0 003\n")
        remove = self.write("remove", "0 002\n")
        out = self.run_cohort(keep=keep, remove=remove)
        import pandas as pd
        self.assertEqual(set(pd.read_csv(out.sample_map, sep="\t", dtype=str).IID), {"001", "003"})

    def test_binary_encoding_and_case_limits(self):
        self.pheno.write_text("FID IID a\n0 001 1\n0 002 2\n0 003 1\n0 004 2\n")
        opts = replace(self.options, trait_type="binary", binary_encoding="12", min_cases=2, min_controls=2)
        out = self.run_cohort(options=opts)
        self.assertEqual(json.loads(out.report.read_text())["traits"]["a"]["cases"], 2)
        self.pheno.write_text("FID IID a\n0 001 1\n0 002 2\n0 003 1\n0 004 3\n")
        with self.assertRaisesRegex(ValueError, "binary"):
            self.run_cohort(options=opts)

    def test_collinear_covariates_and_constant_traits_fail(self):
        cov = self.write("cov", "FID IID x y\n0 001 1 2\n0 002 2 4\n0 003 3 6\n0 004 4 8\n")
        with self.assertRaisesRegex(ValueError, "rank"):
            self.run_cohort(covariates=cov)
        self.pheno.write_text("FID IID a\n0 001 1\n0 002 1\n0 003 1\n0 004 1\n")
        with self.assertRaisesRegex(ValueError, "constant"):
            self.run_cohort()

    def test_localization_error_and_invalid_options(self):
        with self.assertRaisesRegex(ValueError, "localization"):
            require_local_file(Path("gs://bucket/data"))
        with self.assertRaises(ValueError):
            replace(self.options, batch_size=-1)
        with self.assertRaises(ValueError):
            replace(self.options, trait_type="other")

    def test_transform_is_deterministic_and_unselected_missing_is_allowed(self):
        opts = replace(self.options, phenotype_columns=("a",), inverse_normal_transform=True)
        self.pheno.write_text("FID IID a b\n0 001 1 NA\n0 002 2 NA\n0 003 3 NA\n0 004 4 NA\n")
        out = self.run_cohort(options=opts)
        first = out.phenotypes.read_bytes()
        self.assertEqual(self.run_cohort(options=opts).phenotypes.read_bytes(), first)
