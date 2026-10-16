import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from scripts.glow.convert import convert_pgen

class ConversionTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.pgen=self.root/"a.pgen"; self.pgen.write_bytes(b"synthetic")
        self.pvar=self.root/"a.pvar"; self.pvar.write_text("#CHROM POS ID REF ALT\n1 2 v1 A G\n2 3 v2 C T\n")
        self.psam=self.root/"a.psam"; self.psam.write_text("#IID\n001\n002\n003\n")
        self.mapping=self.root/"map.tsv"; self.mapping.write_text("sample_id\tFID\tIID\ns_two\t0\t002\ns_one\t0\t001\n")

    def convert(self, **kw):
        args=dict(pgen=self.pgen,pvar=self.pvar,psam=self.psam,sample_map=self.mapping,
                  variant_ids=None,source_index=0,role="association",out_dir=self.root/"out",threads=1,memory_mb=512)
        args.update(kw)
        return convert_pgen(**args)

    def fake_export(self, args, **kwargs):
        prefix=Path(args[args.index("--out")+1])
        source=Path(args[args.index("--pfile")+1])
        if "--make-pgen" in args:
            import shutil
            order=Path(args[args.index("--indiv-sort")+2]).read_text().splitlines()
            prefix.with_suffix(".pgen").write_bytes(b"ordered")
            shutil.copyfile(source.with_suffix(".pvar"),prefix.with_suffix(".pvar"))
            prefix.with_suffix(".psam").write_text("#FID IID\n"+"\n".join(order)+"\n")
        else:
            order=source.with_suffix(".psam").read_text().splitlines()[1:]
            prefix.with_suffix(".bgen").write_bytes(b"bgen")
            prefix.with_suffix(".sample").write_text("ID_1 ID_2 missing\n0 0 0\n"+"".join("0 "+row.split()[1]+" 0\n" for row in order))
        self.args=args

    def test_export_uses_dosages_and_order_without_modifying_source(self):
        before=self.psam.read_bytes()
        with patch("scripts.glow.convert.subprocess.run",side_effect=self.fake_export):
            out=self.convert()
        self.assertEqual(self.psam.read_bytes(),before)
        self.assertTrue(out.bgen.is_file())
        self.assertIn("ref-first",self.args)
        self.assertIn("bits=16",self.args)
        self.assertIn("id-paste=iid",self.args)
        self.assertEqual(json.loads(out.report.read_text())["source_index"],0)

    def test_bad_metadata_fails_before_plink(self):
        for row in ("X 2 v1 A G","1 2 v1 A G,T","1 2 v1 A ."):
            with self.subTest(row=row):
                self.pvar.write_text("#CHROM POS ID REF ALT\n"+row+"\n")
                with self.assertRaises(ValueError):
                    self.convert()

    def test_missing_retained_sample_and_duplicate_variant_fail(self):
        self.psam.write_text("#IID\n001\n")
        with self.assertRaisesRegex(ValueError,"sample"):
            self.convert()
        self.psam.write_text("#IID\n001\n002\n")
        self.pvar.write_text("#CHROM POS ID REF ALT\n1 2 v1 A G\n2 3 v1 C T\n")
        with self.assertRaisesRegex(ValueError,"duplicate"):
            self.convert()

    def test_sources_with_identical_basenames_keep_separate_outputs(self):
        import shutil
        other = self.root / "second"
        other.mkdir()
        for path in (self.pgen, self.pvar, self.psam):
            shutil.copyfile(path, other / path.name)
        with patch("scripts.glow.convert.subprocess.run", side_effect=self.fake_export):
            first = self.convert(out_dir=self.root / "src0", source_index=0)
            second = self.convert(pgen=other / self.pgen.name, pvar=other / self.pvar.name,
                                  psam=other / self.psam.name, out_dir=self.root / "src1", source_index=1)
        self.assertEqual(first.bgen.name, second.bgen.name)
        self.assertNotEqual(first.bgen, second.bgen)
        self.assertEqual(json.loads(first.report.read_text())["source_index"], 0)
        self.assertEqual(json.loads(second.report.read_text())["source_index"], 1)

    def test_sort_happens_before_bgen_export(self):
        calls = []
        def record(args, **kwargs):
            calls.append(args)
            return self.fake_export(args, **kwargs)
        with patch("scripts.glow.convert.subprocess.run", side_effect=record):
            self.convert()
        self.assertEqual(len(calls), 2)
        self.assertIn("--make-pgen", calls[0])
        self.assertIn("--indiv-sort", calls[0])
        self.assertNotIn("--indiv-sort", calls[1])
        self.assertIn("--export", calls[1])
        self.assertEqual(calls[0][calls[0].index("--out") + 1], calls[1][calls[1].index("--pfile") + 1])

    def test_plink_failure_prints_its_diagnostic_log(self):
        import subprocess
        def fail(args, **kwargs):
            kwargs["stdout"].write("PLINK_DIAGNOSTIC: unsupported export option\n")
            kwargs["stdout"].flush()
            raise subprocess.CalledProcessError(8, args)
        with patch("scripts.glow.convert.subprocess.run", side_effect=fail):
            with self.assertLogs(level="ERROR") as logs:
                with self.assertRaises(subprocess.CalledProcessError):
                    self.convert()
            self.assertIn("PLINK_DIAGNOSTIC", "\n".join(logs.output))
