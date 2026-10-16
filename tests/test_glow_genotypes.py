import unittest
from scripts.glow.genotypes import dosage_from_probabilities, observed_mac, impute_values

class DosageTests(unittest.TestCase):
    def test_fractional_dosage_and_missing_values(self):
        self.assertAlmostEqual(dosage_from_probabilities([.1,.6,.3]),1.2)
        self.assertIsNone(dosage_from_probabilities([]))
        self.assertIsNone(dosage_from_probabilities([0,0,0]))
        self.assertEqual(impute_values([.2,1.8,None]),[.2,1.8,1.0])
        self.assertAlmostEqual(observed_mac([.2,1.8,None]),2)
    def test_invalid_probabilities_and_ploidy_fail(self):
        for values in ([.2,.2,.2],[float("nan"),0,1],[-.1,.6,.5],[0,1]):
            with self.subTest(values=values):
                with self.assertRaises(ValueError): dosage_from_probabilities(values)
        with self.assertRaises(ValueError): dosage_from_probabilities([0,1,0],ploidy=1)
        with self.assertRaises(ValueError): dosage_from_probabilities([0,1,0],phased=True)
    def test_empty_and_invalid_dosages_fail(self):
        with self.assertRaises(ValueError): impute_values([None,None])
        with self.assertRaises(ValueError): observed_mac([3])

class BgenReaderTests(unittest.TestCase):
    def test_explicit_sample_file_and_arbitrary_localized_bgen_name(self):
        from pathlib import Path
        from unittest.mock import MagicMock
        from scripts.glow.genotypes import read_bgen
        spark=MagicMock(); reader=spark.read.format.return_value
        reader.option.return_value=reader
        bgen=Path("/work/name ' %/input_without_extension")
        sample=Path("/work/different ' %/samples.sample")
        result=read_bgen(spark,bgen.as_uri(),sample.as_uri())
        self.assertIs(result,reader.load.return_value)
        reader.load.assert_called_once_with(str(bgen))
        calls={call.args[0]:call.args[1] for call in reader.option.call_args_list}
        self.assertEqual(calls['ignoreExtension'],'true')
        self.assertEqual(calls['sampleFilePath'],sample.as_uri())
        self.assertEqual(calls['sampleIdColumn'],'ID_2')
