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
