import unittest
from types import SimpleNamespace
import pandas as pd
from scripts.glow.analysis import FittedBatch, predict_chromosome, validate_tables, validate_binary_blocks
from scripts.glow.contracts import ScientificOptions

class AnalysisTests(unittest.TestCase):
    def test_sample_order_and_trait_order_must_match(self):
        pheno=pd.DataFrame({"a":[1.,2.],"b":[2.,3.]},index=["s1","s2"])
        cov=pd.DataFrame(index=["s2","s1"])
        with self.assertRaisesRegex(ValueError,"sample"):
            validate_tables(pheno,cov,["s1","s2"])
        with self.assertRaisesRegex(ValueError,"sample"):
            validate_tables(pheno,pd.DataFrame(index=pheno.index),["s2","s1"])

    def test_wrong_chromosome_fails_before_prediction(self):
        fit=FittedBatch(None,None,["s1"],["a"],"quantitative",("1","2"))
        with self.assertRaisesRegex(ValueError,"chromosome"):
            predict_chromosome(fit,"3")

    def test_prediction_reindexes_samples_and_uses_binary_linear_offset(self):
        calls=[]
        def predict(**kw):
            calls.append(kw)
            return pd.DataFrame({"a":[2.,1.]},index=pd.MultiIndex.from_tuples([("s2","1"),("s1","1")]))
        fit=FittedBatch(SimpleNamespace(transform_loco=predict),None,["s1","s2"],["a"],"binary",("1","2"))
        result=predict_chromosome(fit,"1")
        self.assertEqual(result["a"].tolist(),[1.,2.])
        self.assertEqual(calls,[{"chromosomes":["1"],"response":"linear"}])

    def test_binary_training_fold_cannot_lose_a_class(self):
        pheno=pd.DataFrame({"a":[0,0,1,1]},index=["s1","s2","s3","s4"])
        with self.assertRaisesRegex(ValueError,"class"):
            validate_binary_blocks(pheno,{"0":["s1","s2"],"1":["s3","s4"]})
        validate_binary_blocks(pheno,{"0":["s1","s3"],"1":["s2","s4"]})
