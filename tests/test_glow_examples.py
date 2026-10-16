import json
import unittest
from pathlib import Path
import WDL
import yaml
ROOT=Path(__file__).resolve().parents[1]
class GlowExampleTests(unittest.TestCase):
    def test_examples_preserve_file_types_and_execution_mode(self):
        for cloud in (False,True):
            doc=WDL.load(str(ROOT/"workflows"/("glow_regenie_dataproc.wdl" if cloud else "glow_regenie.wdl")))
            namespace=doc.workflow.name
            for kind in ("quantitative","binary"):
                suffix="dataproc." if cloud else ""
                values=json.loads((ROOT/f"examples/glow.{suffix}{kind}.inputs.json").read_text())
                bound=WDL.values_from_json(values,doc.workflow.available_inputs,doc.workflow.required_inputs,namespace=namespace)
                self.assertEqual(bound.resolve("trait_type").value,kind)
                self.assertTrue(all(isinstance(v,WDL.Value.File) for v in bound.resolve("step1_genotypes").value.values()))
                self.assertNotIn(namespace+".project",values)
                self.assertNotIn(namespace+".service_account",values)
                if cloud:
                    for name in ("region","subnetwork_uri","staging_bucket","staging_prefix","run_id","spark_image"):
                        self.assertTrue(bound.resolve(name).value)
                    self.assertEqual(bound.resolve("docker_image").value.split("@")[1],bound.resolve("spark_image").value.split("@")[1])
                else:
                    self.assertFalse(any(key.split(".")[-1] in ("region","spark_image","run_id") for key in values))
    def test_dockstore_has_all_workflows(self):
        descriptor=yaml.safe_load((ROOT/".dockstore.yml").read_text())
        self.assertEqual({item["name"] for item in descriptor["workflows"]},{"regenie","glow-regenie","glow-regenie-dataproc"})
        for item in descriptor["workflows"]:
            self.assertTrue((ROOT/item["primaryDescriptorPath"].lstrip("/")).is_file())
            for path in item["testParameterFiles"]:
                self.assertTrue((ROOT/path.lstrip("/")).is_file())
    def test_dataproc_identity_remains_an_optional_override(self):
        doc=WDL.load(str(ROOT/"workflows/glow_regenie_dataproc.wdl"))
        task=next(task for task in doc.tasks if task.name=="RunGlowAnalysis")
        for name in ("project","service_account"):
            self.assertTrue(next(item for item in task.inputs if item.name==name).type.optional)
