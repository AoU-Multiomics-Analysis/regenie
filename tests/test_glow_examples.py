import json
import unittest
from pathlib import Path
import WDL
import yaml
ROOT=Path(__file__).resolve().parents[1]
class GlowExampleTests(unittest.TestCase):
    def test_both_examples_preserve_file_types_and_cloud_controls(self):
        doc=WDL.load(str(ROOT/"workflows/glow_regenie.wdl"))
        for kind in ("quantitative","binary"):
            values=json.loads((ROOT/f"examples/glow.{kind}.inputs.json").read_text())
            bound=WDL.values_from_json(values,doc.workflow.available_inputs,doc.workflow.required_inputs,namespace="GlowRegenie")
            self.assertEqual(bound.resolve("trait_type").value,kind)
            self.assertTrue(all(isinstance(v,WDL.Value.File) for v in bound.resolve("step1_genotypes").value.values()))
            for name in ("project","region","subnetwork_uri","service_account","staging_bucket","staging_prefix","run_id","docker_image","spark_image"):
                self.assertTrue(bound.resolve(name).value)
            self.assertEqual(values["GlowRegenie.docker_image"].split("@")[1],values["GlowRegenie.spark_image"].split("@")[1])
    def test_dockstore_has_both_workflows(self):
        descriptor=yaml.safe_load((ROOT/".dockstore.yml").read_text())
        self.assertEqual({item["name"] for item in descriptor["workflows"]},{"regenie","glow-regenie"})
        for item in descriptor["workflows"]:
            self.assertTrue((ROOT/item["primaryDescriptorPath"].lstrip("/")).is_file())
            for path in item["testParameterFiles"]:
                self.assertTrue((ROOT/path.lstrip("/")).is_file())
