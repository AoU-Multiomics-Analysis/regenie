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
            for name in ("region","subnetwork_uri","staging_bucket","staging_prefix","run_id","docker_image","spark_image"):
                self.assertTrue(bound.resolve(name).value)
            self.assertNotIn("GlowRegenie.project",values)
            self.assertNotIn("GlowRegenie.service_account",values)
            self.assertEqual(values["GlowRegenie.docker_image"].split("@")[1],values["GlowRegenie.spark_image"].split("@")[1])
    def test_dockstore_has_both_workflows(self):
        descriptor=yaml.safe_load((ROOT/".dockstore.yml").read_text())
        self.assertEqual({item["name"] for item in descriptor["workflows"]},{"regenie","glow-regenie"})
        for item in descriptor["workflows"]:
            self.assertTrue((ROOT/item["primaryDescriptorPath"].lstrip("/")).is_file())
            for path in item["testParameterFiles"]:
                self.assertTrue((ROOT/path.lstrip("/")).is_file())

    def test_project_and_service_account_are_optional_overrides(self):
        doc=WDL.load(str(ROOT/"workflows/glow_regenie.wdl"))
        values=json.loads((ROOT/"examples/glow.quantitative.inputs.json").read_text())
        for name in ("project","service_account"):
            values.pop("GlowRegenie."+name,None)
        WDL.values_from_json(values,doc.workflow.available_inputs,doc.workflow.required_inputs,namespace="GlowRegenie")
        task=next(task for task in doc.tasks if task.name=="RunGlowAnalysis")
        for name in ("project","service_account"):
            self.assertTrue(next(item for item in task.inputs if item.name==name).type.optional)
