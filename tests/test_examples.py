import json
import unittest
from pathlib import Path

import WDL
import yaml

ROOT = Path(__file__).resolve().parents[1]


class ExampleTests(unittest.TestCase):
    def test_image_can_use_workflow_default_or_explicit_override(self):
        document = WDL.load(str(ROOT / "workflows/regenie.wdl"))
        inputs = json.loads((ROOT / "examples/quantitative.inputs.json").read_text())
        inputs.pop("Regenie.docker_image", None)
        WDL.values_from_json(inputs, document.workflow.available_inputs,
                             document.workflow.required_inputs, namespace="Regenie")
        inputs["Regenie.docker_image"] = "registry.example.org/regenie@sha256:" + "a" * 64
        values = WDL.values_from_json(inputs, document.workflow.available_inputs,
                                      document.workflow.required_inputs, namespace="Regenie")
        self.assertEqual(values.resolve("docker_image").value, inputs["Regenie.docker_image"])

    def test_workflow_inputs_accept_missing_null_or_blank_phenotype_columns(self):
        document = WDL.load(str(ROOT / "workflows/regenie.wdl"))
        original = json.loads((ROOT / "examples/quantitative.inputs.json").read_text())
        for selection in (None, "", "   "):
            with self.subTest(selection=selection):
                inputs = dict(original, **{"Regenie.phenotype_columns": selection})
                WDL.values_from_json(inputs, document.workflow.available_inputs,
                                     document.workflow.required_inputs, namespace="Regenie")
        inputs = {key: value for key, value in original.items() if key != "Regenie.phenotype_columns"}
        WDL.values_from_json(inputs, document.workflow.available_inputs,
                             document.workflow.required_inputs, namespace="Regenie")

    def test_examples_supply_required_inputs_and_preserve_file_types(self):
        document = WDL.load(str(ROOT / "workflows/regenie.wdl"))
        for path in (ROOT / "examples").glob("*.json"):
            with self.subTest(example=path.name):
                values = WDL.values_from_json(json.loads(path.read_text()), document.workflow.available_inputs,
                                              document.workflow.required_inputs, namespace="Regenie")
                genotypes = values.resolve("step1_genotypes")
                self.assertTrue(all(isinstance(value, WDL.Value.File) for value in genotypes.value.values()))
                associations = values.resolve("association_genotypes")
                self.assertGreater(len(associations.value), 0)
                self.assertTrue(all(isinstance(file, WDL.Value.File) for item in associations.value for file in item.value.values()))

    def test_dockstore_descriptor_points_to_existing_files(self):
        descriptor = yaml.safe_load((ROOT / ".dockstore.yml").read_text())
        workflow = descriptor["workflows"][0]
        self.assertTrue((ROOT / workflow["primaryDescriptorPath"].lstrip("/")).is_file())
        for path in workflow["testParameterFiles"]:
            self.assertTrue((ROOT / path.lstrip("/")).is_file())


if __name__ == "__main__":
    unittest.main()
