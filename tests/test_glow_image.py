"""Check build configuration invariants that affect publication and runtime."""
import unittest
from pathlib import Path
import yaml
ROOT=Path(__file__).resolve().parents[1]

class ImageTests(unittest.TestCase):
    def test_environment_uses_exact_direct_dependency_pins(self):
        for name in ("environment.yml","build-environment.yml"):
            values=yaml.safe_load((ROOT/"envs/glow"/name).read_text())
            self.assertEqual(values["channels"][:2],["conda-forge","bioconda"])
            for package in values["dependencies"]:
                self.assertIsInstance(package,str)
                self.assertEqual(package.count("="),1,package)
                self.assertNotIn("*",package)
    def test_automatic_builds_are_limited_to_included_python_scripts(self):
        workflow=yaml.safe_load((ROOT/".github/workflows/glow-image.yml").read_text())
        triggers=workflow.get("on",workflow.get(True))
        self.assertIn("workflow_dispatch",triggers)
        for event in ("push","pull_request"):
            paths=triggers[event]["paths"]
            self.assertTrue(paths)
            self.assertTrue(all(p.startswith("scripts/glow/") and p.endswith(".py") for p in paths))
        steps=workflow["jobs"]["build-test-publish"]["steps"]
        publish=next(i for i,s in enumerate(steps) if s.get("name")=="Publish tested production image")
        test=next(i for i,s in enumerate(steps) if s.get("name")=="Run synthetic Spark tests")
        self.assertLess(test,publish)
        self.assertNotIn("always()",steps[publish].get("if",""))

    def test_ci_executors_have_at_least_one_gibibyte_heap(self):
        from scripts.glow.smoke_test import spark_command
        for master in ("local-cluster[1,1,1536]", "local-cluster[2,1,1536]"):
            command = spark_command(master)
            self.assertEqual(command[command.index("--executor-memory") + 1], "1g")
