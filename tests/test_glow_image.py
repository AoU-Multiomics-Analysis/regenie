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

    def test_local_image_contains_spark_and_is_tested_before_publication(self):
        dockerfile=(ROOT/'envs/glow/Dockerfile').read_text()
        self.assertIn('FROM production AS terra-local',dockerfile)
        self.assertIn('pyspark=3.5.1 openjdk=17.0.15',dockerfile)
        self.assertIn('FROM terra-local AS ci-test',dockerfile)
        workflow=yaml.safe_load((ROOT/'.github/workflows/glow-image.yml').read_text())
        steps=workflow['jobs']['build-test-publish']['steps']
        publish=next(s for s in steps if s.get('name')=='Publish tested production image')['run']
        self.assertIn('docker push "$LOCAL_IMAGE_TAG"',publish)
        self.assertIn('docker tag glow-local "$LOCAL_IMAGE_TAG"',publish)

    def test_image_change_gate_skips_wdl_edits_and_detects_included_scripts(self):
        import os
        import subprocess
        import tempfile
        workflow=yaml.safe_load((ROOT/'.github/workflows/glow-image.yml').read_text())
        steps=workflow['jobs']['build-test-publish']['steps']
        gate=next(step for step in steps if step.get('id')=='changes')['run']
        with tempfile.TemporaryDirectory() as directory:
            repo=Path(directory); (repo/'scripts/glow').mkdir(parents=True)
            script=repo/'scripts/glow/app.py'; script.write_text('one')
            wdl=repo/'workflow.wdl'; wdl.write_text('one')
            def git(*args):
                return subprocess.check_output(['git','-c','user.name=Test','-c','user.email=test@example.invalid',*args],cwd=repo,text=True).strip()
            git('init','-q'); git('add','.'); git('commit','-qm','initial'); before=git('rev-parse','HEAD')
            for changed,expected in ((wdl,'false'),(script,'true')):
                changed.write_text('two'); git('add','.'); git('commit','-qm','change'); after=git('rev-parse','HEAD')
                output=repo/'output'; output.unlink(missing_ok=True)
                result=subprocess.run(['bash','-e','-o','pipefail','-c',gate],cwd=repo,text=True,capture_output=True,
                                      env=dict(os.environ,EVENT_NAME='pull_request',BEFORE_SHA=before,AFTER_SHA=after,GITHUB_OUTPUT=str(output)))
                self.assertEqual(result.returncode,0,result.stderr)
                self.assertIn('build='+expected,output.read_text()); before=after
