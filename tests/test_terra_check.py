import importlib.util
import tempfile
import unittest
from pathlib import Path

import WDL

ROOT = Path(__file__).resolve().parents[1]


class TerraCheckTests(unittest.TestCase):
    def setUp(self):
        script = ROOT / "tools/check_terra_wdl.py"
        self.assertTrue(script.is_file(), "Terra scope checker is missing")
        spec = importlib.util.spec_from_file_location("check_terra_wdl", script)
        self.checker = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.checker)

    def check(self, source):
        return self.checker.workflow_write_errors(WDL.parse_document("version 1.0\n" + source))

    def test_workflow_declarations_cannot_create_files(self):
        self.assertEqual(len(self.check('workflow w { File bad = write_lines(["a"]) }')), 1)

    def test_scatter_conditional_call_and_output_expressions_are_checked(self):
        source = '''
        task t { input { File f } command { cat ~{f} } }
        workflow w {
          scatter (i in [1, 2]) {
            if (i == 1) { call t { input: f = write_lines(["a"]) } }
          }
          output { File bad = write_json({"x": 1}) }
        }
        '''
        self.assertEqual(len(self.check(source)), 2)

    def test_task_scope_and_task_commands_can_create_files(self):
        source = '''
        task t {
          input { Array[String] items }
          File list = write_lines(items)
          command { cat ~{list} }
          output { File f = stdout() }
        }
        workflow w { call t { input: items = ["a"] } }
        '''
        self.assertEqual(self.check(source), [])

    def test_comments_and_string_literals_are_not_false_positives(self):
        self.assertEqual(self.check('workflow w { String note = "write_lines(x)" # write_json(x)\n }'), [])

    def output_check(self, source):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.wdl"
            path.write_text("version 1.0\n" + source)
            return self.checker.manifest_output_errors(WDL.load(str(path)))

    def test_file_arrays_and_nested_file_structs_cannot_come_from_manifests(self):
        for output in ('Array[File] files = read_lines("paths.txt")',
                       'Array[Group] files = read_json("paths.json")',
                       'Map[String, Array[Group]] metadata = read_json("paths.json")\n Array[Group] files = metadata["groups"]'):
            with self.subTest(output=output):
                errors = self.output_check('''struct Group { Array[File] files }
                    task t { command { echo test } output { ''' + output + ''' } }
                    workflow w { call t }''')
                self.assertTrue(errors)
                self.assertIn("File outputs", errors[-1])

    def test_manifest_metadata_and_explicit_file_outputs_are_allowed(self):
        self.assertEqual(self.output_check('''struct Batch { Int index String names }
            task t { command { echo test } output {
                Array[Batch] batches = read_json("batches.json")
                Array[File] files = glob("results/*.loco")
                File manifest = "batches.json"
            } }
            workflow w { call t }'''), [])

    def test_production_task_outputs_do_not_depend_on_paths_inside_manifests(self):
        self.assertEqual(self.checker.manifest_output_errors(WDL.load(str(ROOT / "workflows/regenie.wdl"))), [])


if __name__ == "__main__":
    unittest.main()
