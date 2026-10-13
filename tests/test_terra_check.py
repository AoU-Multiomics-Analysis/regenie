import importlib.util
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


if __name__ == "__main__":
    unittest.main()
