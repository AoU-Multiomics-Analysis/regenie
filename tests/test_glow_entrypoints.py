"""Spark executes these files directly, outside Python module mode."""
import os
import subprocess
import sys
import unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

class EntryPointTests(unittest.TestCase):
    def test_spark_entry_files_accept_direct_execution(self):
        for name in ("spark_main.py", "reference_test.py"):
            with self.subTest(name=name):
                result = subprocess.run([sys.executable, str(ROOT / "scripts/glow" / name), "--help"],
                                        env=dict(os.environ, PYTHONPATH=str(ROOT)),
                                        capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("usage:", result.stdout)
