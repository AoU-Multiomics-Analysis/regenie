"""Exercise the real WDL command boundary with cloud File values."""
import os
import shlex
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
import WDL
from test_wdl import LocalStdLib

ROOT=Path(__file__).resolve().parents[1]

class GlowWdlTests(unittest.TestCase):
    def test_workflow_has_typed_files_and_no_workflow_writes(self):
        doc=WDL.load(str(ROOT/"workflows/glow_regenie.wdl"))
        from importlib.util import spec_from_file_location,module_from_spec
        spec=spec_from_file_location("check",ROOT/"tools/check_terra_wdl.py")
        checker=module_from_spec(spec); spec.loader.exec_module(checker)
        self.assertEqual(checker.workflow_write_errors(doc),[])
        self.assertEqual(checker.manifest_output_errors(doc),[])
        self.assertEqual(doc.workflow.name,"GlowRegenie")

    def test_command_time_localization_and_shell_quoting(self):
        doc=WDL.load(str(ROOT/"workflows/glow_regenie.wdl"))
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            localized=root/"localized ' $(touch INJECTED)"; localized.mkdir()
            for name in ("data.pgen","data.pvar","data.psam","phenotypes.tsv","sample_map.tsv","covariates.tsv",
                         "batches.tsv","model.bgen","model.sample","variants.tsv"):
                (localized/name).write_text("data")
            for task in doc.tasks:
                with self.subTest(task=task.name):
                    supplied={}
                    for declaration in task.inputs:
                        name=declaration.name; kind=declaration.type
                        if isinstance(kind,WDL.Type.File):
                            supplied[name]="gs://bucket/"+("data."+name if name in ("pgen","pvar","psam") else "phenotypes.tsv")
                        elif isinstance(kind,WDL.Type.Array) and isinstance(kind.item_type,WDL.Type.File):
                            supplied[name]=["gs://bucket/model.bgen","gs://bucket/variants.tsv"]
                        elif declaration.expr is None and not kind.optional:
                            if isinstance(kind,WDL.Type.String): supplied[name]="value ' $(touch INJECTED)"
                            elif isinstance(kind,WDL.Type.Int): supplied[name]=1
                    env=WDL.values_from_json(supplied,task.available_inputs)
                    mapped=WDL.Env.Bindings()
                    for binding in env:
                        value=WDL.Value.rewrite_paths(binding.value,lambda f:str(localized/f.value.rsplit("/",1)[-1]))
                        mapped=mapped.bind(binding.name,value)
                    env=mapped; std=LocalStdLib("1.0",localized)
                    for declaration in [*task.inputs,*task.postinputs]:
                        if not env.has_binding(declaration.name):
                            value=declaration.expr.eval(env,std) if declaration.expr else WDL.Value.Null()
                            env=env.bind(declaration.name,value)
                    command=task.command.eval(env,std).value
                    if task.name=="RunGlowAnalysis":
                        self.assertNotIn("--project",command)
                        self.assertNotIn("--service-account",command)
                    probe=root/"probe.py"
                    probe.write_text("""import pathlib,sys
args=sys.argv[1:]
for index,flag in enumerate(args[:-1]):
    if flag in ('--pgen','--pvar','--psam','--sample-map','--phenotypes','--covariates','--phenotype-batches','--model-bgen','--model-sample','--model-variants'):
        assert pathlib.Path(args[index+1]).is_file(), (flag,args[index+1])
    if flag.endswith('-list'):
        rows=pathlib.Path(args[index+1]).read_text().splitlines()
        assert rows and all(pathlib.Path(row).is_file() for row in rows)
print('LOCALIZATION_OK')
""")
                    for module in ("prepare","convert","controller"):
                        command=command.replace("python3 -m scripts.glow."+module,shlex.quote(sys.executable)+" "+shlex.quote(str(probe)))
                    work=root/task.name; work.mkdir()
                    result=subprocess.run(["bash","-c",command],cwd=work,text=True,capture_output=True)
                    self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                    self.assertIn("LOCALIZATION_OK",result.stdout)
                    self.assertFalse((work/"INJECTED").exists())
