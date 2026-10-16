"""Check local Spark resources, File handling, and verified output collection."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from scripts.glow.contracts import ScientificOptions, InputArtifact, OutputArtifact
from scripts.glow.storage import checksum, publish_outputs
from scripts.glow import local_runner

class LocalRunnerTests(unittest.TestCase):
    def test_resource_limits_and_single_vm_command(self):
        for values in ({'cpu':0}, {'memory_gb':4,'driver_memory_gb':4}, {'parallelism':0}, {'shuffle_partitions':0}):
            with self.subTest(values=values), self.assertRaises(ValueError):
                local_runner.LocalOptions(**values)
        options=local_runner.LocalOptions(cpu=8,memory_gb=32,driver_memory_gb=16,parallelism=20,shuffle_partitions=30)
        command=local_runner.spark_command(options,Path('/work/scratch'))
        self.assertEqual(command[command.index('--master')+1],'local[8]')
        self.assertEqual(command[command.index('--driver-memory')+1],'16g')
        self.assertIn('spark.default.parallelism=20',command)
        self.assertIn('spark.sql.shuffle.partitions=30',command)
        self.assertIn('spark.sql.execution.arrow.pyspark.enabled=true',command)
        self.assertNotIn('--executor-memory',command)

    def test_cloud_paths_and_misaligned_arrays_fail_before_spark(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); good=root/'good'; good.write_text('data')
            for name,content in [('bgen','gs://bucket/input.bgen\n'),('sample',str(good)+'\n'),('variants',str(good)+'\n')]:
                (root/name).write_text(content)
            argv=[]
            for name in ('model-bgen','model-sample','model-variants','sample-map','phenotypes','covariates','phenotype-batches'):
                argv+=['--'+name,str(good)]
            for name in ('bgen','sample','variants'): argv+=['--association-'+name+'-list',str(root/name)]
            argv+=['--trait-type','quantitative','--out-dir',str(root/'out')]
            with patch.object(local_runner,'run_local') as run:
                with self.assertRaisesRegex(ValueError,'localization'): local_runner.main(argv)
                run.assert_not_called()
            (root/'bgen').write_text(str(good)+'\n'+str(good)+'\n')
            with self.assertRaisesRegex(ValueError,'differ in length'): local_runner.main(argv)

    def test_subprocess_receives_file_uris_and_outputs_are_collected(self):
        with tempfile.TemporaryDirectory(prefix="glow ' ") as directory:
            root=Path(directory); source=root/'source'; source.write_text('data')
            inputs=[InputArtifact(name,source) for name in ('model/bgen','model/sample','model/variants','sample_map','phenotypes','covariates','phenotype_batches','association/0/bgen','association/0/sample','association/0/variants')]
            def fake_run(command,**kwargs):
                def value(flag): return command[command.index(flag)+1]
                self.assertTrue(value('--model-bgen').startswith('file:///'))
                self.assertIn('%27',value('--model-bgen'))
                prefix=value('--output-prefix'); outputs=[]
                for kind,sub in [('association','associations'),('prediction','predictions'),('metrics','metrics')]:
                    path=root/'artifact'; path.write_text(kind)
                    from scripts.glow.storage import upload_file
                    uri=prefix+'/results/'+sub+'/result.txt'; upload_file(None,uri,path)
                    outputs.append(OutputArtifact(kind,uri,sub+'/result.txt',checksum(path),path.stat().st_size))
                publish_outputs(None,prefix,value('--identity-sha256'),outputs,{'test':True})
                return type('Result',(),{'returncode':0})()
            with patch.object(local_runner.subprocess,'run',side_effect=fake_run):
                local_runner.run_local(local_runner.LocalOptions(cpu=2,memory_gb=4,driver_memory_gb=1),ScientificOptions('quantitative'),inputs,root/'out')
            report=json.loads((root/'out/execution.json').read_text())
            self.assertEqual(report['state'],'SUCCEEDED')
            self.assertEqual(report['spark_master'],'local[2]')
            self.assertEqual((root/'out/associations/result.txt').read_text(),'association')
            self.assertTrue((root/'out/output_manifest.json').is_file())

    def test_failed_spark_process_fails_task_and_writes_report(self):
        import subprocess
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); source=root/'input'; source.write_text('data')
            names=('model/bgen','model/sample','model/variants','sample_map','phenotypes','covariates','phenotype_batches','association/0/bgen','association/0/sample','association/0/variants')
            inputs=[InputArtifact(name,source) for name in names]
            with patch.object(local_runner.subprocess,'run',side_effect=subprocess.CalledProcessError(2,['spark-submit'])):
                with self.assertRaises(subprocess.CalledProcessError):
                    local_runner.run_local(local_runner.LocalOptions(),ScientificOptions('quantitative'),inputs,root/'out')
            self.assertEqual(json.loads((root/'out/execution.json').read_text())['state'],'FAILED')
            self.assertFalse((root/'out/output_manifest.json').exists())
