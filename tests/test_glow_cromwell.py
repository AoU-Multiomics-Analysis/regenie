"""Check the production graph with Cromwell's File mapper and no cloud calls."""
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

@unittest.skipUnless(os.environ.get("CROMWELL_JAR"),"Set CROMWELL_JAR to run the engine localization test")
class GlowCromwellTests(unittest.TestCase):
    def test_local_source_arrays_and_generated_outputs_are_localized(self):
        self.run_probe(local=True)

    def test_dataproc_source_arrays_and_generated_outputs_are_localized(self):
        self.run_probe(local=False)

    def run_probe(self,local):
        with tempfile.TemporaryDirectory(prefix="glow-cromwell-") as directory:
            root=Path(directory); source=root/"source"; source.mkdir()
            for index in range(2):
                for extension in ("pgen","pvar","psam"):
                    (source/f"source{index}.{extension}").write_text(f"source {index}")
            for name in ("variants","phenotypes","covariates","keep","remove"):
                (source/name).write_text(name)
            probe=root/"probe.py"
            probe.write_text("""import json,pathlib,sys
args=sys.argv[1:]
def value(flag): return args[args.index(flag)+1]
mode=value('--probe-task')
for index,flag in enumerate(args[:-1]):
    if flag in ('--pgen','--pvar','--psam','--sample-map','--phenotypes','--covariates','--keep','--remove','--variant-ids','--model-bgen','--model-sample','--model-variants','--phenotype-batches'):
        path=pathlib.Path(args[index+1])
        assert path.is_file() and '/inputs/' in str(path), (flag,str(path))
out=pathlib.Path(value('--out-dir')); out.mkdir()
if mode=='prepare':
    for name in ('sample_map.tsv','phenotypes.tsv','covariates.tsv','phenotype_batches.tsv','preparation.json'):
        (out/name).write_text(name)
elif mode=='convert':
    index=value('--source-index'); role=value('--role')
    for name in ('genotypes.bgen','genotypes.sample','variants.tsv','conversion.json','conversion.log'):
        (out/name).write_text(role+index)
else:
    for kind in ('bgen','sample','variants'):
        rows=pathlib.Path(value('--association-'+kind+'-list')).read_text().splitlines()
        assert len(rows)==2
        assert [pathlib.Path(p).read_text() for p in rows]==['association0','association1']
        assert all('/inputs/' in p for p in rows)
    for sub in ('associations','predictions','metrics','logs'): (out/sub).mkdir()
    (out/'associations'/'source0.tsv.gz').write_text('source0')
    (out/'associations'/'source1.tsv.gz').write_text('source1')
    (out/'predictions'/'loco.tsv.gz').write_text('prediction')
    (out/'metrics'/'counts.json').write_text('{}')
    (out/'logs'/'application.log').write_text('log')
    for name in ('execution.json','output_manifest.json'): (out/name).write_text('{}')
print('COMMAND_FILE_LOCALIZATION_OK',mode)
""")
            for name in ("glow_regenie.wdl","glow_regenie_dataproc.wdl"):
                wdl=(ROOT/"workflows"/name).read_text()
                wdl=re.sub(r"(?m)^\s+docker: docker_image\n","\n",wdl)
                for mode in ("prepare","convert","controller","local_runner"):
                    wdl=wdl.replace("python3 -m scripts.glow."+mode,shlex.quote(sys.executable)+" "+shlex.quote(str(probe))+" --probe-task "+mode)
                (root/name).write_text(wdl)
            workflow=root/("glow_regenie.wdl" if local else "glow_regenie_dataproc.wdl")
            namespace="GlowRegenie" if local else "GlowRegenieDataproc"
            genotypes=[{ext:str(source/f"source{i}.{ext}") for ext in ("pgen","pvar","psam")} for i in range(2)]
            values={"step1_genotypes":genotypes[0],"step1_variant_ids":str(source/"variants"),
                    "association_genotypes":genotypes,"phenotypes":str(source/"phenotypes"),
                    "covariates":str(source/"covariates"),"keep_samples":str(source/"keep"),"remove_samples":str(source/"remove"),
                    "trait_type":"quantitative","docker_image":"test","spark_image":"test","project":"project-one",
                    "region":"us-central1","subnetwork_uri":"test","service_account":"test","staging_bucket":"test",
                    "staging_prefix":"gs://test/run","run_id":"test-run"}
            if local:
                for name in ("spark_image","project","region","subnetwork_uri","service_account","staging_bucket","staging_prefix","run_id"):
                    values.pop(name)
            inputs=root/"inputs.json"; inputs.write_text(json.dumps({namespace+"."+k:v for k,v in values.items()}))
            config=root/"cromwell.conf"; config.write_text('include required(classpath("application"))\nbackend.default="Local"\ncall-caching.enabled=false\n')
            metadata=root/"metadata.json"
            result=subprocess.run(["java","-Dconfig.file="+str(config),"-jar",os.environ["CROMWELL_JAR"],"run",str(workflow),"-i",str(inputs),"-m",str(metadata)],
                                  cwd=root,text=True,capture_output=True,timeout=300)
            self.assertEqual(result.returncode,0,(result.stdout+result.stderr)[-12000:])
            output=json.loads(metadata.read_text())
            self.assertEqual(output["status"],"Succeeded")
            self.assertEqual(len(output["outputs"][namespace+".associations"]),2)
            self.assertEqual(len(output["outputs"][namespace+".association_sets"]),2)
            logs=output["outputs"][namespace+".conversion_logs"]
            self.assertEqual([Path(p).read_text() for p in logs],["model0","association0","association1"])
            self.assertTrue(all(Path(p).is_file() for p in output["outputs"][namespace+".associations"]))
