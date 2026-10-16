import tempfile
import unittest
from pathlib import Path
from dataclasses import replace
from google.cloud import dataproc_v1 as d
from google.api_core.exceptions import NotFound, DeadlineExceeded
from scripts.glow.contracts import ScientificOptions, InputArtifact, StagedArtifact
from scripts.glow.controller import CloudOptions, make_batch, run_controller
from test_glow_storage import Client, StorageTests

DIGEST="a"*64
def cloud():
    return CloudOptions(project="project-one",region="us-central1",subnetwork_uri="projects/project-one/regions/us-central1/subnetworks/default",
                        service_account="worker@project-one.iam.gserviceaccount.com",staging_bucket="test-bucket",
                        staging_prefix="gs://bucket/base",run_id="test-run",terra_image="ghcr.io/org/glow@sha256:"+DIGEST,
                        spark_image="us-central1-docker.pkg.dev/project-one/repo/glow@sha256:"+DIGEST)

class ControllerTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name); self.path=self.root/"input"; self.path.write_bytes(b"test")
        self.names=["model/bgen","model/sample","model/variants","sample_map","phenotypes","covariates","phenotype_batches","association/0/bgen","association/0/sample","association/0/variants","application/spark_main.py"]
    def test_request_uses_typed_runtime_and_named_arguments(self):
        staged=[StagedArtifact(n,"gs://bucket/base/run/inputs/"+n,"b"*64,4,1) for n in self.names]
        batch=make_batch(cloud(),staged,ScientificOptions("binary"),"c"*64)
        self.assertEqual(batch.runtime_config.version,"1.2")
        self.assertEqual(batch.environment_config.execution_config.ttl.total_seconds(),21600)
        self.assertIn("--association-bgen",list(batch.pyspark_batch.args))
        self.assertIn("--trait-type",list(batch.pyspark_batch.args))
        self.assertNotIn("--config",list(batch.pyspark_batch.args))
        self.assertEqual(batch.runtime_config.properties["spark.dynamicAllocation.maxExecutors"],"4")
    def test_invalid_cloud_settings_fail_before_clients(self):
        for changes in ({"run_id":"../bad"},{"spark_image":"ghcr.io/org/glow@sha256:"+DIGEST},{"poll_seconds":1},
                        {"driver_cores":1},{"driver_memory_mb":50000},{"staging_prefix":"gs://bucket"},
                        {"spark_image":"us-central1-docker.pkg.dev/project-one/repo/glow@sha256:"+"b"*64}):
            with self.subTest(changes=changes):
                with self.assertRaises(ValueError): replace(cloud(),**changes)
    def test_uncertain_create_is_not_submitted_twice(self):
        client=Client()
        from scripts.glow.storage import publish_outputs
        from scripts.glow.contracts import OutputArtifact
        import hashlib
        class BatchClient:
            def __init__(self): self.batch=None; self.creates=0
            def get_batch(self,**kw):
                if self.batch is None: raise NotFound("not found")
                return self.batch
            def create_batch(self,request,**kw):
                self.creates+=1
                self.batch=d.Batch(request["batch"])
                self.batch.state=d.Batch.State.SUCCEEDED
                self.batch.operation="projects/project-one/regions/us-central1/operations/one"
                args=list(self.batch.pyspark_batch.args)
                identity=args[args.index("--identity-sha256")+1]
                uri="gs://bucket/base/test-run"
                outputs=[]
                for kind,directory in (("association","associations"),("prediction","predictions"),("metrics","metrics")):
                    remote=uri+"/results/"+kind
                    client.bucket("bucket").blob("base/test-run/results/"+kind).upload_from_string(b"result")
                    outputs.append(OutputArtifact(kind,remote,directory+"/"+kind,hashlib.sha256(b"result").hexdigest(),6))
                publish_outputs(client,uri,identity,outputs,{})
                raise DeadlineExceeded("lost response")
        remote=BatchClient()
        result=run_controller(cloud(),ScientificOptions("quantitative"),[InputArtifact(n,self.path) for n in self.names],
                              self.root/"out",remote,client,sleeper=lambda _:None)
        self.assertEqual(len(result),3)
        self.assertEqual(remote.creates,1)
        again=run_controller(cloud(),ScientificOptions("quantitative"),[InputArtifact(n,self.path) for n in self.names],
                              self.root/"out2",remote,client,sleeper=lambda _:None)
        self.assertEqual(len(again),3)
        self.assertEqual(remote.creates,1)
