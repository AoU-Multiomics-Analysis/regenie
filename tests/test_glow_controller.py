import tempfile
import unittest
from pathlib import Path
from dataclasses import replace
from google.cloud import dataproc_v1 as d
from google.api_core.exceptions import NotFound, DeadlineExceeded
from scripts.glow.contracts import ScientificOptions, InputArtifact, StagedArtifact
from scripts.glow.controller import CloudOptions, make_batch, run_controller
from test_glow_storage import Client

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
        self.assertEqual(batch.runtime_config.properties["spark.sql.execution.arrow.pyspark.enabled"],"true")
        self.assertEqual(batch.runtime_config.properties["spark.sql.execution.arrow.pyspark.fallback.enabled"],"false")
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

    def active_remote(self, state=d.Batch.State.RUNNING):
        class Remote:
            def __init__(self):
                self.batch = None
                self.creates = 0
                self.cancellations = []
                self.fail_cancel = False
                self.outage = False
            def get_batch(self, **kwargs):
                if self.outage:
                    from google.api_core.exceptions import ServiceUnavailable
                    raise ServiceUnavailable("status unavailable")
                if self.batch is None: raise NotFound("missing")
                return self.batch
            def create_batch(self, request, **kwargs):
                self.creates += 1
                self.batch = d.Batch(request["batch"])
                self.batch.state = state
                self.batch.operation = "projects/project-one/regions/us-central1/operations/one"
            def cancel_operation(self, request, **kwargs):
                self.cancellations.append(request["name"])
                if self.fail_cancel: raise RuntimeError("cancel unavailable")
        return Remote()

    def run_remote(self, remote, **kwargs):
        return run_controller(cloud(), ScientificOptions("quantitative"),
                              [InputArtifact(n, self.path) for n in self.names],
                              self.root / "out", remote, Client(), **kwargs)

    def test_stop_signal_with_status_outage_requests_cancellation(self):
        import signal
        remote = self.active_remote()
        def sleep(seconds):
            if seconds == 30:
                remote.outage = True
                signal.getsignal(signal.SIGTERM)(signal.SIGTERM, None)
        with self.assertRaisesRegex(RuntimeError, "stopped"):
            self.run_remote(remote, sleeper=sleep)
        self.assertEqual(remote.creates, 1)
        self.assertEqual(remote.cancellations, [remote.batch.operation])
        report = __import__("json").loads((self.root / "out/execution.json").read_text())
        self.assertEqual(report["cancellation"], "requested")

    def test_deadline_records_cancellation_failure_and_normal_polling(self):
        remote = self.active_remote()
        remote.fail_cancel = True
        now = [0]
        sleeps = []
        def sleep(seconds):
            sleeps.append(seconds)
            if seconds == 30: now[0] = 24000
        with self.assertRaisesRegex(RuntimeError, "timed out"):
            self.run_remote(remote, sleeper=sleep, clock=lambda: now[0])
        self.assertEqual(sleeps, [1, 30])
        self.assertEqual(len(remote.cancellations), 1)
        report = __import__("json").loads((self.root / "out/execution.json").read_text())
        self.assertIn("cancel unavailable", report["cancellation"])

    def test_terminal_failure_does_not_collect_results(self):
        remote = self.active_remote(d.Batch.State.FAILED)
        with self.assertRaisesRegex(RuntimeError, "FAILED"):
            self.run_remote(remote, sleeper=lambda _: None)
        report = __import__("json").loads((self.root / "out/execution.json").read_text())
        self.assertEqual(report["state"], "FAILED")
        self.assertIn("error", report)
        self.assertFalse((self.root / "out/output_manifest.json").exists())

    def test_cli_uses_detected_identity_when_overrides_are_absent(self):
        from unittest.mock import patch
        from dataclasses import fields
        from scripts.glow.controller import main, CloudOptions
        arguments=[]
        options=cloud()
        for field in fields(CloudOptions):
            if field.type is str and field.name not in ("project","service_account"):
                arguments.extend(["--"+field.name.replace("_","-"),getattr(options,field.name)])
        arguments.extend(["--trait-type","quantitative","--out-dir",str(self.root/"out")])
        for name in ("model-bgen","model-sample","model-variants","sample-map","phenotypes","covariates","phenotype-batches"):
            arguments.extend(["--"+name,str(self.path)])
        listing=self.root/"files.txt"; listing.write_text(str(self.path)+"\n")
        for kind in ("bgen","sample","variants"):
            arguments.extend(["--association-"+kind+"-list",str(listing)])
        values={"project/project-id":options.project,"instance/service-accounts/default/email":options.service_account}
        with patch("scripts.glow.controller.read_instance_metadata",side_effect=values.__getitem__), \
             patch("scripts.glow.controller.run_controller",return_value=[]) as run, \
             patch("scripts.glow.controller.setup_logging"):
            self.assertEqual(main(arguments),0)
        self.assertEqual(run.call_args.args[0],options)
