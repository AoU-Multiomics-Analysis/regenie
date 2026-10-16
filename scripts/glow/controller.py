"""Submit one identifiable Spark batch and collect verified results."""
import argparse
import json
import re
import signal
import time
import uuid
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from google.cloud import dataproc_v1 as dataproc
from google.cloud import storage
from google.api_core.exceptions import AlreadyExists, DeadlineExceeded, NotFound, ServiceUnavailable, InternalServerError
from .contracts import GLOW_COMMIT, LOG, InputArtifact, add_scientific_arguments, options_from_args, require_local_file, scientific_arguments, setup_logging
from .storage import input_identity, stage_inputs, collect_outputs, split_gs

@dataclass(frozen=True)
class CloudOptions:
    project: str
    region: str
    subnetwork_uri: str
    service_account: str
    staging_bucket: str
    staging_prefix: str
    run_id: str
    terra_image: str
    spark_image: str
    driver_cores: int = 4
    driver_memory_mb: int = 8192
    driver_overhead_mb: int = 4096
    executor_cores: int = 4
    executor_memory_mb: int = 8192
    executor_overhead_mb: int = 4096
    min_executors: int = 2
    max_executors: int = 4
    driver_disk_gb: int = 250
    executor_disk_gb: int = 250
    ttl_seconds: int = 21600
    poll_seconds: int = 30

    def __post_init__(self):
        if not re.fullmatch(r"[a-z][a-z0-9-]{4,28}[a-z0-9]",self.project):
            raise ValueError("Invalid Google project ID")
        if not re.fullmatch(r"[a-z]+-[a-z]+[0-9]+",self.region):
            raise ValueError("Invalid Google region")
        if not re.fullmatch(r"[a-z][a-z0-9-]{3,39}",self.run_id):
            raise ValueError("run_id must contain 4-40 lowercase letters, digits, or hyphens")
        subnet=re.fullmatch(r"(?:https://www.googleapis.com/compute/v1/)?projects/[^/]+/regions/([^/]+)/subnetworks/[^/]+",self.subnetwork_uri)
        if not subnet or subnet[1]!=self.region: raise ValueError("Subnetwork region differs from batch region")
        if not re.fullmatch(r"[a-z0-9._-]+@(?:[a-z0-9-]+\.iam|developer)\.gserviceaccount\.com",self.service_account):
            raise ValueError("Invalid workload service account")
        if not re.fullmatch(r"[a-z0-9][a-z0-9._-]{1,220}[a-z0-9]",self.staging_bucket):
            raise ValueError("Invalid staging bucket name")
        split_gs(self.staging_prefix)
        from .storage import validate_relative
        validate_relative(self.staging_prefix.split("/",3)[3].rstrip("/"))
        images=[]
        for image in (self.terra_image,self.spark_image):
            match=re.fullmatch(r"[^\s@]+@sha256:([a-f0-9]{64})",image)
            if not match: raise ValueError("Image references must contain a SHA-256 digest")
            images.append(match[1])
        if images[0]!=images[1]: raise ValueError("Terra and Spark image content digests differ")
        if not re.match(r"(?:[a-z0-9-]+-docker\.pkg\.dev|(?:[a-z0-9-]+\.)?gcr\.io)/",self.spark_image):
            raise ValueError("Spark image must be in a supported Google registry")
        for role in ("driver","executor"):
            cores=getattr(self,role+"_cores")
            memory=getattr(self,role+"_memory_mb"); overhead=getattr(self,role+"_overhead_mb")
            if cores not in (4,8,16) or min(memory,overhead)<=0 or not 1024<= (memory+overhead)/cores <=7424:
                raise ValueError("Invalid Spark cores or memory")
            if getattr(self,role+"_disk_gb")<250: raise ValueError("Spark disk must be at least 250 GiB")
        if not 2<=self.min_executors<=self.max_executors<=2000 or not 600<=self.ttl_seconds<=1209600 or self.poll_seconds<30:
            raise ValueError("Invalid executor bounds, lifetime, or polling interval")

def make_batch(options,staged,science,identity_sha256):
    by_name={item.name:item.uri for item in staged}
    args=[]
    for name,flag in (("model/bgen","model-bgen"),("model/sample","model-sample"),("model/variants","model-variants"),
                      ("sample_map","sample-map"),("phenotypes","phenotypes"),("covariates","covariates"),("phenotype_batches","phenotype-batches")):
        args += ["--"+flag,by_name[name]]
    indices=sorted({int(name.split("/")[1]) for name in by_name if name.startswith("association/")})
    if indices!=list(range(len(indices))) or not indices: raise ValueError("Association sources must be nonempty and contiguous")
    for index in indices:
        for kind in ("bgen","sample","variants"): args += ["--association-"+kind,by_name[f"association/{index}/{kind}"]]
    run_uri=options.staging_prefix.rstrip("/")+"/"+options.run_id
    args += ["--output-prefix",run_uri,"--identity-sha256",identity_sha256]+scientific_arguments(science)
    properties={"spark.dynamicAllocation.enabled":"true","spark.dynamicAllocation.minExecutors":str(options.min_executors),
                "spark.dynamicAllocation.maxExecutors":str(options.max_executors),"spark.executor.instances":str(options.min_executors),
                "spark.pyspark.python":"/opt/conda/bin/python","spark.pyspark.driver.python":"/opt/conda/bin/python"}
    for role in ("driver","executor"):
        properties["spark."+role+".cores"]=str(getattr(options,role+"_cores"))
        properties["spark."+role+".memory"]=str(getattr(options,role+"_memory_mb"))+"m"
        properties["spark."+role+".memoryOverhead"]=str(getattr(options,role+"_overhead_mb"))+"m"
        properties["spark.dataproc."+role+".disk.size"]=str(getattr(options,role+"_disk_gb"))+"g"
    return dataproc.Batch(
        pyspark_batch=dataproc.PySparkBatch(main_python_file_uri=by_name["application/spark_main.py"],args=args,
                                          jar_file_uris=["file:///opt/glow/jars/glow-assembly.jar"]),
        runtime_config=dataproc.RuntimeConfig(version="1.2",container_image=options.spark_image,properties=properties),
        environment_config=dataproc.EnvironmentConfig(execution_config=dataproc.ExecutionConfig(
            service_account=options.service_account,subnetwork_uri=options.subnetwork_uri,staging_bucket=options.staging_bucket,
            ttl={"seconds":options.ttl_seconds})))

def request_matches(actual,expected):
    runtime=actual.runtime_config
    version_ok=runtime.version=="1.2" or runtime.version.startswith("1.2.")
    execution=actual.environment_config.execution_config
    wanted=expected.environment_config.execution_config
    return (version_ok and runtime.container_image==expected.runtime_config.container_image and
            actual.pyspark_batch==expected.pyspark_batch and
            all(runtime.properties.get(k)==v for k,v in expected.runtime_config.properties.items()) and
            execution.service_account==wanted.service_account and execution.subnetwork_uri==wanted.subnetwork_uri and
            execution.staging_bucket==wanted.staging_bucket and execution.ttl==wanted.ttl)

def run_controller(options,science,inputs,out_dir,batch_client=None,storage_client=None,*,clock=time.monotonic,sleeper=time.sleep):
    client=batch_client or dataproc.BatchControllerClient(client_options={"api_endpoint":options.region+"-dataproc.googleapis.com:443"})
    files=storage_client or storage.Client(project=options.project)
    identity=input_identity(inputs,{"science":asdict(science),"cloud":asdict(options),"glow_commit":GLOW_COMMIT})
    run_uri=options.staging_prefix.rstrip("/")+"/"+options.run_id
    staged=stage_inputs(files,run_uri,identity,inputs)
    expected=make_batch(options,staged,science,identity["sha256"])
    parent=f"projects/{options.project}/locations/{options.region}"
    name=parent+"/batches/"+options.run_id+"-glow"
    request_id=str(uuid.uuid5(uuid.NAMESPACE_URL,parent+"/"+options.run_id+"/"+identity["sha256"]))
    out_dir=Path(out_dir); out_dir.mkdir(parents=True,exist_ok=True)
    report={"batch":name,"run_uri":run_uri,"identity_sha256":identity["sha256"],"requested_runtime":"1.2","state":"UNKNOWN"}
    stop=[False]; handlers={}
    def on_signal(signum,frame): stop[0]=True
    for signum in (signal.SIGTERM,signal.SIGINT):
        handlers[signum]=signal.signal(signum,on_signal)
    batch=None
    try:
        for attempt in range(6):
            try:
                batch=client.get_batch(name=name,retry=None,timeout=60)
                break
            except NotFound:
                if stop[0]: raise RuntimeError("Controller stopped before submission")
                try:
                    LOG.info("Submit batch %s",name)
                    client.create_batch(request={"parent":parent,"batch":expected,"batch_id":options.run_id+"-glow","request_id":request_id},retry=None,timeout=60)
                except (AlreadyExists,DeadlineExceeded,ServiceUnavailable,InternalServerError,TimeoutError):
                    LOG.info("Submission outcome uncertain; resolve the same batch ID")
            except (DeadlineExceeded,ServiceUnavailable,InternalServerError):
                LOG.info("Batch lookup will be retried")
            sleeper(min(2**attempt,30))
        if batch is None: raise RuntimeError("Cannot resolve batch after bounded submission retries; inspect "+name)
        deadline=clock()+options.ttl_seconds+1800
        previous=None
        while True:
            if not request_matches(batch,expected): raise ValueError("Existing batch differs from the current run identity/request")
            state=dataproc.Batch.State(batch.state).name
            report.update(state=state,resolved_runtime=batch.runtime_config.version,operation=batch.operation,uuid=batch.uuid)
            (out_dir/"execution.json").write_text(json.dumps(report,indent=2)+"\n")
            if state!=previous:
                LOG.info("Batch %s state %s",name,state); previous=state
            if state=="SUCCEEDED": return collect_outputs(files,run_uri,identity["sha256"],out_dir)
            if state in ("FAILED","CANCELLED"): raise RuntimeError("Remote batch "+state+": "+batch.state_message)
            if state not in ("PENDING","RUNNING","CANCELLING"): raise RuntimeError("Unknown remote batch state "+state)
            if stop[0] or clock()>deadline:
                try:
                    if not batch.operation: raise RuntimeError("Batch has no cancellable operation")
                    client.cancel_operation(request={"name":batch.operation},retry=None,timeout=60)
                    report["cancellation"]="requested"
                except Exception as error:
                    report["cancellation"]="failed: "+str(error)
                raise RuntimeError("Controller stopped or timed out; remote lifetime still applies")
            sleeper(options.poll_seconds)
            for attempt in range(4):
                try:
                    batch=client.get_batch(name=name,retry=None,timeout=60); break
                except (DeadlineExceeded,ServiceUnavailable,InternalServerError):
                    if attempt==3: raise
                    sleeper(min(2**attempt,30))
    except Exception as error:
        report["error"]=str(error)
        raise
    finally:
        (out_dir/"execution.json").write_text(json.dumps(report,indent=2)+"\n")
        for signum,handler in handlers.items(): signal.signal(signum,handler)

def read_file_list(path):
    paths=[Path(line) for line in require_local_file(path).read_text().splitlines()]
    if not paths: raise ValueError("Empty association input list")
    return [require_local_file(p) for p in paths]

def main(argv=None):
    setup_logging(); parser=argparse.ArgumentParser(description=__doc__)
    for field in fields(CloudOptions):
        name="--"+field.name.replace("_","-")
        if field.type is str: parser.add_argument(name,required=True)
        else: parser.add_argument(name,type=int,default=field.default)
    for name in ("model-bgen","model-sample","model-variants","sample-map","phenotypes","covariates","phenotype-batches","out-dir"):
        parser.add_argument("--"+name,type=Path,required=True)
    for name in ("bgen","sample","variants"):
        parser.add_argument("--association-"+name+"-list",type=Path,required=True)
    add_scientific_arguments(parser)
    args=parser.parse_args(argv)
    cloud=CloudOptions(**{f.name:getattr(args,f.name) for f in fields(CloudOptions)})
    artifacts=[InputArtifact(name,require_local_file(getattr(args,flag))) for name,flag in (
        ("model/bgen","model_bgen"),("model/sample","model_sample"),("model/variants","model_variants"),("sample_map","sample_map"),
        ("phenotypes","phenotypes"),("covariates","covariates"),("phenotype_batches","phenotype_batches"))]
    groups=[read_file_list(getattr(args,"association_"+kind+"_list")) for kind in ("bgen","sample","variants")]
    if len({len(group) for group in groups})!=1: raise ValueError("Association input lists differ in length")
    for index,triplet in enumerate(zip(*groups)):
        artifacts.extend(InputArtifact(f"association/{index}/{kind}",path) for kind,path in zip(("bgen","sample","variants"),triplet))
    artifacts.append(InputArtifact("application/spark_main.py",Path(__file__).with_name("spark_main.py")))
    run_controller(cloud,options_from_args(args),artifacts,args.out_dir)
    return 0
if __name__=="__main__":
    raise SystemExit(main())
