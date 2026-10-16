import hashlib
import io
import tempfile
import unittest
from pathlib import Path
from google.api_core.exceptions import PreconditionFailed, NotFound
from scripts.glow.contracts import InputArtifact, OutputArtifact
from scripts.glow.storage import input_identity, stage_inputs, publish_outputs, collect_outputs

class Blob:
    def __init__(self,bucket,key): self.bucket=bucket; self.key=key; self.metadata={}
    @property
    def size(self): return len(self.bucket.data[self.key][0])
    @property
    def generation(self): return self.bucket.data[self.key][1]
    def exists(self): return self.key in self.bucket.data
    def reload(self):
        if not self.exists(): raise NotFound("missing")
    def upload_from_string(self,data,if_generation_match=None,**kw):
        if if_generation_match==0 and self.exists(): raise PreconditionFailed("exists")
        if isinstance(data,str): data=data.encode()
        self.bucket.data[self.key]=(data,self.bucket.data.get(self.key,(None,0))[1]+1)
    def upload_from_filename(self,path,**kw): self.upload_from_string(Path(path).read_bytes(),**kw)
    def download_as_bytes(self):
        if not self.exists(): raise NotFound("missing")
        return self.bucket.data[self.key][0]
    def download_to_filename(self,path): Path(path).write_bytes(self.download_as_bytes())
    def open(self,mode): return io.BytesIO(self.download_as_bytes())

class Bucket:
    def __init__(self): self.data={}
    def blob(self,key): return Blob(self,key)
class Client:
    def __init__(self): self.buckets={}
    def bucket(self,name): return self.buckets.setdefault(name,Bucket())

class StorageTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name); self.client=Client(); self.uri="gs://bucket/base/run"
        self.path=self.root/"input"; self.path.write_bytes(b"input data")
    def test_staging_is_immutable_and_independent_of_local_path(self):
        identity=input_identity([InputArtifact("model/bgen",self.path)],{"trait_type":"binary"})
        staged=stage_inputs(self.client,self.uri,identity,[InputArtifact("model/bgen",self.path)])
        moved=self.root/"moved"; moved.write_bytes(self.path.read_bytes())
        self.assertEqual(identity,input_identity([InputArtifact("model/bgen",moved)],{"trait_type":"binary"}))
        self.assertEqual(staged,stage_inputs(self.client,self.uri,identity,[InputArtifact("model/bgen",moved)]))
        other=input_identity([InputArtifact("model/bgen",self.path)],{"trait_type":"quantitative"})
        with self.assertRaises(ValueError): stage_inputs(self.client,self.uri,other,[InputArtifact("model/bgen",self.path)])

    def outputs(self):
        outputs=[]
        for kind,directory,name in (("association","associations","src0_chr1_batch0.tsv.gz"),("prediction","predictions","chr1_batch0.tsv.gz"),("metrics","metrics","counts.json")):
            uri=self.uri+"/results/"+name
            self.client.bucket("bucket").blob("base/run/results/"+name).upload_from_string(b"result")
            outputs.append(OutputArtifact(kind,uri,directory+"/"+name,hashlib.sha256(b"result").hexdigest(),6))
        return outputs

    def test_complete_outputs_are_verified_and_collected(self):
        outputs=self.outputs()
        with self.assertRaises(Exception): collect_outputs(self.client,self.uri,"identity",self.root/"out")
        publish_outputs(self.client,self.uri,"identity",outputs,{"glow":"2.0.3"})
        paths=collect_outputs(self.client,self.uri,"identity",self.root/"out")
        self.assertEqual(len(paths),3)
        self.assertTrue(all(p.read_bytes()==b"result" for p in paths))
        self.client.bucket("bucket").data["base/run/results/counts.json"]=(b"bad",2)
        with self.assertRaises(ValueError): collect_outputs(self.client,self.uri,"identity",self.root/"out2")

    def test_path_escape_and_neighbor_prefix_are_rejected(self):
        from dataclasses import replace
        for changes in ({"relative_path":"../escape"},{"relative_path":"/escape"},{"uri":"gs://other/base/run/data"},{"uri":self.uri+"2/data"}):
            outputs=self.outputs(); outputs[0]=replace(outputs[0],**changes)
            with self.assertRaises(ValueError): publish_outputs(self.client,self.uri,"identity",outputs,{})

    def test_duplicate_output_paths_fail(self):
        outputs=self.outputs()
        with self.assertRaises(ValueError): publish_outputs(self.client,self.uri,"identity",outputs+[outputs[0]],{})
