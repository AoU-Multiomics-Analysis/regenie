"""Immutable staging and checked transfer of remote analysis results."""
import hashlib
import json
import os
from dataclasses import asdict
from pathlib import Path, PurePosixPath
from urllib.parse import urlsplit, unquote
from google.api_core.exceptions import PreconditionFailed
from .contracts import InputArtifact, StagedArtifact, OutputArtifact, LOG, require_local_file

DIRECTORIES={"association":"associations","prediction":"predictions","metrics":"metrics","log":"logs"}

def canonical_bytes(value):
    return json.dumps(value,sort_keys=True,separators=(",",":"),allow_nan=False).encode()

def checksum(path):
    with Path(path).open("rb") as stream: return hash_stream(stream)

def hash_stream(stream):
    digest=hashlib.sha256()
    while chunk:=stream.read(8*1024*1024): digest.update(chunk)
    return digest.hexdigest()

def split_gs(uri):
    parts=urlsplit(uri)
    if parts.scheme!="gs" or not parts.netloc or parts.query or parts.fragment or not parts.path.lstrip("/"):
        raise ValueError("Invalid Cloud Storage object URI")
    return parts.netloc,parts.path.lstrip("/")

def local_uri(uri):
    parts=urlsplit(uri)
    if parts.scheme!="file" or parts.netloc or parts.query or parts.fragment:
        raise ValueError("Invalid local test URI")
    return Path(unquote(parts.path))

def spark_path(uri):
    """Hadoop Path(String) expects unescaped local path components."""
    return str(local_uri(uri)) if uri.startswith('file:') else uri

def remote_blob(client,uri):
    bucket,key=split_gs(uri)
    return client.bucket(bucket).blob(key)

def read_bytes(client,uri):
    return local_uri(uri).read_bytes() if uri.startswith("file:") else remote_blob(client,uri).download_as_bytes()

def object_hash(client,uri):
    if uri.startswith("file:"): return checksum(local_uri(uri))
    with remote_blob(client,uri).open("rb") as stream: return hash_stream(stream)

def object_size(client,uri):
    if uri.startswith("file:"): return local_uri(uri).stat().st_size
    blob=remote_blob(client,uri); blob.reload(); return int(blob.size)

def write_new(client,uri,data):
    try:
        if uri.startswith("file:"):
            path=local_uri(uri); path.parent.mkdir(parents=True,exist_ok=True)
            with path.open("xb") as stream: stream.write(data)
        else:
            remote_blob(client,uri).upload_from_string(data,if_generation_match=0)
    except (FileExistsError,PreconditionFailed):
        if read_bytes(client,uri)!=data: raise ValueError("Run object already exists with different content")

def upload_file(client,uri,path):
    if uri.startswith("file:"):
        target=local_uri(uri); target.parent.mkdir(parents=True,exist_ok=True)
        with target.open("xb") as dst,Path(path).open("rb") as src:
            import shutil
            shutil.copyfileobj(src,dst)
    else:
        remote_blob(client,uri).upload_from_filename(str(path),if_generation_match=0)

def input_identity(inputs,settings):
    names=[artifact.name for artifact in inputs]
    if len(set(names))!=len(names): raise ValueError("Duplicate input artifact name")
    records=[]
    for artifact in sorted(inputs,key=lambda item:item.name):
        validate_relative(artifact.name)
        path=require_local_file(artifact.path)
        records.append({"name":artifact.name,"sha256":checksum(path),"size_bytes":path.stat().st_size})
    identity={"schema_version":1,"settings":settings,"inputs":records}
    return dict(identity,sha256=hashlib.sha256(canonical_bytes(identity)).hexdigest())

def stage_inputs(client,run_uri,identity,inputs):
    write_new(client,run_uri.rstrip("/")+"/identity.json",canonical_bytes(identity))
    expected={row["name"]:row for row in identity["inputs"]}
    if set(expected)!={artifact.name for artifact in inputs}: raise ValueError("Identity input names differ")
    staged=[]
    for artifact in inputs:
        record=expected[artifact.name]
        if checksum(artifact.path)!=record["sha256"]: raise ValueError("Input changed after identity calculation")
        uri=run_uri.rstrip("/")+"/inputs/"+artifact.name
        blob=remote_blob(client,uri)
        try: blob.upload_from_filename(str(artifact.path),if_generation_match=0)
        except PreconditionFailed: pass
        blob.reload()
        if int(blob.size)!=record["size_bytes"] or object_hash(client,uri)!=record["sha256"]:
            raise ValueError("Staged input checksum mismatch")
        staged.append(StagedArtifact(artifact.name,uri,record["sha256"],record["size_bytes"],int(blob.generation)))
        LOG.info("Staged input %s (%d bytes)",artifact.name,record["size_bytes"])
    return staged

def validate_relative(name):
    path=PurePosixPath(name)
    if not name or path.is_absolute() or any(part in ("",".","..") for part in name.split("/")) or "\\" in name or "\n" in name or "\r" in name:
        raise ValueError("Unsafe relative output/input path")
    return path

def validate_outputs(run_uri,outputs):
    if not outputs or not {"association","prediction","metrics"} <= {out.kind for out in outputs}:
        raise ValueError("Missing required output categories")
    names=set(); uris=set()
    for out in outputs:
        relative=validate_relative(out.relative_path)
        if out.kind not in DIRECTORIES or relative.parts[0]!=DIRECTORIES[out.kind] or len(relative.parts)!=2:
            raise ValueError("Invalid output kind or destination")
        if out.relative_path in names or out.uri in uris: raise ValueError("Duplicate output destination or URI")
        if not out.uri.startswith(run_uri.rstrip("/")+"/results/"):
            raise ValueError("Output URI is outside the run results prefix")
        if len(out.sha256)!=64 or any(c not in "0123456789abcdef" for c in out.sha256) or out.size_bytes<0:
            raise ValueError("Invalid output checksum or size")
        names.add(out.relative_path); uris.add(out.uri)

def publish_outputs(client,run_uri,identity_sha256,outputs,software):
    validate_outputs(run_uri,outputs)
    for out in outputs:
        if object_size(client,out.uri)!=out.size_bytes or object_hash(client,out.uri)!=out.sha256:
            raise ValueError("Output checksum mismatch before completion")
    manifest=canonical_bytes({"schema_version":1,"identity_sha256":identity_sha256,"software":software,"outputs":[asdict(out) for out in outputs]})
    write_new(client,run_uri.rstrip("/")+"/outputs.json",manifest)
    write_new(client,run_uri.rstrip("/")+"/COMPLETE",hashlib.sha256(manifest).hexdigest().encode())

def collect_outputs(client,run_uri,identity_sha256,out_dir):
    manifest=read_bytes(client,run_uri.rstrip("/")+"/outputs.json")
    if read_bytes(client,run_uri.rstrip("/")+"/COMPLETE").decode()!=hashlib.sha256(manifest).hexdigest():
        raise ValueError("Output completion marker does not match manifest")
    decoded=json.loads(manifest)
    if decoded.get("schema_version")!=1 or decoded.get("identity_sha256")!=identity_sha256:
        raise ValueError("Output identity or manifest version mismatch")
    outputs=[OutputArtifact(**row) for row in decoded["outputs"]]
    validate_outputs(run_uri,outputs)
    out_dir=Path(out_dir); out_dir.mkdir(parents=True,exist_ok=True)
    temp=out_dir/".downloads"; temp.mkdir(exist_ok=True)
    paths=[]
    for index,out in enumerate(outputs):
        target=out_dir/out.relative_path
        if not target.resolve().is_relative_to(out_dir.resolve()): raise ValueError("Output path escapes destination")
        partial=temp/(str(index)+".partial")
        if out.uri.startswith("file:"):
            import shutil
            shutil.copyfile(local_uri(out.uri),partial)
        else: remote_blob(client,out.uri).download_to_filename(str(partial))
        if partial.stat().st_size!=out.size_bytes or checksum(partial)!=out.sha256:
            partial.unlink(missing_ok=True); raise ValueError("Downloaded output checksum mismatch")
        target.parent.mkdir(parents=True,exist_ok=True)
        os.replace(partial,target); paths.append(target)
    (out_dir/"output_manifest.json").write_bytes(manifest)
    LOG.info("Collected %d verified output files",len(paths))
    return paths

def list_uris(client,prefix):
    if prefix.startswith("file:"):
        root=local_uri(prefix)
        return sorted(p.as_uri() for p in root.rglob("*") if p.is_file())
    bucket,key=split_gs(prefix.rstrip("/")+"/")
    return sorted("gs://"+bucket+"/"+blob.name for blob in client.list_blobs(bucket,prefix=key))
