"""Build the fixed upstream Glow release on the GitHub image builder."""
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path
from .contracts import GLOW_COMMIT, setup_logging, LOG

def file_hash(path):
    with Path(path).open("rb") as stream: return hashlib.file_digest(stream,"sha256").hexdigest()

def main(argv=None):
    setup_logging(); parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir",type=Path,required=True)
    parser.add_argument("--out-dir",type=Path,required=True)
    parser.add_argument("--sbt-launcher",type=Path,required=True)
    args=parser.parse_args(argv); source=args.source_dir
    actual=subprocess.check_output(["git","rev-parse","HEAD"],cwd=source,text=True).strip()
    if actual!=GLOW_COMMIT: raise ValueError("Glow source commit differs from the approved release")
    if '2.0.3' not in (source/"version.sbt").read_text(): raise ValueError("Unexpected upstream release version")
    args.out_dir.mkdir(parents=True,exist_ok=True)
    LOG.info("Build fixed Glow source %s",GLOW_COMMIT)
    subprocess.run(["java","-XX:ActiveProcessorCount=2","-Xmx3g","-Dsbt.supershell=false","-Dsbt.log.noformat=true",
                    "-jar",str(args.sbt_launcher),"python/generateFunctions","core/assembly"],cwd=source,check=True)
    subprocess.run([sys.executable,"-m","build","--wheel","--no-isolation","--outdir",str(args.out_dir),str(source/"python")],check=True)
    wheels=list(args.out_dir.glob("*.whl"))
    jars=list((source/"core/target/scala-2.12").glob("glow-spark3-assembly-2.0.3.jar"))
    if len(wheels)!=1 or len(jars)!=1: raise ValueError("Expected one Glow wheel and one Scala 2.12 assembly")
    with zipfile.ZipFile(wheels[0]) as wheel:
        names=[n for n in wheel.namelist() if n.endswith(".dist-info/METADATA")]
        if len(names)!=1 or "\nVersion: 2.0.3\n" not in wheel.read(names[0]).decode():
            raise ValueError("Built wheel has the wrong version")
    jar=args.out_dir/"glow-assembly.jar"; shutil.copyfile(jars[0],jar)
    with zipfile.ZipFile(jar) as archive:
        if "org/apache/spark/SparkContext.class" in archive.namelist():
            raise ValueError("Glow assembly must not package the Spark runtime")
    receipt={"source_commit":GLOW_COMMIT,"version":"2.0.3","scala":"2.12.19","spark":"3.5.1","sbt":"1.9.9",
             "jar_sha256":file_hash(jar),"wheel_sha256":file_hash(wheels[0])}
    (args.out_dir/"build.json").write_text(json.dumps(receipt,indent=2)+"\n")
    return 0
if __name__=="__main__": raise SystemExit(main())
