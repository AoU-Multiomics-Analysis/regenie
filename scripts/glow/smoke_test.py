"""Run real PLINK conversion and one/multiple-worker Glow smoke tests."""
import argparse
import gzip
import hashlib
import json
import logging
import os
import shutil
import subprocess
import sys
from dataclasses import replace
from pathlib import Path
import numpy as np
import pandas as pd
from .contracts import SPARK_ARROW_PROPERTIES, ScientificOptions, scientific_arguments, setup_logging, LOG
from .prepare import prepare_cohort
from .convert import convert_pgen
from .storage import local_uri

def run(command,log):
    LOG.info("Run %s; evidence %s",Path(str(command[0])).name,log)
    with Path(log).open("w") as stream:
        result=subprocess.run([str(value) for value in command],stdout=stream,stderr=subprocess.STDOUT)
    if result.returncode:
        print(Path(log).read_text()[-18000:],file=sys.stderr)
        raise RuntimeError(f"Command failed with exit {result.returncode}: {log}")

def generate_data(root):
    rng=np.random.default_rng(20261015)
    samples=[f"{index:04d}" for index in range(240)]
    dosages=[]
    for variant in range(200):
        probability=.1+(variant%8)*.1
        values=rng.binomial(2,probability,len(samples)).astype(float)
        values=np.clip(values+.25*((np.arange(len(samples))+variant)%3-1),0,2)
        values[(np.arange(len(samples))+variant)%29==0]=np.nan
        dosages.append(values)
    dosages=np.asarray(dosages)
    dosages[196]=0
    dosages[197]=0; dosages[197,0]=1
    dosages[198]=np.nan
    dosages[199]=1
    vcf=root/"data.vcf"
    with vcf.open("w") as stream:
        stream.write('##fileformat=VCFv4.2\n##contig=<ID=1>\n##contig=<ID=2>\n')
        stream.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        stream.write('##FORMAT=<ID=DS,Number=1,Type=Float,Description="ALT dosage">\n')
        stream.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+"\t".join(samples)+"\n")
        for index,values in enumerate(dosages):
            calls=["./.:." if np.isnan(value) else ("0/0","0/1","1/1")[int(np.rint(value))]+f":{value:.6f}" for value in values]
            stream.write(f"{1+index//100}\t{1000+(index%100)*10}\tv{index+1}\tA\tG\t.\tPASS\t.\tGT:DS\t"+"\t".join(calls)+"\n")
    run(["plink2","--vcf",vcf,"dosage=DS","--double-id","--make-pgen","--out",root/"whole","--threads","2"],root/"import.log")
    psam=root/"whole.psam"; original=pd.read_csv(psam,sep=r"\s+",dtype=str)
    pd.DataFrame({"#IID":original.IID}).to_csv(psam,sep="\t",index=False)
    age=rng.uniform(20,75,len(samples))
    qt1=.4*np.nan_to_num(dosages[0],nan=.5)+.02*age+rng.normal(size=len(samples))
    qt2=.3*np.nan_to_num(dosages[110],nan=.5)+rng.normal(size=len(samples))
    probability=1/(1+np.exp(-(-.5+.3*np.nan_to_num(dosages[5],nan=.5)+.01*age)))
    bt1=rng.binomial(1,probability); bt2=rng.binomial(1,.45,len(samples))
    pheno=pd.DataFrame({"FID":"0","IID":samples,"qt1":qt1,"qt2":qt2,"bt1":bt1,"bt2":bt2,"bt1_12":bt1+1,"bt2_12":bt2+1})
    pheno.to_csv(root/"phenotypes.tsv",sep="\t",index=False)
    pd.DataFrame({"FID":"0","IID":samples,"age":age}).to_csv(root/"covariates.tsv",sep="\t",index=False)
    (root/"common.txt").write_text("".join(f"v{i+1}\n" for i in range(196)))
    np.savez(root/"expected.npz",samples=np.asarray(samples),dosages=dosages)

def spark_command(master):
    workers=1 if master=="local[2]" else int(master.split("[",1)[1].split(",",1)[0])
    arrow_args=[argument for key,value in SPARK_ARROW_PROPERTIES.items() for argument in ("--conf",key+"="+value)]
    return ["spark-submit","--master",master,"--driver-memory","1g","--executor-memory","1g",*arrow_args,
            "--conf","spark.dynamicAllocation.enabled=false","--conf",f"spark.cores.max={workers}",
            "--conf","spark.executor.cores=1","--conf","spark.executor.memoryOverhead=384m",
            "--conf","spark.rpc.netty.dispatcher.numThreads=2","--conf","spark.shuffle.io.serverThreads=2",
            "--conf","spark.shuffle.io.clientThreads=2","--conf","spark.executorEnv.OPENBLAS_NUM_THREADS=1",
            "--conf","spark.executorEnv.OMP_NUM_THREADS=1","--conf","spark.ui.enabled=false",
            "--jars","/opt/glow/jars/glow-assembly.jar"]

def read_associations(prefix,collected=False):
    manifest=json.loads((prefix/("output_manifest.json" if collected else "outputs.json")).read_text())
    frames=[pd.read_csv(prefix/row["relative_path"] if collected else local_uri(row["uri"]),sep="\t") for row in manifest["outputs"] if row["kind"]=="association"]
    frame=pd.concat(frames,ignore_index=True)
    if frame.empty or not frame.groupby("phenotype").pvalue.apply(lambda x:np.isfinite(x).any()).all():
        raise ValueError("No finite association results for a phenotype")
    return frame.sort_values(["source_index","contigName","variant_id","phenotype"]).reset_index(drop=True)

def main(argv=None):
    setup_logging(); parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir",type=Path,required=True)
    parser.add_argument("--spark-master",default=None)
    args=parser.parse_args(argv); root=args.work_dir.resolve(); root.mkdir(parents=True,exist_ok=True)
    logging.getLogger().addHandler(logging.FileHandler(root/"smoke.log"))
    for name in ("build.json", "environment.explicit.txt", "python-packages.txt"):
        shutil.copyfile(Path("/opt/glow") / name, root / name)
    run(["micromamba", "list", "--explicit", "--name", "base"], root / "ci-environment.explicit.txt")
    run(["plink2", "--version"], root / "plink-version.txt")
    run(["java", "-version"], root / "java-version.txt")
    if not os.access(Path(os.environ["SPARK_HOME"]) / "work", os.W_OK):
        raise RuntimeError("CI Spark worker directory must be writable by the image user")
    data=root/"data"; data.mkdir(); generate_data(data)
    base=ScientificOptions("quantitative",phenotype_columns=("qt1","qt2"),sample_block_count=3,variants_per_block=50,batch_size=1)
    profiles=[("quantitative",base),("binary01",replace(base,trait_type="binary",phenotype_columns=("bt2","bt1"),batch_size=2)),
              ("binary12",replace(base,trait_type="binary",phenotype_columns=("bt2_12","bt1_12"),binary_encoding="12",batch_size=2))]
    prepared={}
    for name,options in profiles:
        prepared[name]=prepare_cohort(data/"whole.psam",data/"phenotypes.tsv",data/"covariates.tsv",None,None,options,root/(name+"_prepared"))
    model=convert_pgen(data/"whole.pgen",data/"whole.pvar",data/"whole.psam",prepared["quantitative"].sample_map,
                       data/"common.txt",0,"model",root/"model",2,1024)
    association=convert_pgen(data/"whole.pgen",data/"whole.pvar",data/"whole.psam",prepared["quantitative"].sample_map,
                             None,0,"association",root/"association",2,1024)
    for cohort in prepared.values():
        if cohort.sample_map.read_bytes()!=prepared["quantitative"].sample_map.read_bytes(): raise ValueError("Synthetic sample sets differ")
    reference=root/"reference"
    run(spark_command("local[2]")+["/opt/regenie-spark/scripts/glow/reference_test.py",
        "--bgen",association.bgen,"--sample",association.sample,"--variants",association.variants,
        "--sample-map",prepared["quantitative"].sample_map,"--expected",data/"expected.npz","--out-dir",reference],root/"reference.log")
    results={}
    for name,options in profiles:
        masters=[args.spark_master] if args.spark_master else (["local-cluster[1,1,1536]"] if name=="binary12" else ["local-cluster[1,1,1536]","local-cluster[2,1,1536]"])
        for master in masters:
            workers=master.split("[",1)[1].split(",",1)[0] if master.startswith("local-cluster") else "local"
            output=root/(name+"_workers"+workers)
            cohort=prepared[name]
            command=spark_command(master)+["/opt/regenie-spark/scripts/glow/spark_main.py",
                "--model-bgen",model.bgen.as_uri(),"--model-sample",model.sample.as_uri(),"--model-variants",model.variants.as_uri(),
                "--association-bgen",association.bgen.as_uri(),"--association-sample",association.sample.as_uri(),"--association-variants",association.variants.as_uri(),
                "--sample-map",cohort.sample_map.as_uri(),"--phenotypes",cohort.phenotypes.as_uri(),"--covariates",cohort.covariates.as_uri(),
                "--phenotype-batches",cohort.batches.as_uri(),"--output-prefix",output.as_uri(),
                "--identity-sha256",hashlib.sha256((name+master).encode()).hexdigest()]+scientific_arguments(options)
            run(command,root/(name+"_workers"+workers+".log"))
            result=read_associations(output)
            if name in results:
                pd.testing.assert_frame_equal(results[name],result,check_dtype=False,check_exact=False,rtol=1e-6,atol=1e-8)
            else: results[name]=result
    # Exercise the exact single-VM CLI used by the WDL, including file lists and collection.
    localized=root/"localized ' % [files]"; localized.mkdir()
    localized_model={}; localized_association={}
    for kind in ("bgen","sample","variants"):
        for role,source,paths in (("model",model,localized_model),("association",association,localized_association)):
            target=localized/(role+" ' % [data]."+kind)
            shutil.copyfile(getattr(source,kind),target); paths[kind]=target
    for name,options in profiles:
        cohort=prepared[name]; output=root/(name+"_terra_local")
        command=[sys.executable,"-m","scripts.glow.local_runner","--cpu","2","--memory-gb","4","--driver-memory-gb","1",
                 "--parallelism","2","--shuffle-partitions","2","--out-dir",str(output)]
        for flag,path in (("model-bgen",localized_model["bgen"]),("model-sample",localized_model["sample"]),("model-variants",localized_model["variants"]),
                          ("sample-map",cohort.sample_map),("phenotypes",cohort.phenotypes),("covariates",cohort.covariates),("phenotype-batches",cohort.batches)):
            command += ["--"+flag,str(path)]
        for kind in ("bgen","sample","variants"):
            file_list=root/(name+"_"+kind+".txt"); file_list.write_text(str(localized_association[kind])+"\n")
            command += ["--association-"+kind+"-list",str(file_list)]
        run(command+scientific_arguments(options),root/(name+"_terra_local.log"))
        pd.testing.assert_frame_equal(results[name],read_associations(output,collected=True),check_dtype=False,check_exact=False,rtol=1e-6,atol=1e-8)
    first=results["binary01"].copy(); second=results["binary12"].copy()
    second["phenotype"]=second.phenotype.str.removesuffix("_12")
    pd.testing.assert_frame_equal(first,second,check_dtype=False,check_exact=False,rtol=1e-6,atol=1e-8)
    (root/"summary.json").write_text(json.dumps({"status":"passed","samples":240,"variants":200,
                  "profiles":list(results),"single_vm_cli":True,"localized_special_paths":True,"separate_worker_processes":args.spark_master is None,"rtol":1e-6,"atol":1e-8},indent=2)+"\n")
    LOG.info("Synthetic Spark tests passed")
    return 0
if __name__=="__main__": raise SystemExit(main())
