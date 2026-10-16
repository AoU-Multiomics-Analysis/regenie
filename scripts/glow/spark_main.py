"""Run model fitting and association testing inside a Spark batch."""
import argparse
import gzip
import importlib.metadata
import json
import logging
import os
import shutil
import tempfile
from dataclasses import asdict
from pathlib import Path
import pandas as pd
from .contracts import GLOW_COMMIT, LOG, OutputArtifact, add_scientific_arguments, options_from_args, setup_logging
from .storage import checksum, collect_outputs, list_uris, local_uri, object_hash, object_size, publish_outputs, remote_blob, upload_file

def local_copy(client,uri,path):
    if uri.startswith("file:"): shutil.copyfile(local_uri(uri),path)
    elif uri.startswith("gs:"): remote_blob(client,uri).download_to_filename(str(path))
    else: raise ValueError("Expected an explicitly staged gs:// input or local CI file URI")
    return path

def parse_arguments(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("model-bgen","model-sample","model-variants","sample-map","phenotypes","covariates","phenotype-batches","output-prefix","identity-sha256"):
        parser.add_argument("--"+name,required=True)
    for kind in ("bgen","sample","variants"):
        parser.add_argument("--association-"+kind,action="append",required=True)
    add_scientific_arguments(parser)
    args=parser.parse_args(argv)
    if len({len(args.association_bgen),len(args.association_sample),len(args.association_variants)})!=1:
        raise ValueError("Association input arrays differ in length")
    if len(args.identity_sha256)!=64: raise ValueError("Invalid run identity hash")
    return args

def main(argv=None):
    setup_logging(); args=parse_arguments(argv); options=options_from_args(args)
    import glow
    from pyspark.sql import SparkSession, functions as F
    from google.cloud import storage
    from .analysis import fit_batch,predict_chromosome,test_chromosome,release_fit
    from .genotypes import read_dosages,retain_usable_variants
    client=storage.Client() if args.output_prefix.startswith("gs:") else None
    if client is None and not args.output_prefix.startswith("file:"):
        raise ValueError("Unsupported output prefix")
    spark=glow.register(SparkSession.builder.appName("GlowRegenie").getOrCreate())
    if not spark.version.startswith("3.5.") or importlib.metadata.version("glow.py")!="2.0.3":
        raise ValueError("Expected Spark 3.5 and Glow 2.0.3")
    scala=spark.sparkContext._jvm.scala.util.Properties.versionNumberString()
    if not scala.startswith("2.12."): raise ValueError("Glow JAR requires Scala 2.12")
    receipt=json.loads(Path("/opt/glow/build.json").read_text())
    if receipt["source_commit"]!=GLOW_COMMIT: raise ValueError("Glow source identity mismatch")
    run_uri=args.output_prefix.rstrip("/")
    counts={"sources":[],"batches":[]}
    outputs=[]
    with tempfile.TemporaryDirectory(prefix="glow-driver-") as directory:
        work=Path(directory); log_path=work/"application.log"
        handler=logging.FileHandler(log_path); handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
        logging.getLogger().addHandler(handler)
        try:
            small={}
            for name in ("sample_map","phenotypes","covariates","phenotype_batches"):
                small[name]=local_copy(client,getattr(args,name),work/(name+".tsv"))
            ids=pd.read_csv(small["sample_map"],sep="\t",dtype=str,usecols=["sample_id"]).sample_id.tolist()
            covariates=pd.read_csv(small["covariates"],sep="\t",index_col="sample_id")
            batches=pd.read_csv(small["phenotype_batches"],sep="\t",dtype={"trait_name":str})
            if covariates.index.tolist()!=ids: raise ValueError("Prepared covariate sample order mismatch")
            def load_source(bgen,sample,variants,index,role):
                sample_path=local_copy(client,sample,work/(role+"_sample_"+str(index)+".sample"))
                sample_table=pd.read_csv(sample_path,sep=r"\s+",skiprows=[1],dtype=str)
                if sample_table.ID_2.tolist()!=ids: raise ValueError("BGEN sample file differs from the prepared cohort")
                frame=read_dosages(spark,bgen,sample_path.as_uri(),variants,ids,index)
                frame,report=retain_usable_variants(frame,None if role=="model" else options.min_mac)
                path=run_uri+"/work/"+role+"_"+str(index)
                frame.write.mode("errorifexists").parquet(path)
                counts["sources"].append(dict(report,role=role,source_index=index))
                return spark.read.parquet(path)
            model=load_source(args.model_bgen,args.model_sample,args.model_variants,0,"model")
            associations=[load_source(b,s,v,i,"association") for i,(b,s,v) in enumerate(zip(args.association_bgen,args.association_sample,args.association_variants))]
            chromosomes=sorted({row.contigName for frame in associations for row in frame.select("contigName").distinct().collect()},key=int)
            model_chromosomes={row.contigName for row in model.select("contigName").distinct().collect()}
            if not set(chromosomes)<=model_chromosomes: raise ValueError("Association chromosome is absent from model panel")
            def upload_local(kind,path,basename):
                from .storage import DIRECTORIES
                uri=run_uri+"/results/"+DIRECTORIES[kind]+"/"+basename
                upload_file(client,uri,path)
                outputs.append(OutputArtifact(kind,uri,DIRECTORIES[kind]+"/"+basename,checksum(path),path.stat().st_size))
            for index,table in batches.groupby("batch_index",sort=True):
                traits=table.sort_values("trait_index").trait_name.tolist()
                pheno=pd.read_csv(small["phenotypes"],sep="\t",usecols=["sample_id"]+traits,index_col="sample_id")[traits]
                if pheno.index.tolist()!=ids: raise ValueError("Prepared phenotype sample order mismatch")
                fit=fit_batch(model,pheno,covariates,options)
                try:
                    validation=[json.loads(row) for row in fit.regression.cv_df.toJSON().collect()]
                    counts["batches"].append({"batch_index":int(index),"traits":traits,"samples":len(ids),"model_validation":validation})
                    for chromosome in chromosomes:
                        LOG.info("Predict and test batch %d chromosome %s",index,chromosome)
                        offsets=predict_chromosome(fit,chromosome)
                        prediction=work/"prediction.tsv.gz"
                        offsets.to_csv(prediction,sep="\t",compression="gzip")
                        upload_local("prediction",prediction,f"chr{chromosome}_batch{index}.loco.tsv.gz")
                        for source_index,frame in enumerate(associations):
                            if not frame.filter(F.col("contigName")==chromosome).limit(1).count(): continue
                            result=test_chromosome(frame,pheno,covariates,offsets,chromosome,options)
                            result=result.drop("values").withColumn("sample_count",F.lit(len(ids))).withColumn("effect_allele",F.col("alternateAllele"))
                            prefix=run_uri+f"/results/association_parts/src{source_index}/chr{chromosome}/batch{index}"
                            result.write.mode("errorifexists").option("header",True).option("sep","\t").option("nullValue","NA").option("compression","gzip").csv(prefix)
                            parts=[uri for uri in list_uris(client,prefix) if uri.endswith(".csv.gz")]
                            if not parts: raise ValueError("Association output has no parts")
                            for part_index,uri in enumerate(parts):
                                basename=f"src{source_index}_chr{chromosome}_batch{index}_part{part_index:05d}.tsv.gz"
                                outputs.append(OutputArtifact("association",uri,"associations/"+basename,object_hash(client,uri),object_size(client,uri)))
                        del offsets
                finally: release_fit(fit)
            software={"glow":"2.0.3","spark":spark.version,"scala":scala,"python":os.sys.version.split()[0],
                      "glow_commit":GLOW_COMMIT,"application_commit":os.environ.get("PROJECT_REVISION","unknown")}
            counts["software"]=software; counts["science"]=asdict(options)
            metrics=work/"analysis.json"; metrics.write_text(json.dumps(counts,indent=2,allow_nan=False)+"\n")
            upload_local("metrics",metrics,"analysis.json")
            LOG.info("Analysis finished; verify all outputs")
            handler.flush(); upload_local("log",log_path,"application.log")
            publish_outputs(client,run_uri,args.identity_sha256,outputs,software)
        finally:
            logging.getLogger().removeHandler(handler); handler.close(); spark.stop()
    return 0

if __name__=="__main__":
    raise SystemExit(main())
