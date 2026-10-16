"""Export a validated PGEN source to dosage-preserving BGEN."""
import argparse
import io
import json
import subprocess
from dataclasses import dataclass
from pathlib import Path
import pandas as pd
from .contracts import LOG, require_local_file, setup_logging
from .prepare import read_ids

@dataclass(frozen=True)
class ConvertedFiles:
    bgen: Path
    sample: Path
    variants: Path
    report: Path
    log: Path

def canonical_chromosome(value):
    text=str(value).removeprefix("chr")
    if not text.isdigit() or not 1 <= int(text) <= 22:
        raise ValueError("Only autosomes 1-22 are supported")
    return str(int(text))

def read_variants(path):
    lines=[line for line in require_local_file(path).read_text().splitlines() if not line.startswith("##")]
    table=pd.read_csv(io.StringIO("\n".join(lines)),sep=r"\s+",dtype=str,keep_default_na=False)
    table=table.rename(columns={"#CHROM":"chromosome","CHROM":"chromosome","POS":"position","ID":"variant_id","REF":"referenceAllele","ALT":"alternateAllele"})
    needed=["chromosome","position","variant_id","referenceAllele","alternateAllele"]
    if not set(needed)<=set(table.columns) or table.empty:
        raise ValueError("Invalid or empty PVAR")
    if table.variant_id.duplicated().any() or table.variant_id.isin(["", "."]).any():
        raise ValueError("Missing or duplicate variant ID")
    table["chromosome"]=table.chromosome.map(canonical_chromosome)
    for column in ("referenceAllele","alternateAllele"):
        if table[column].str.contains(",").any() or table[column].isin(["", "."]).any():
            raise ValueError("Only biallelic variants with defined alleles are supported")
    if (pd.to_numeric(table.position,errors="raise")<1).any():
        raise ValueError("Invalid variant position")
    return table[needed]

def convert_pgen(pgen,pvar,psam,sample_map,variant_ids,source_index,role,out_dir,threads,memory_mb):
    for path in (pgen,pvar,psam,sample_map):
        require_local_file(path)
    if role not in ("model","association") or source_index<0 or threads<1 or memory_mb<1:
        raise ValueError("Invalid converter settings")
    if role=="model" and variant_ids is None:
        raise ValueError("Model conversion requires a variant list")
    mapping=pd.read_csv(sample_map,sep="\t",dtype=str,keep_default_na=False)
    if mapping.empty or mapping.sample_id.duplicated().any() or mapping.duplicated(["FID","IID"]).any():
        raise ValueError("Invalid sample map")
    samples=read_ids(psam)
    wanted=list(zip(mapping.FID,mapping.IID))
    if not set(wanted)<=set(samples.index):
        raise ValueError("Genotype source lacks a retained sample")
    variants=read_variants(pvar)
    selected=None
    if variant_ids is not None:
        selected={line.strip() for line in require_local_file(variant_ids).read_text().splitlines() if line.strip()}
        variants=variants.loc[variants.variant_id.isin(selected)]
        if not selected or variants.empty or (role=="model" and set(variants.variant_id)!=selected):
            raise ValueError("Variant selection is empty or model IDs are absent")
    out_dir=Path(out_dir); out_dir.mkdir(parents=True,exist_ok=True)
    inputs=out_dir/"inputs"; inputs.mkdir(exist_ok=True)
    prefix=inputs/"source"
    for suffix,path in (("pgen",pgen),("pvar",pvar)):
        target=prefix.with_suffix("."+suffix)
        if not target.exists():
            target.symlink_to(Path(path).resolve())
    ids=dict(zip(wanted,mapping.sample_id))
    local=samples.reset_index(drop=True).copy()
    local["IID"]=[ids.get(pair,"excluded_"+str(i)) for i,pair in enumerate(samples.index)]
    local["FID"]="0"
    if "SEX" not in local: local["SEX"]="0"
    local=local[["FID","IID"]+[c for c in local if c not in ("FID","IID")]]
    local=local.rename(columns={"FID":"#FID"})
    local.to_csv(prefix.with_suffix(".psam"),sep="\t",index=False)
    order=inputs/"sample_order.txt"
    order.write_text("".join("0 "+token+"\n" for token in mapping.sample_id))
    export=out_dir/"genotypes"
    command=["plink2","--pfile",str(prefix),"--keep",str(order),"--indiv-sort","file",str(order),
             "--threads",str(threads),"--memory",str(memory_mb)]
    if variant_ids is not None: command += ["--extract",str(Path(variant_ids).resolve())]
    command += ["--export","bgen-1.2","bits=16","ref-first","id-paste=iid","--out",str(export)]
    log=out_dir/"conversion.log"
    LOG.info("Export %s source %d: %d samples, %d variants",role,source_index,len(mapping),len(variants))
    with log.open("w") as stream:
        subprocess.run(command,check=True,stdout=stream,stderr=subprocess.STDOUT)
    sample=export.with_suffix(".sample")
    actual=pd.read_csv(require_local_file(sample),sep=r"\s+",skiprows=[1],dtype=str)
    if actual.ID_2.tolist()!=mapping.sample_id.tolist():
        raise ValueError("Exported sample order differs from the prepared cohort")
    variants["source_index"]=source_index
    variant_file=out_dir/"variants.tsv"; variants.to_csv(variant_file,sep="\t",index=False)
    report=out_dir/"conversion.json"
    report.write_text(json.dumps({"source_index":source_index,"role":role,"samples":len(mapping),"variants":len(variants),"chromosomes":sorted(set(variants.chromosome),key=int)},indent=2)+"\n")
    return ConvertedFiles(require_local_file(export.with_suffix(".bgen")),sample,variant_file,report,log)

def main(argv=None):
    setup_logging(); parser=argparse.ArgumentParser(description=__doc__)
    for name in ("pgen","pvar","psam","sample-map","out-dir"):
        parser.add_argument("--"+name,type=Path,required=True)
    parser.add_argument("--variant-ids",type=Path)
    parser.add_argument("--role",choices=("model","association"),required=True)
    parser.add_argument("--source-index",type=int,default=0)
    parser.add_argument("--threads",type=int,default=4)
    parser.add_argument("--memory-mb",type=int,default=12000)
    args=parser.parse_args(argv)
    convert_pgen(args.pgen,args.pvar,args.psam,args.sample_map,args.variant_ids,args.source_index,args.role,args.out_dir,args.threads,args.memory_mb)
    return 0
if __name__=="__main__":
    raise SystemExit(main())
