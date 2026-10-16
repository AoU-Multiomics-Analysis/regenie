"""Independent small-data checks, executed only in the CI Spark image."""
import argparse
import json
import importlib.metadata
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm
from scripts.glow.genotypes import read_dosages, retain_usable_variants

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("bgen","sample","variants","sample-map","expected","out-dir"):
        parser.add_argument("--"+name,required=True,type=Path)
    args=parser.parse_args(argv)
    import glow
    from pyspark.sql import SparkSession
    assert importlib.metadata.version("glow.py")=="2.0.3"
    assert "/opt/regenie-spark" not in str(Path(glow.__file__))
    spark=glow.register(SparkSession.builder.appName("GlowReference").getOrCreate())
    try:
        mapping=pd.read_csv(args.sample_map,sep="\t",dtype=str)
        expected=np.load(args.expected,allow_pickle=False)
        by_iid={identifier:i for i,identifier in enumerate(expected["samples"])}
        indices=[by_iid[iid] for iid in mapping.IID]
        frame=read_dosages(spark,args.bgen.as_uri(),args.sample.as_uri(),args.variants.as_uri(),mapping.sample_id.tolist(),0)
        for row in frame.select("variant_id","values").collect():
            actual=np.asarray([np.nan if v is None else v for v in row.values],dtype=float)
            desired=expected["dosages"][int(row.variant_id[1:])-1,indices]
            np.testing.assert_allclose(actual,desired,rtol=0,atol=1e-4,equal_nan=True)
        usable,counts=retain_usable_variants(frame,5)
        assert counts["all_missing"]==1,counts
        assert counts["constant"]==2,counts
        assert counts["mac_filtered"]==1,counts
        assert counts["retained"]==196,counts
        n=60
        ids=[f"s{i}" for i in range(n)]
        cov=np.linspace(-1,1,n)
        x=np.asarray([(i*7%19)/10 for i in range(n)])
        offset=.1*np.cos(np.arange(n))
        y=.4*x+.2*cov+np.sin(np.arange(n)*1.7)
        pheno=pd.DataFrame({"qt":y},index=ids)
        covariates=pd.DataFrame({"age":cov},index=ids)
        offsets=pd.DataFrame({"qt":offset},index=ids)
        genotypes=spark.createDataFrame([("1",1,"v",x.tolist())],"contigName string, start long, variant_id string, values array<double>")
        output=glow.gwas.linear_regression(genotypes,pheno,covariates,offset_df=offsets).first()
        fit=sm.OLS(y-offset,np.column_stack([np.ones(n),cov,x])).fit()
        np.testing.assert_allclose([output.effect,output.stderror,output.pvalue],
                                   [fit.params[-1],fit.bse[-1],fit.pvalues[-1]],rtol=1e-6,atol=1e-8)
        binary=np.asarray([(i*13%23)>10 for i in range(n)],dtype=float)
        labels=pd.DataFrame({"bt":binary},index=ids)
        binary_offsets=pd.DataFrame({"bt":offset},index=ids)
        result=glow.gwas.logistic_regression(genotypes,labels,covariates,offset_df=binary_offsets,correction="none").first()
        design=np.column_stack([np.ones(n),cov])
        null=sm.GLM(binary,design,family=sm.families.Binomial(),offset=offset).fit()
        mu=null.predict(); weights=mu*(1-mu)
        residual=x-design@np.linalg.solve(design.T@(weights[:,None]*design),design.T@(weights*x))
        statistic=(residual@(binary-mu))**2/np.sum(weights*residual**2)
        np.testing.assert_allclose([result.chisq,result.pvalue],[statistic,stats.chi2.sf(statistic,1)],rtol=1e-6,atol=1e-8)
        corrected=glow.gwas.logistic_regression(genotypes,labels,covariates,offset_df=binary_offsets,
                                                correction="approx-firth",pvalue_threshold=1.0).first()
        assert np.isfinite(corrected.pvalue) and 0<=corrected.pvalue<=1
        args.out_dir.mkdir(parents=True,exist_ok=True)
        (args.out_dir/"reference.json").write_text(json.dumps({"dosage_tolerance":1e-4,"genotype_counts":counts,
                    "linear_reference":"statsmodels OLS","binary_reference":"independent null score calculation",
                    "numerical_rtol":1e-6,"numerical_atol":1e-8},indent=2)+"\n")
    finally: spark.stop()
    return 0
if __name__=="__main__": raise SystemExit(main())
