"""Validate BGEN probabilities and retain observed allele counts."""
import math
from .storage import spark_path

def dosage_from_probabilities(values, ploidy=2, phased=False):
    if not values or all(value==0 for value in values):
        return None
    if ploidy != 2 or phased or len(values)!=3:
        raise ValueError("Unsupported phased, nondiploid, or multiallelic genotype")
    if any(not math.isfinite(v) or v<0 or v>1 for v in values) or abs(sum(values)-1)>1e-4:
        raise ValueError("Invalid genotype probabilities")
    return float(values[1]+2*values[2])

def observed_mac(values):
    called=[float(v) for v in values if v is not None]
    if any(not math.isfinite(v) or not 0<=v<=2 for v in called):
        raise ValueError("Invalid genotype dosage")
    total=sum(called)
    return min(total,2*len(called)-total)

def impute_values(values):
    called=[v for v in values if v is not None]
    if not called:
        raise ValueError("No observed genotype dosages")
    observed_mac(values)
    mean=sum(called)/len(called)
    return [mean if v is None else v for v in values]

def read_dosages(spark,bgen_uri,sample_uri,variants_uri,expected_sample_ids,source_index):
    from pyspark.sql import functions as F
    source=spark.read.format("bgen").option("emitHardCalls","false").option("useBgenIndex","false").load(spark_path(bgen_uri))
    if not expected_sample_ids or len(set(expected_sample_ids))!=len(expected_sample_ids):
        raise ValueError("Invalid expected sample IDs")
    first=source.select(F.transform("genotypes",lambda g:g["sampleId"]).alias("ids")).first()
    if first is None or first.ids!=expected_sample_ids:
        raise ValueError("BGEN sample order differs from prepared sample order")
    # Scalar validation is also used to check the reader's probability representation.
    probe=source.select(F.slice("genotypes",1,10).alias("g")).first()
    for genotype in probe.g:
        dosage_from_probabilities(genotype.posteriorProbabilities,genotype.ploidy,genotype.phased)
    def invalid(g):
        probs=g["posteriorProbabilities"]
        total=F.aggregate(probs,F.lit(0.0),lambda acc,p:acc+p)
        called=(F.size(probs)>0)&(total!=0)
        bad_number=F.exists(probs,lambda p:p.isNull()|F.isnan(p)|(p<0)|(p>1))
        return called & ((g["ploidy"]!=2)|g["phased"]|(F.size(probs)!=3)|bad_number|(F.abs(total-1)>1e-4))
    bad=source.filter((F.size("alternateAlleles")!=1)|(F.size("genotypes")!=len(expected_sample_ids))|F.exists("genotypes",invalid))
    if bad.limit(1).count():
        raise ValueError("Unsupported or invalid BGEN genotype record")
    expected_ids=spark.createDataFrame([(expected_sample_ids,)],"expected_ids array<string>")
    ordered=source.crossJoin(F.broadcast(expected_ids))
    if ordered.filter(F.transform("genotypes",lambda g:g["sampleId"])!=F.col("expected_ids")).limit(1).count():
        raise ValueError("BGEN sample order changes between variants")
    def dosage(g):
        p=g["posteriorProbabilities"]
        return F.when((F.size(p)==3)&((p[0]+p[1]+p[2])>0),p[1]+2*p[2]).otherwise(F.lit(None).cast("double"))
    frame=source.select(F.regexp_replace("contigName","^chr","").alias("contigName"),"start",
                        # PLINK writes the PVAR ID as BGEN rsID. Glow returns
                        # names=[SNPID, rsID]; SNPID is empty in this export.
                        F.element_at("names",2).alias("variant_id"),"referenceAllele",
                        F.element_at("alternateAlleles",1).alias("alternateAllele"),
                        F.lit(source_index).alias("source_index"),F.transform("genotypes",dosage).alias("values"))
    metadata=spark.read.option("header",True).option("sep","\t").csv(spark_path(variants_uri)).select(
        "variant_id",F.col("chromosome").alias("expected_chr"),F.col("position").cast("long").alias("expected_pos"),
        F.col("referenceAllele").alias("expected_ref"),F.col("alternateAllele").alias("expected_alt"))
    joined=frame.join(metadata,"variant_id","full")
    mismatched=joined.filter(F.col("expected_chr").isNull()|F.col("contigName").isNull()|
                            (F.col("contigName")!=F.col("expected_chr"))|((F.col("start")+1)!=F.col("expected_pos"))|
                            (F.col("referenceAllele")!=F.col("expected_ref"))|(F.col("alternateAllele")!=F.col("expected_alt")))
    if mismatched.limit(1).count() or frame.groupBy("variant_id").count().filter("count > 1").limit(1).count():
        raise ValueError("BGEN variant metadata differs from exported PVAR")
    return frame.withColumn("observed_count",F.size(F.filter("values",lambda v:v.isNotNull()))).withColumn(
        "alternate_allele_count",F.aggregate("values",F.lit(0.0),lambda acc,v:acc+F.coalesce(v,F.lit(0.0)))).withColumn(
        "minor_allele_count",F.least(F.col("alternate_allele_count"),2*F.col("observed_count")-F.col("alternate_allele_count")))

def retain_usable_variants(frame,min_mac=None):
    from pyspark.sql import functions as F
    total=frame.count()
    called=frame.filter(F.col("observed_count")>0)
    missing=total-called.count()
    usable=called.filter(F.array_min("values")!=F.array_max("values"))
    constant=total-missing-usable.count()
    if min_mac is not None: usable=usable.filter(F.col("minor_allele_count")>=min_mac)
    mean=F.col("alternate_allele_count")/F.col("observed_count")
    usable=usable.withColumn("values",F.transform("values",lambda v:F.coalesce(v,mean)))
    count=usable.count()
    if not count: raise ValueError("No usable variants after genotype filtering")
    return usable,{"total":total,"all_missing":missing,"constant":constant,"retained":count,"mac_filtered":total-missing-constant-count}
