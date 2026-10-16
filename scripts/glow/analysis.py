"""Fit Glow models and validate the chromosome/sample offset boundary."""
from dataclasses import dataclass
import numpy as np
import pandas as pd
from .contracts import LOG

@dataclass
class FittedBatch:
    regression: object
    reduction: object
    sample_ids: list[str]
    traits: list[str]
    trait_type: str
    chromosomes: tuple[str, ...]

def validate_tables(phenotypes,covariates,sample_ids):
    if not phenotypes.index.is_unique or not phenotypes.columns.is_unique:
        raise ValueError("Duplicate sample or trait identifiers")
    if phenotypes.index.tolist()!=sample_ids or covariates.index.tolist()!=sample_ids:
        raise ValueError("Phenotype/covariate sample order differs from genotypes")
    if not np.isfinite(phenotypes.to_numpy(float)).all() or not np.isfinite(covariates.to_numpy(float)).all():
        raise ValueError("Nonfinite analysis table")
    design=np.column_stack([np.ones(len(sample_ids)),covariates.to_numpy(float)])
    if np.linalg.matrix_rank(design)!=design.shape[1]:
        raise ValueError("Covariate design has deficient rank")

def validate_binary_blocks(phenotypes,blocks):
    for held_out in blocks.values():
        train=phenotypes.drop(index=held_out)
        if any(set(train[column])!={0,1} for column in train):
            raise ValueError("Binary training subset lacks a class")

def fit_batch(model_frame,phenotypes,covariates,options):
    from pyspark.sql import functions as F
    from glow.wgr import block_variants_and_samples,RidgeReduction,RidgeRegression,LogisticRidgeRegression
    sample_ids=phenotypes.index.tolist()
    validate_tables(phenotypes,covariates,sample_ids)
    chromosomes=tuple(sorted((row.contigName for row in model_frame.select("contigName").distinct().collect()),key=int))
    if len(chromosomes)<2:
        raise ValueError("Whole-genome model requires at least two chromosomes for LOCO")
    variants=model_frame.withColumn("alternateAlleles",F.array("alternateAllele")).withColumn("names",F.array("variant_id"))
    blocks,mapping=block_variants_and_samples(variants,sample_ids,options.variants_per_block,options.sample_block_count)
    if options.trait_type=="binary":
        validate_binary_blocks(phenotypes,mapping)
    for held_out in mapping.values():
        train=covariates.drop(index=held_out)
        design=np.column_stack([np.ones(len(train)),train.to_numpy(float)])
        if np.linalg.matrix_rank(design)!=design.shape[1]:
            raise ValueError("Training covariate design has deficient rank")
    LOG.info("Fit %d traits using %d sample blocks",len(phenotypes.columns),len(mapping))
    reduction=RidgeReduction(blocks,phenotypes,mapping,cov_df=covariates,label_type=options.trait_type)
    reduction.fit_transform()
    kind=LogisticRidgeRegression if options.trait_type=="binary" else RidgeRegression
    regression=kind.from_ridge_reduction(reduction)
    regression.fit()
    return FittedBatch(regression,reduction,sample_ids,list(phenotypes.columns),options.trait_type,chromosomes)

def predict_chromosome(fit,chromosome):
    if chromosome not in fit.chromosomes:
        raise ValueError("Test chromosome is absent from the model panel")
    kwargs={"chromosomes":[chromosome]}
    if fit.trait_type=="binary": kwargs["response"]="linear"
    result=fit.regression.transform_loco(**kwargs)
    expected=pd.MultiIndex.from_tuples([(identifier,chromosome) for identifier in fit.sample_ids],names=["sample_id","contigName"])
    if not isinstance(result.index,pd.MultiIndex) or set(result.index)!=set(expected) or list(result.columns)!=fit.traits or not result.index.is_unique:
        raise ValueError("LOCO prediction sample, chromosome, or phenotype mismatch")
    result=result.reindex(expected)
    if not np.isfinite(result.to_numpy(float)).all():
        raise ValueError("Nonfinite LOCO prediction")
    return result

def test_chromosome(frame,phenotypes,covariates,offsets,chromosome,options):
    import glow
    validate_tables(phenotypes,covariates,phenotypes.index.tolist())
    expected=[(identifier,chromosome) for identifier in phenotypes.index]
    if not isinstance(offsets.index,pd.MultiIndex) or offsets.index.tolist()!=expected or list(offsets.columns)!=list(phenotypes.columns):
        raise ValueError("LOCO offset sample, chromosome, or phenotype order mismatch")
    if not np.isfinite(offsets.to_numpy(float)).all():
        raise ValueError("Invalid LOCO offset")
    selected=frame.filter(frame.contigName==chromosome)
    kwargs=dict(genotype_df=selected,phenotype_df=phenotypes,covariate_df=covariates,offset_df=offsets,
                contigs=[chromosome],values_column="values",intersect_samples=False)
    if options.trait_type=="binary":
        return glow.gwas.logistic_regression(**kwargs,correction="approx-firth",pvalue_threshold=options.firth_threshold)
    return glow.gwas.linear_regression(**kwargs)

def release_fit(fit):
    for owner in (fit.regression,fit.reduction):
        for name in ("model_df","cv_df","reduced_block_df","block_df"):
            obj=getattr(owner,name,None)
            if obj is not None and hasattr(obj,"unpersist"):
                obj.unpersist()
