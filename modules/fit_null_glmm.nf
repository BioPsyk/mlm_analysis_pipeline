#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process fit_null_glmm {
    label 'big_mem'
    
    input:
        tuple val(bfile),
            path(bed),
            path(bim),
            path(fam),
            path(pheno_cov),
            val(pheno_name),
            val(cov_cols),
            val(id_col),
            val(outcome),
            val(n_threads),
            val(loco),
            val(inv_normalize_qt),
            path(null_glmm_script)

    output:
        tuple path("${bfile}_${pheno_name}.rda"),
            path("${bfile}_${pheno_name}.varianceRatio.txt")

    script:
    if(outcome == "binary") {
        """
        Rscript step1_fitNULLGLMM.R \
            --plinkFile=${bfile} \
            --phenoFile=${pheno_cov} \
            --phenoCol=${pheno_name} \
            --covarColList=${cov_cols} \
            --sampleIDColinphenoFile=${id_col} \
            --traitType=${outcome} \
            --outputPrefix=${bfile}_${pheno_name} \
            --nThreads=${n_threads} \
            --LOCO=${loco} \
            --IsOverwriteVarianceRatioFile=TRUE
        """
    }
    else {
        """
        Rscript step1_fitNULLGLMM.R \
            --plinkFile=${bfile} \
            --phenoFile=${pheno_cov} \
            --phenoCol=${pheno_name} \
            --covarColList=${cov_cols} \
            --sampleIDColinphenoFile=${id_col} \
            --traitType=${outcome} \
            --invNormalize=${inv_normalize_qt} \
            --outputPrefix=${bfile}_${pheno_name} \
            --nThreads=${n_threads} \
            --LOCO=${loco} \
            --tauInit=1,0
    """
    }
}