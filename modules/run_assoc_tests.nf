#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process run_assoc_tests {
    label 'mod_mem'
    publishDir launchDir

    input:
        tuple val(chr),
            path(vcf),
            path(vcf_idx_tbi),
            path(vcf_idx_csi),
            val(cohort),
            val(population),
            val(snps),
            val(dose_col),
            val(min_maf),
            val(min_mac),
            val(pheno_name),
            path(saige_null_glmm),
            path(saige_variance_ratio),
            path(assoc_test_script_path)

    output:
        path("${cohort}_${population}_${pheno_name}.SAIGE.vcf.dosage.txt")

    script:
    """
    Rscript step2_SPAtests.R \
        --vcfFile=$vcf \
        --vcfFileIndex=$vcf_idx \
        --vcfField=$dose_col \
        --chrom=$chr \
        --minMAF=$min_maf \
        --minMAC=$min_mac \
        --GMMATmodelFile=$saige_null_glmm \
        --varianceRatioFile=$saige_variance_ratio \
        --SAIGEOutputFile=${cohort}_${population}_${pheno_name}.SAIGE.vcf.dosage.txt \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE
    """
}