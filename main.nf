#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
import groovy.json.JsonSlurper

include { fit_null_glmm } from './modules/fit_null_glmm.nf'
include { run_assoc_tests } from './modules/run_assoc_tests.nf'
include { plot_assoc } from './modules/plot_assoc.nf'

def help_message() {
    log.info """
    A pipeline to perform mixed linear model association analysis
    Author: Vivek Appadurai | vivke.appadurai@regionh.dk

    Usage: nextflow run main.nf

    Options:

    --bfile <genotypes_pruned> [A pruned set of genotypes to construct a GRM in SAIGE] (Default: iPSYCH2012 EurUnrel)
    --genotypes <file.json> [A .json file of per chromosome VCF files to perform association analysis] (Default: iPSYCH2012.json)
    --pheno_cov <file.pheno> [A file of phenotype and covariates] (Default: Phenotype: skizo2015I | Covariates: Age, Gender, first 10 PCs)
    --binary <T/F> [Is the outcome dichotomous?] (Default: T)
    --phenotype <Phenotype Name in the pheno_cov file, will also be used as output prefix> (Default: skizo2015I)
    --help prints this message
    """
}

if(params.help) {
    help_message()
    exit 0
}

log.info """
================================================================================================
IBP - MLM - ASSOCIATION PIPELINE V1.0 - NF
================================================================================================
PLINK genotypes for GRM                  : $params.bfile
VCF genotypes for association tests      : $params.genotypes
File of phenotype to test and covariates : $params.pheno_cov
Is the outcome binary?                   : $params.binary
Name of the phenotype                    : $params.phenotype
================================================================================================
"""

String vcf_files = new File(params.genotypes).text
def vcf_dict     = new JsonSlurper.parseText(vcf_files) 

vcf_geno_ch = Channel.of(1..22) 
    | map {a -> [a, 
        vcf_dict[a.toString()]."vcf", 
        vcf_dict[a.toString()]."idx",
        vcf_dict."meta"."cohort",
        vcf_dict."meta"."population"]

plink_geno_ch = Channel.fromFilePairs(params.bfile + ".{bed,bim,fam}")

workflow {

    // Fit the null GLMM in SAIGE

    plink_geno_ch \
    | combine(params.pheno_cov) \
    | combine(params.phenotype) \
    | combine(Channel.of("Age,Gender,AgeSq,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10")) \
    | combine(Channel.of("IID")) \
    | combine(Channel.of(params.binary)) \
    | combine(Channel.of(target_prefix)) \
    | fit_null_glmm() \
    | set { saige_null_glmm_ch }

    // Perform single variant association tests

    Channel.of(1..22) \
    | combine(vcf_geno_ch, by: 0) \
    | combine(Channel.of("DS")) \
    | combine(Channel.of("0.001")) \
    | combine(Channel.of("1")) \
    | combine(Channel.of(params.sample)) \
    | combine(saige_null_glmm_ch)
    | run_assoc_tets() \
    | collectFile(name: "${target_prefix}_${params.binary}.assoc",
    keepHeader: true,
    skip: 1) \
    | set { assoc_out_ch }

    // Make Manhattan and qq-plots

    assoc_out_ch \
    combine(Channel.of("${target_prefix}_${params.binary}")) \
    plot_assoc()
}