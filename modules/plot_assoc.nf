#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process plot_assoc {
    label 'low_mem'
    publishDir launchDir

    input:
        tuple path(assoc),
            val(out_prefix),
            path(plot_assoc_path)

    output:
        path("${out_prefix}_Manhattan.png")
        path("${out_prefix}_QQ.png")

    script:
    """
    Rscript ./plot_assoc.R $assoc $out_prefix
    """
}