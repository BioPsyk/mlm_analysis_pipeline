#!/usr/bin/env Rscript

require(dplyr, quietly = TRUE)
require(qqman, quietly = TRUE)
require(ggplot2, quietly = TRUE)
require(data.table, quietly = TRUE)

args = commandArgs(trailingOnly = TRUE)
assoc = fread(args[1], header = TRUE) 
out_prefix = args[2]

if("Is.SPA.converge" %in% colnames(assoc)) {
    assoc = assoc %>% filter(Is.SPA.converge == 1)
}

if("Is.converge" %in% colnames(assoc)) {
    assoc = assoc %>% filter(Is.converge == 1)
}
assoc = assoc %>% select(CHR, SNPID, POS, `p.value`)
colnanes(assoc) = c("CHR", "SNP", "BP", "P")

png(paste(out_prefix, "Manhattan.png", sep = "_"),
    width = 8, 
    height = 5, 
    units = "in", 
    res = 300)
manhattan(assoc, annotatePval = 5e-8)
dev.off()

png(paste(out_prefix, "QQ.png", sep = "_"), 
    width = 5, 
    height = 5, 
    units = "in", 
    res = 300)
qq(assoc$P, cex = 1, pch = 18)
dev.off()