#!/usr/bin/env Rscript

require(dplyr, quietly = TRUE)
require(qqman, quietly = TRUE)
require(ggplot2, quietly = TRUE)

args       = commandArgs(trailingOnly = TRUE)
assoc      = read.table(args[1], 
                        header = TRUE, 
                        stringsAsFactors = FALSE)
out_prefix = args[2]

assoc = assoc %>% filter(Is.converge == "TRUE")
assoc = assoc %>% select(CHR, SNPID, POS, p.value)
colnames(assoc) = c("CHR", "SNP", "BP", "P")

png(paste(out_prefix, "Manhattan.png", sep = "_"),
    width = 8, 
    height = 5, 
    units = "in", 
    res = 300)
manhattan(assoc)
dev.off()

png(paste(out_prefix, "QQ.png", sep = "_"), 
    width = 5, 
    height = 5, 
    units = "in", 
    res = 300)
qq(assoc$P)
dev.off()