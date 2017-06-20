#!/opt/local/bin/RScript

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

library(vcfR)
vcf <- read.vcfR(args[1], verbose = FALSE)
chrom <- create.chromR(name="Supercontig", vcf=vcf, verbose=FALSE)
gt <- extract.gt(chrom, element="GT")
dim(gt)

library(ape)
library(gplots)
gt_dist <- dist.gene(t(gt), method = "pairwise")
m_gt_dist <- as.matrix(gt_dist)
pdf(args[2])
heatmap.2(m_gt_dist, trace="none", margins = c(10, 10), main = "Distance between samples based on SNPs", cexRow=0.45, cexCol=0.45 )
dev.off()
