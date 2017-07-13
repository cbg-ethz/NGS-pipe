#!/opt/local/bin/RScript

library(ape)
library(gplots)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

df <- read.table(args[1], header = TRUE)
gt_dist <- dist.gene(t(df), method = "pairwise")
m_gt_dist <- as.matrix(gt_dist)
pdf(args[2])
heatmap.2(m_gt_dist, trace="none", margins = c(10, 10), main = "Sample similarity", cexRow=0.45, cexCol=0.45 )
dev.off()
