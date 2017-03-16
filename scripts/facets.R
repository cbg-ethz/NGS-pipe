
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

library(facets)
set.seed(1234)
rcmat = readSnpMatrix(args[1])
xx = preProcSample(rcmat)
oo=procSample(xx,cval=150)
oo$dipLogR
fit=emcncf(oo)
head(fit$cncf)
fit$purity
fit$ploidy

pdf(args[3])
plotSample(x=oo,emfit=fit)
dev.off()

fit2=emcncf2(oo)
head(fit2$cncf)

write.table(fit2$cncf, args[2], quote=FALSE, row.names=FALSE, sep="\t")
