
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

purityFile = paste(args[2],"purity.txt",sep=".")
ploidyFile = paste(args[2],"ploidy.txt",sep=".")

write(fit$purity, file = purityFile)
write(fit$ploidy, file = ploidyFile)

pdf(args[3])
plotSample(x=oo,emfit=fit)
dev.off()

fit2=emcncf(oo)
head(fit2$cncf)

write.table(fit2$cncf, args[2], quote=FALSE, row.names=FALSE, sep="\t")
