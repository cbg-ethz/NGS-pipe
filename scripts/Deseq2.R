##################################################
#' Wrapper for DESeq 2
#'
#' @author: Simon Dirmeier
#' @email: simon.dirmeier@bsse.ethz.ch
##################################################

#' Call the DESeq2 tool and test for diff expression
#'
#' @param args  command line arguments
#' @return void



  args<-commandArgs(TRUE)
  in.file <- out.file <- NA
  if (length(args) == 2){
    in.file <- args[[1]]
    out.file <- args[[2]]
  } else{
    stop("\nUSAGE:\n\tRscript deseq.R <Tabular In File> <Output Pattern>")
  }

  library(DESeq2) # load DeSeq2
  dat <- read.table(in.file,header=F)
  expression <- as.matrix(apply(dat[-1,-1], 2, as.double))
  condition <- as.factor(as.vector(unname(unlist(dat[1, -1]))))
  genes <- as.vector(unname(unlist(dat[-1 ,1])))
  dds <- DESeqDataSetFromMatrix(expression, DataFrame(condition), ~ condition)
  diff.exp.analysis <- DESeq(dds)
  diff.exp.analysis.results <- results(diff.exp.analysis)
  sorted.indexes <- order(diff.exp.analysis.results$padj)
  diff.exp.analysis.results <- diff.exp.analysis.results[sorted.indexes,]
  
  ## from here on code is copied from hans
  pdf(paste(out.file,".dispersion.pdf", sep=''))
  plotDispEsts(diff.exp.analysis)
  dev.off()

  print(paste(out.file,".ma.pdf", sep=''))
  pdf(paste(out.file,".ma.pdf", sep=''))
  plotMA(diff.exp.analysis)
  dev.off()

  pdf(paste(out.file,".filtering.pdf", sep=''))
  plot(diff.exp.analysis.results$baseMean, pmin(-log10(diff.exp.analysis.results$pvalue),50), log="x", xlab="mean of normalized counts", ylab=expression(-log[10](pvalue)))
  abline(v=10,col="red",lwd=1)
  dev.off()

  ## end copy code
  
  diff.exp.analysis.results$padj <- p.adjust(diff.exp.analysis.results$pvalue, method="BH")
  res.frame <- as.data.frame(diff.exp.analysis.results)
  sorted.genes <- genes[sorted.indexes]
  res.frame$genes <- sorted.genes
  write.table(res.frame, file=paste(out.file,".results.tsv", sep=''), quote=F, sep="\t", row.names=F)



