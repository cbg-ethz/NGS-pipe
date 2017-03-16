#!/cluster/work/bewi/modules/R/3.3.0/bin/Rscript
.libPaths()


# read in command line arguments
args <- commandArgs(TRUE)
bed_file <- args[1] # without first two rows and without last column, i.e. directly starting with the ranges, e.g. chr1 1000 2000
tumor<- args[2]
control <- args[3]
outDir <- args[4]

# set default values for optional parameters
minBaseQ=25
estDispersion=FALSE
overdisp=100
estDirichlet=FALSE
alternative="two.sided"

# now parse over the input arguments to check whether there are optional parameters specified
if(length(args) > 4){
	for(i in 5:length(args)){
		if(args[i]=="--minBaseQ"){
			minBaseQ <- as.numeric(args[i+1])
		} else if (args[i]=="--estimateDispersion"){
			estDispersion <- TRUE
		} else if (args[i]=="--overdispersion"){
			overdisp <- as.numeric(args[i+1])
		} else if (args[i]=="--estimateDirichlet"){
			estDirichlet <- TRUE
		} else if (args[i]=="--alternative"){
			alternative <- args[i+1]
		}
	}
}

# get one string that summarizes all the parameters set - and only specifically mentions them if they are non-default values
##parameter_setting <- paste("minBaseQ",minBaseQ,"_estDispersion",estDispersion,"_overdisp",overdisp,"_estDirichlet",estDirichlet,"_alternative_",alternative,sep="")
# the alternative will always be written into the name of the file, because it is always very important to know
parameter_setting <- paste("alternative_",alternative,sep="")
if (minBaseQ != 25){
	parameter_setting <- paste(parameter_setting,"_minBaseQ",minBaseQ,sep="")
}
if(estDispersion){
	parameter_setting <- paste(parameter_setting,"_estDispersion",estDispersion,sep="")
}
if(overdisp != 100){
	parameter_setting <- paste(parameter_setting,"_overdisp",overdisp,sep="")
}
if(estDirichlet){
	parameter_setting <- paste(parameter_setting,"_estDirichlet",estDirichlet,sep="")
}
cat(paste("parameter_setting = ",parameter_setting,"\n",sep=""))

# specify the output files which will be generated, and create output directory
dir.create(outDir, showWarnings = FALSE)
#summary_deepSNV_object <- paste("", outDir ,"/deepSNV_", sub(".bam","",basename(tumor)) ,"_",sub(".bam","",basename(control)),"_",parameter_setting,".RData", sep="")
#fileToWriteSNVs_withDel <-  paste("",outDir,"/" , sub(".bam","",basename(tumor)),"_",sub(".bam","",basename(control)),"_",parameter_setting,".deepSNV_withDel.txt",sep="")
#fileToWriteSNVs <- paste("",outDir,"/" , sub(".bam","",basename(tumor)),"_",sub(".bam","",basename(control)),"_",parameter_setting,".deepSNV.txt",sep="")
summary_deepSNV_object <- paste("", outDir, "/", sub(".bam","",basename(tumor)) ,"_vs_",sub(".bam","",basename(control)),".RData", sep="")
fileToWriteSNVs_withDel <-  paste("",outDir,"/" , sub(".bam","",basename(tumor)),"_vs_",sub(".bam","",basename(control)),"_withDel.txt",sep="")
fileToWriteSNVs <- paste("",outDir,"/" , sub(".bam","",basename(tumor)),"_vs_",sub(".bam","",basename(control)),".txt",sep="")

if(file.exists(fileToWriteSNVs)){
	cat("Output file already exists. Quitting.\n")
	quit()
}

# load the necessary packages and print the version numbers
library(deepSNV)
library(Rsamtools)
sessionInfo()

# create a GRangesList object from the bed file regions
if ( ! file.exists(bed_file) ){
	stop(paste("Cannot find the bed file ",bed_file," .", sep="")) 
}

if ( ! require("rtracklayer"))
{
	library(BiocInstaller)
	biocLite("rtracklayer") 
	library("rtracklayer")
}

## This command uses the package rtracklayer ###
## It converts the bed file into a GRanges object ###
#all_exons <- import.bed(bed_file,asRangedData=FALSE) # older deepSNV version (1.12.)
all_exons <- import.bed(bed_file) # (deepSNV 1.18.)

## Now, we want to convert this into a GRangesList object ###
## Each chromosome will be its own GRanges object ###
chr_names <- unique(seqnames(all_exons)); 

# fill in the chromosomes 1-22, X, Y
all_exons_list=list()
for(i in 1:length(chr_names)){
	start_index <- (which(as.character((seqnames(all_exons)))==chr_names[i])[1])
	length_chr <- length(which(as.character((seqnames(all_exons)))==chr_names[i]))
	end_index <-  (which(as.character((seqnames(all_exons)))==chr_names[i])[length_chr])
	all_exons_list[[i]] <- all_exons[start_index:end_index]
}
### import.bed automatically took care of that the 0-based bed file is converted into a 1-based GRanges list. So the ranges are correct.

print("Exons loaded.")
print("Working with the following parameters for deepSNV:")
print("bed_file:")
print(bed_file)
print("control bam:")
print(control)
print("tumor bam:")
print(tumor)
print("output directory:")
print(outDir)

nr_chr <- length(all_exons_list)
#__ Run deepSNV ____________________________________________________________________________________________________________________________________________________________________
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
deepSNV_TU_NO=list()
summary_deepSNV_TU_NO_all=list()
print("Number of chromosomes:")
print(nr_chr)
#__ Loop over all chromosomes -- 23=chromosome X, 24=chromosome Y __________________________________________________________________________________________________________________
n=0 # count how many tests where performed
for (i in 1:nr_chr){
	print("i=")
	print(i)
	print(tumor)
	print(control)
	print(alternative)
	print(overdisp)
	print(head(all_exons_list[[i]]))
	print(minBaseQ)
	deepSNV_TU_NO[[i]] = deepSNV(test=tumor, control=control, alternative = alternative, model="betabin", over.dispersion=overdisp, regions=all_exons_list[[i]], q=minBaseQ)
	print("deepSNV done")
        ## alternative=two.sided: Test if the error rate in the tumor is different (greater or smaller) than in the control
	if (estDispersion){
		deepSNV_TU_NO[[i]] <- estimateDispersion(deepSNV_TU_NO[[i]]) # If test is a deepSNV object, automatically taken from the corresponding slot if unspecified
		print("Estimated the dispersion")
	}
	if (estDirichlet){
		deepSNV_TU_NO[[i]] <- estimateDirichlet(deepSNV_TU_NO[[i]])
		print("Estimated the Dirichlet prior")
	}

	# Create the summary, using adjust.method NULL
	print(paste("Create summary for chromosome ",i,".",sep="")) 
	summary_deepSNV_TU_NO_all[[i]]=summary(deepSNV_TU_NO[[i]],sig.level=1-10^(-7),adjust.method=NULL)
	print(paste("Summary ready for chromosome ",i,".",sep=""))
	nr_tests=dim(summary_deepSNV_TU_NO_all[[i]])[1]
	print(paste("Number of tests on chromosome ",i,": ",nr_tests,sep=""))
	n=n+nr_tests
}
print(paste("Computed n as: ",n,sep=""))
## n should be the number of positions we are analyzing times 4 (all the bases and "-", i.e. {A,C,G,T,-} without the consensus base)

## choosing sig.level=0.6 or higher gives an error message with the deepSNV version 1.12.0 or more recent versions and alternative="greater"
## with alternative="two.sided", it works with sig.level=1-10^(-7) which is at least very close to 1,
## the number of SNVs left in the list, i.e. the number of tests performed, affects the correction for multiple testing,
## but we remove only the tests with pvalue of 1. We can argue that it is ok to correct only for the number of tests which are left now, because the other tests where not "testable"
## If we consider for example that VarScan2 only tests if certain conditions are met, it is fine

#__ Save the summary and the deepSNV object as intermediate result ________________________________________________________________________________________________________________
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
save(deepSNV_TU_NO, summary_deepSNV_TU_NO_all, file=summary_deepSNV_object)
#__________________________________________________________________________________________________________________________________________________________________________________

print("Concatenate all the summary tables...")
#__ Concatenate the summary tables from each chromosome ____________________________________________________________________________________________________________________________
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
summary_all = summary_deepSNV_TU_NO_all[[1]]
if(nr_chr > 1){
	for(i in 2:nr_chr){
		summary_all=rbind(summary_all, summary_deepSNV_TU_NO_all[[i]])
	}
}
#___________________________________________________________________________________________________________________________________________________________________________________

print("Adjust for multiple testing.")
#__ Adjust for multiple testing ____________________________________________________________________________________________________________________________________________________
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
adjusted_pValues <- p.adjust(summary_all$p.val, method="BH", n) ## Moritz recommended BH for exome data
#___________________________________________________________________________________________________________________________________________________________________________________


##__ Create final output lists  _____________________________________________________________________________________________________________________________________________
##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
output <- summary_all
output[,16] <- adjusted_pValues
names(output)[5] <- "raw.pVal"
names(output)[16] <- "adjusted.pVal"
##output_hc <- output[which(output[,16]<1),]
output_woDels <- output[-which(output == "-",arr.ind=TRUE)[,1],] ## write also an output file that does not contain any deletions
#___________________________________________________________________________________________________________________________________________________________________________________


#__ Save the final output table with the snvs ______________________________________________________________________________________________________________________________________
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
write.table(output, file = fileToWriteSNVs_withDel,sep = "\t", na  = "NA", dec = ".", row.names = FALSE, col.names = TRUE, quote=FALSE)
write.table(output_woDels, file = fileToWriteSNVs, sep = "\t", na  = "NA", dec = ".", row.names = FALSE, col.names = TRUE, quote=FALSE)
#___________________________________________________________________________________________________________________________________________________________________________________
print("Finished successfully")
