
args <- commandArgs(TRUE)
if(length(args)<3){
	cat(paste("Usage: rank_combination.R <combined_out_file.txt> <tool_1.vcf> <tool_2.vcf> ... <tool_n.vcf>\n",sep=""))
	cat("Command line arguments:\n")
	cat("<combined_out_file.txt>\t\t\t\t= the output file with the combined and ranked variants to be generated\n")
	cat("<tool_1.vcf> <tool_2.vcf> ... <tool_n.vcf>\t= the vcf files of the tools to be combined. Should be at least two files.\n")
	cat("Assumptions:\n")
	cat("- The vcf files have a sixth column with the confidence score from the variant callers for ranking.\n")
	cat("- In the case of MuTect, the seventh column is considered, which has the \"ACCEPT\" or \"REJECT\" label.\n")
	quit()
}


# functions
# normalize the scores to be within 0-250
normalize_scores <- function(raw_scores){
	max_score=max(raw_scores)
	#cat(paste("The maximum score is: ",max_score,"\n",sep=""))
	normalized_scores=(raw_scores/max_score)*250
	normalized_scores
}


## # This will be used to define a unique id chr/pos/ref/var or now rather chr___pos___ref___var, because there was a variant reported from varscan with variant allele "C/G"!
## Assumption: There is no chromosome, position, reference or variant allele that contains the substring "___"
#MutationID_separator="/"
MutationID_separator="___"

outFileName=args[1]
tool1=args[2]
toolsToCombine=list()
MuTectIn=FALSE
cntTools=0
numTools=0
typeoftest=""
for(i in 2:length(args)){
	currTool=args[i]
        # check whether the file is empty, and if so, skip it
        tryCatch({
          test_read <- read.table(args[i])
        
          # if there are lines for input, we continue counting it as an input vcf
	  #cat(paste("Found " ,currTool, " as input.\n",sep=""))
	  numTools=numTools+1
	  #if(grepl("muTect",basename(currTool)) || grepl("mutect",basename(currTool)) || grepl("MuTect",basename(currTool)) || grepl("MUTECT",basename(currTool)) ){
	  if(grepl("muTect",currTool) || grepl("mutect",currTool) || grepl("MuTect",currTool) || grepl("MUTECT",currTool) ){
		MuTectIn=TRUE
		muTectTool=currTool
		#cat("MuTect is among the tools.\n")
		# MuTect does not have a confidence score for ranking -> we consider the 7th column where it says "ACCEPT" or "REJECT"
	  } else {
		cntTools=cntTools+1
		toolsToCombine[[cntTools]]=currTool
		# All tools other than MuTect are assumed to have the confidence score in their 6th column, like in a vcf file
          }
        }, error=function(cond) {
          # the file is probably empty except for the header lines that start with a "#" 
          message("Error: The current vcf could not be read in. Here's the original error message:")
          message(cond)
          message("\nThe current vcf will be skipped.")           
        })
}

ActualNumTools=numTools
if(MuTectIn){
	numTools=numTools-1 # numTools is number of tools without MuTect
}

cat(paste("There were ",ActualNumTools," tools provided for combination:\n",sep=""))
for(i in 1:numTools){
	cat(paste(toolsToCombine[[i]],"\n",sep=""))
}
if(MuTectIn){
	cat(paste(muTectTool,"\n",sep=""))
}


fileToGenerate=paste(sub(".vcf","",sub(".txt","",outFileName)),".txt",sep="")
if(file.exists(fileToGenerate)){
	cat(paste(fileToGenerate, " already exists.\n"))
	quit()
}


filestocombine=unlist(toolsToCombine)
if(MuTectIn){
	filestocombine=c(filestocombine,muTectTool)
}

kk=numTools ## number of tools without MuTect
toolstotest=c(seq(1,numTools)) # for looping over the elements of "filestocombine"

debugginglimit<- -1 # set to negative number to do the entire files
usefulcols<-c(2,6)

orig<-vector("list", ActualNumTools)
kk=numTools
realKK=ActualNumTools

for(toolno in toolstotest){ # we take the results from the different tools
  score <- read.table(filestocombine[toolno],nrows=debugginglimit)
  scoresimp<-score[,usefulcols] # select the useful columns
  scoresimp[,1]<-paste(score[,1],score[,2],score[,4],score[,5],sep=MutationID_separator) # define a unique id chr/pos/ref/var or now rather chr___pos___ref___var because there was a variant reported from varscan with variant allele "C/G"! 
  scoresimp[,2]<-rank(-score[,usefulcols[2]]) # replace the scores by a rank
  names(scoresimp)<-c("id","rank") # rename the columns
  orig[[toolno]]<-scoresimp # store the result
  cat(paste("tool ",toolno," read in.\n",sep=""))
}

if(MuTectIn){
	toolno<-ActualNumTools # muTect will be the last tool in the vector
	usefulcols<-c(2,7) # column seven has the accept/reject value
	score<-read.table(filestocombine[toolno],nrows=debugginglimit)
	scoresimp<-score[,usefulcols] # select the useful columns
	scoresimp[,1]<-paste(score[,1],score[,2],score[,4],score[,5],sep=MutationID_separator) # define a unique id
	scoresimp[,2]<-rank(score[,usefulcols[2]]=="REJECT") # replace the scores by a rank
	names(scoresimp)<-c("id","rank") # rename the columns
	orig[[toolno]]<-scoresimp # store the result
	cat(paste("tool ",toolno," read in.\n",sep=""))
}


if(MuTectIn){
	comby<-orig[ActualNumTools]
	for(toolno in (ActualNumTools-1):1){
		comby<-merge(comby,orig[[toolstotest[toolno]]],by="id",all=TRUE,suffixes=c("",toolno)) # join data frames
		cat(paste("tool ",toolstotest[toolno]," merged.\n",sep=""))
	}
} else {
	comby<-orig[toolstotest[1]]
	for(toolno in 2:kk){
		comby<-merge(comby,orig[[toolstotest[toolno]]],by="id",all=TRUE,suffixes=c("",toolno)) # join data frames
		cat(paste("tool ",toolstotest[toolno]," merged.\n",sep=""))
	}
}

NN<-nrow(comby)+1
comby[is.na(comby)] <- NN # set NAs to end of the list

cat("data combined\n")

# transform to a cube

rescalecube<-comby
for(i in 1:(realKK)){ # transform according to equation in notes
  rescalecube[,i+1]<-comby[,i+1]/(2*(NN))
}

cat("transformed to cube\n")

# transform to a normal
rescaley<-rescalecube
for(i in 1:(realKK)){ # treat the p-type values as the CDF of a normal 
  rescaley[,i+1]<-qnorm(rescalecube[,i+1])
}
cat("transformed to normal\n")

# correlations and covariances
covcororig<-cov.wt(rescaley[,1:(realKK)+1],cor=TRUE,center=FALSE) # the center=FALSE option creates the matrix of second moments
covmat<-covcororig$cov
cormat<-covcororig$cor

cat("The original matrix of second moments:\n")
print(covmat)

decomp<-svd(covmat) # provides a
Ddiag <- decomp$d # diagonal vector
O <- decomp$u # and an ortogonal matrix

newmat<-cormat
newmat[,]<-(sum(cormat)/(realKK)-1)/(realKK-1)
diag(newmat)<-1
desiredcovmat<-t(t(newmat*sqrt(diag(covmat)))*sqrt(diag(covmat)))

cat("A matrix of second moments with equal correlations:\n")
print(desiredcovmat)

decompd<-svd(desiredcovmat) # provides a
Ddiagd <- decompd$d # diagonal vector
Od <- decompd$u # and an ortogonal matrix
sqrtcovmat<-O%*%diag(sqrt(Ddiag))%*%t(O)
sqrtdesiredcovmat<-Od%*%diag(sqrt(Ddiagd))%*%t(Od)
transmat<-solve(sqrtcovmat)%*%sqrtdesiredcovmat

cat("A possible transformation matrix:\n")
print(transmat)# a transformation matrix with balanced correlations

# the problem is that the direction of each column of O is undefined
centreofmass<-rep(0,realKK) # find the center of mass
for(i in 1:(realKK)){
  centreofmass[i]<-mean(rescaley[,i+1])
}

transmatnew<-transmat
signtest<-t(as.matrix(centreofmass))%*%transmat #check whether the new center of mass is in the lower quadrant
for(i in 1:(realKK)){
  if(signtest[i]>0){
    transmatnew[,i]<--transmat[,i] # if not swap it
  }
}
cat("Check that the transformed mean is in the lower quadrant\n")
print(t(as.matrix(centreofmass))%*%transmatnew)

decorry<-as.matrix(rescaley[,1:(realKK)+1])%*%transmatnew # this removes the correlation 

cat("data decorrelated\n")

cat("new correlations\n")
newcovs<-cov.wt(decorry,center=FALSE,cor=TRUE)
print(newcovs)

#save(covcororig,newcovs,file=paste(sub(".vcf","",sub(".txt","",outFileName)),"_covariances.rData",sep=""))

# going back to the cube space

decorcube<-comby
for(i in 1:(realKK)){ # now we go back to the CDF of a normal
  decorcube[,i+1]<-pnorm(decorry[,i])
}

cat("transformed to cube\n")

sumrank<-rank(-rowSums(decorcube[,-1])) # this is the sum of the transformed scores then ranked
prodrank<-rank(-rowSums(log(decorcube[,-1]))) # this is effectively the product of the transformed scores then ranked

decorcube$prodrank<-prodrank

#print(dim(decorcube))

if(ActualNumTools==3){
	# decorcube is a matrix with ActualNumTools+2 columns
	toprintasVCFfile<-cbind(decorcube,decorcube[,2:5]) # just a matrix with 9 columns, because we want the VCF format

} else if(ActualNumTools==2){
	# decorcube is a matrix with ActualNumTools+2 columns
	toprintasVCFfile<-cbind(decorcube,decorcube[,2:4],decorcube[,3:4]) # just a matrix with 9 columns, because we want the VCF format
	
} else if(ActualNumTools==4){
	# decorcube is a matrix with ActualNumTools+2 columns
	toprintasVCFfile<-cbind(decorcube,decorcube[,3:5]) # just a matrix with 9 columns, because we want the VCF format
	
} else if(ActualNumTools==5){
	# decorcube is a matrix with ActualNumTools+2 columns
	toprintasVCFfile<-cbind(decorcube,decorcube[,4:5]) # just a matrix with 9 columns, because we want the VCF format
	
} else if(ActualNumTools==6){
	# decorcube is a matrix with ActualNumTools+2 columns
	toprintasVCFfile<-cbind(decorcube,decorcube[,3]) # just a matrix with 9 columns, because we want the VCF format
	
} else if(ActualNumTools==7){
	# decorcube is a matrix with ActualNumTools+2 columns
	toprintasVCFfile<-cbind(decorcube) # just a matrix with 9 columns, because we want the VCF format
	
} else if(ActualNumTools==8){
	# decorcube is a matrix with ActualNumTools+2 columns
	toprintasVCFfile<-decorcube[,-2] # just a matrix with 9 columns, because we want the VCF format
} else if(ActualNumTools==9){
	# decorcube is a matrix with ActualNumTools+2 columns
	toprintasVCFfile<-decorcube[,c(-2,-3)] # just a matrix with 9 columns, because we want the VCF format
} else {
	stop("Error: This number of tools is not implemented!")
}

names(toprintasVCFfile)<-c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT") # rename the columns
toprintasVCFfile[,c(3,7:9)]<-"." ## now we make the vcf format by setting the 3rd, 7th,8th, and 9th column to "."
toprintasVCFfile[,6]<-normalize_scores(sumrank)

toprintasVCFfile[,c(1:2,4:5)]<-matrix(unlist(strsplit(comby[,1], MutationID_separator)),ncol=4,byrow=TRUE) ## the id is of the format "chr10/100003638/T/C", and here we split it up and assign the different fields to column 1, 2, 4 and 5
names(toprintasVCFfile)<-c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT") # rename the columns

#write.table(toprintasVCFfile,file=paste(sub(".vcf","",sub(".txt","",outFileName)),"_sum.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
#cat("sum data saved\n")

toprintasVCFfile[,6]<-normalize_scores(prodrank)
write.table(toprintasVCFfile,file=fileToGenerate,quote=FALSE,sep="\t",row.names=FALSE)

cat("prod combination saved.\n")


