#This script and all other scripts and files that are being used for the DNAm age estimation are heavily based on the script from Horvath:
#Horvath, S. DNA methylation age of human tissues and cell types. Genome Biol 14, 3156 (2013). https://doi.org/10.1186/gb-2013-14-10-r115

#Changes were made to adapt the script to the gds-file format for DNAm data and to include it in the DNAm snakemake workflow

args = commandArgs(trailingOnly=TRUE)
logfile = args[[1]]
workinggds = args[[2]]
outputfile = args[[3]]

###this script can be used to predict the sex and epigenetic age of samples based on DNAm data###

set.seed(147)

sink(logfile, append=TRUE, split=TRUE)

print("<<<<<<<<<<<<<<<<<<<<")

print(date())
print("Run Horvath sex and epigenetic age prediction:")

sink()

library("bigmelon")
library("WGCNA")
library("sqldf")
library("impute")

#Normalization step
# Comment regarding the following normalization method based on BMIQ.R
# The original BMIQ function from Teschendorff 2013 (Bioinformatics. 2013 Jan 15;29(2):189-96) 
# adjusts for the type-2 bias in Illumina Infinium 450k data.
# Later functions and edits were provided by yours truly, Steve Horvath.
# I changed the code so that one can calibrate methylation data to a gold standard.
# Specifically, I took version v_1.2 by Teschendorff  and fixed minor issues. 
# Also I made the code more robust e.g. by changing the optimization algorithm.
# Toward this end, I used the method="Nelder-Mead" in optim()

source("scripts/epigeneticage/bmiq.r")

#Comment: The file AdditionalFile24NORMALIZATION.R.txt contains R function which will only be invoked in Step 3 below.
#Age transformation and probe annotation functions
trafo= function(x,adult.age=20) { x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y }
anti.trafo= function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }
probeAnnotation21kdatMethUsed=read.csv("scripts/epigeneticage/13059_2013_3156_MOESM22_ESM.csv")
probeAnnotation27k=read.csv("scripts/epigeneticage/13059_2013_3156_MOESM21_ESM.csv")
datClock=read.csv("scripts/epigeneticage/13059_2013_3156_MOESM23_ESM.csv")
#Read in the DNA methylation data (beta values)
# For a small file, e.g. measured on the 27k platform you could just use read.csv. 
# But for large files, e.g. those measured on the 450K platform, I recommend you use read.csv.sql.
#dat0=read.csv("18062020_lcbc1_rawbeta_unfiltered.csv") ;
#YS: loading data from gds file and adding commands to adapt the dataframe structure to the input expected
gfile <- openfn.gds(workinggds)

dat0 <- read.gdsn(index.gdsn(gfile, "backup/betas"))
colnames(dat0) <- read.gdsn(index.gdsn(gfile, "Sample_ID"))
row.names(dat0) <- read.gdsn(index.gdsn(gfile, "Probe_ID"))

dat0 <- as.data.frame(dat0)
dat0 <- cbind(row.names(dat0), dat0)
row.names(dat0) <- c()
colnames(dat0)[1] <- "ProbeID"


nSamples=dim(dat0)[[2]]-1
nProbes= dim(dat0)[[1]]
# the following command may not be needed. But it is sometimes useful when you use read.csv.sql
dat0[,1]= gsub(x=dat0 [,1],pattern="\"",replacement="") 
#Create a log file which will be output into your directory
# The code looks a bit complicated because it serves to create a log file (for error checks etc).
# It will automatically create a log file.
file.remove("scripts/epigeneticage/LogFile.txt")
file.create("scripts/epigeneticage/LogFile.txt")
DoNotProceed=FALSE
cat(paste( "The methylation data set contains", nSamples, "samples (e.g. arrays) and ", nProbes, " probes."),file="LogFile.txt")
if (nSamples==0) {DoNotProceed=TRUE; cat(paste( "\n ERROR: There must be a data input error since there seem to be no samples.\n Make sure that you input a comma delimited file (.csv file)\n that can be read using the R command read.csv.sql . Samples correspond to columns in that file  ."), file="LogFile.txt",append=TRUE) } 
if (nProbes==0) {DoNotProceed=TRUE; cat(paste( "\n ERROR: There must be a data input error since there seem to be zero probes.\n Make sure that you input a comma delimited file (.csv file)\n that can be read using the R command read.csv.sql  CpGs correspond to rows.")   , file="LogFile.txt",append=TRUE) } 
if (  nSamples > nProbes  ) { cat(paste( "\n MAJOR WARNING: It worries me a lot that there are more samples than CpG probes.\n Make sure that probes correspond to rows and samples to columns.\n I wonder whether you want to first transpose the data and then resubmit them? In any event, I will proceed with the analysis."),file="LogFile.txt",append=TRUE) }
if (  is.numeric(dat0[,1]) ) { DoNotProceed=TRUE; cat(paste( "\n Error: The first column does not seem to contain probe identifiers (cg numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file contains probe identifiers such as cg00000292. Instead it contains ", dat0[1:3,1]  ),file="LogFile.txt",append=TRUE)  } 
if (  !is.character(dat0[,1]) ) {  cat(paste( "\n Major Warning: The first column does not seem to contain probe identifiers (cg numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file contains CpG probe identifiers such as cg00000292. Instead it contains ", dat0[1:3,1]  ),file="LogFile.txt",append=TRUE)  } 
datout=data.frame(Error=c("Input error. Please check the log file for details","Please read the instructions carefully."), Comment=c("", "email Steve Horvath."))
if ( ! DoNotProceed ) {
  nonNumericColumn=rep(FALSE, dim(dat0)[[2]]-1)
  for (i in 2:dim(dat0)[[2]] ){ nonNumericColumn[i-1]=! is.numeric(dat0[,i]) }
  if (  sum(nonNumericColumn) >0 ) { cat(paste( "\n MAJOR WARNING: Possible input error. The following samples contain non-numeric beta values: ", colnames(dat0)[-1][ nonNumericColumn], "\n Hint: Maybe you use the wrong symbols for missing data. Make sure to code missing values as NA in the Excel file. To proceed, I will force the entries into numeric values but make sure this makes sense.\n" ),file="LogFile.txt",append=TRUE)  } 
  XchromosomalCpGs=as.character(probeAnnotation27k$Name[probeAnnotation27k$Chr=="X"])
  selectXchromosome=is.element(dat0[,1], XchromosomalCpGs )
  selectXchromosome[is.na(selectXchromosome)]=FALSE
  meanXchromosome=rep(NA, dim(dat0)[[2]]-1)
  if (   sum(selectXchromosome) >=500 )  {
    meanXchromosome= as.numeric(apply( as.matrix(dat0[selectXchromosome,-1]),2,mean,na.rm=TRUE)) }
  if (  sum(is.na(meanXchromosome)) >0 ) { cat(paste( "\n \n Comment: There are lots of missing values for X chromosomal probes for some of the samples. This is not a problem when it comes to estimating age but I cannot predict the gender of these samples.\n " ),file="LogFile.txt",append=TRUE)  } 
  
  match1=match(probeAnnotation21kdatMethUsed$Name , dat0[,1])
  if  ( sum( is.na(match1))>0 ) { 
    missingProbes= probeAnnotation21kdatMethUsed$Name[!is.element( probeAnnotation21kdatMethUsed$Name , dat0[,1])]    
    #YS: changed DoNotProceed = TRUE to FALSE
    DoNotProceed=FALSE; cat(paste( "\n \n Input error: You forgot to include the following ", length(missingProbes), " CpG probes (or probe names):\n ", paste( missingProbes, sep="",collapse=", ")),file="LogFile.txt",append=TRUE)  } 
  
  #STEP 2: Restrict the data to 21k probes and ensure they are numeric
  match1=match(probeAnnotation21kdatMethUsed$Name , dat0[,1])
  #changed stop to cat, directing the error message to the logfile and continuing the process
  if  ( sum( is.na(match1))>0 ) cat(paste(sum( is.na(match1)), "CpG probes cannot be matched"), file="LogFile.txt",append=TRUE)
  dat1= dat0[match1,]
  asnumeric1=function(x) {as.numeric(as.character(x))}
  dat1[,-1]=apply(as.matrix(dat1[,-1]),2,asnumeric1)
  
  #STEP 3: Create the output file called datout
  set.seed(1)
  # Do you want to normalize the data (recommended)?
  normalizeData=TRUE
  source("scripts/epigeneticage/analysis.r")
  # STEP 4: Output the results 
  if (  sum(  datout$Comment  != "" )   ==0 ) { cat(paste( "\n The individual samples appear to be fine. "),file="LogFile.txt",append=TRUE)  } 
  if (  sum(  datout$Comment != "" )   >0 ) { cat(paste( "\n Warnings were generated for the following samples.\n", datout[,1][datout$Comment != ""], "\n Hint: Check the output file for more details."),file="LogFile.txt",append=TRUE)  } 
} 
# output the results into the directory
write.table(datout, outputfile, row.names=F, sep="," )

closefn.gds(gfile)