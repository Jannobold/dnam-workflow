args = commandArgs(trailingOnly=TRUE)
logfile = args[[1]]
idatdir = args[[2]]
outputfile = args[[3]]
sheetname = args[[4]]

###this script can used to estimate cell type compositions of buccal and brain data###

set.seed(147)

sink(logfile, append=TRUE, split=TRUE)

print("<<<<<<<<<<<<<<<<<<<<")

print(date())
print("Run cell type estimation of brain samples:")

library("minfi")

#create samplesheet
files <- list.files(idatdir)
files <- files[seq(from=1,to=length(files)-1,by=2)]
Sample_Name <- sub("_Grn.idat", "", files)
Sentrix_ID <- sub("_R.*_Grn.idat", "", files)

for(i in 1:length(Sentrix_ID)){
	if (nchar(Sentrix_ID[i]) > 12){
		Sentrix_ID[i] <- sub(".*_", "", Sentrix_ID[i])
		Sentrix_ID[i] <- sub(".*/", "", Sentrix_ID[i])
	}
}

Sentrix_Position <- sub(".*_", "", Sample_Name)
samplesheet <- data.frame(Sample_Name, Sentrix_ID, Sentrix_Position)

write.table(samplesheet, sheetname, sep=",", dec=".", col.names=TRUE, row.names=FALSE)

sink(stdout(), type = "message")
targets <- read.metharray.sheet(idatdir)
RGset <- read.metharray.exp(targets = targets)
RGset@annotation=c(array='IlluminaHumanMethylationEPIC', annotation='ilm10b2.hg19')

cellcounts <- estimateCellCounts(RGset, compositeCellType = "DLPFC", processMethod = "auto", probeSelect = "auto", cellTypes = c( "NeuN_neg", "NeuN_pos"), referencePlatform = "IlluminaHumanMethylation450k", returnAll = TRUE, meanPlot = FALSE, verbose = TRUE)

write.table(cellcounts$counts, outputfile, sep=",", dec=".", col.names=TRUE, row.names=TRUE)

print("<<<<<<<<<<<<<<<<<<<<")

sink()