args = commandArgs(trailingOnly=TRUE)
logfile = args[[1]]
workinggds = args[[2]]
outputfile = args[[3]]
idatdir = args[[4]]
sheetname = args[[5]]
tissue = args[[6]]

###this script can used to estimate cell type compositions of buccal and brain data###

set.seed(147)

sink(logfile, append=TRUE, split=TRUE)

print("<<<<<<<<<<<<<<<<<<<<")

print(date())
print("Run cell type estimation:")

library(EpiDISH)

if(tissue == "buccal"){
  print("For buccal cells:")
  data(centEpiFibIC.m)
  } else if(tissue == "saliva"){
  print("For saliva cells:")
  data(centEpiFibIC.m)
  } else if(tissue == "blood"){
  print("For blood cells:")
  data(centDHSbloodDMC.m)
  } else{
  print("Tissue type not found.")
  }

library(bigmelon)

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

gfile <- openfn.gds(workinggds)
betas <- read.gdsn(index.gdsn(gfile, 'betas'))
row.names(betas) <- read.gdsn(index.gdsn(gfile, 'fData/Probe_ID'))
colnames(betas) <- read.gdsn(index.gdsn(gfile, 'pData/barcode'))

if(tissue == "buccal"){
  out.l <- epidish(beta.m = betas, ref.m = centEpiFibIC.m, method = "RPC")
  print("Number of CpGs epithelial cell type estimation is based on (out of 716):")
  print(nrow(out.l$dataREF))
  write.table(out.l$estF, outputfile, sep=",", dec=".", col.names=TRUE, row.names=TRUE)
  } else if(tissue == "saliva"){
  out.l <- epidish(beta.m = betas, ref.m = centEpiFibIC.m, method = "RPC") 
  print("Number of CpGs epithelial cell type estimation is based on (out of 716):")
  print(nrow(out.l$dataREF))
  write.table(out.l$estF, outputfile, sep=",", dec=".", col.names=TRUE, row.names=TRUE)
    
  } else if(tissue == "blood"){
  out.l <- epidish(beta.m = betas, ref.m = centDHSbloodDMC.m, method = "RPC") 
  print("Number of CpGs blood cell type estimation is based on (out of 333):")
  print(nrow(out.l$dataREF))
  write.table(out.l$estF, outputfile, sep=",", dec=".", col.names=TRUE, row.names=TRUE)
    
  } else{
  print("Tissue type not found.")
  write.table("Tissue type not found", outputfile, sep=",", dec=".", col.names=TRUE, row.names=TRUE)
}

closefn.gds(gfile)

print("<<<<<<<<<<<<<<<<<<<<")

sink()
