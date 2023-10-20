library("ewastools")
library("bigmelon")
setwd("PATH TO PROJECT DIRECTORY")

set.seed(147)

gfile <- openfn.gds("gdsfiles/PROJECTNAME_raw_filtered.gds")
beta <- read.gdsn(index.gdsn(gfile, 'backup/betas'))
colnames(beta) <- read.gdsn(index.gdsn(gfile, 'Sample_ID'))
row.names(beta) <- read.gdsn(index.gdsn(gfile, 'Probe_ID'))

closefn.gds(gfile)

beta <- tail(beta, n=59)

#this file contains the raw data for the SNP probes from the EPIC array
write.table(beta, "OUTPUTFILENAME.csv", sep=",", dec=".", col.names=TRUE, row.names=TRUE)

#PATH TO MANIFEST is where you copy-pasted / uploaded the EPIC_SNPprobes_manifest.csv file (from this repository)
snps <- read.csv("PATH TO MANIFEST/EPIC_SNPprobes_manifest.csv", sep=",", dec=".", header=TRUE)

#The following code determines the genotypes based on the raw data extracted above
beta <- as.matrix(beta)

genotypes <- call_genotypes(beta, learn = FALSE, maxiter = 50)
genotypes <- genotypes$gamma

low <- genotypes[[1]]
medium <- genotypes[[2]]
high <- genotypes[[3]]

colnames(beta) <- gsub("X", "", colnames(beta))
colnames(low) <- colnames(beta)
colnames(medium) <- colnames(beta)
colnames(high) <- colnames(beta)
row.names(low) <- row.names(beta)
row.names(medium) <- row.names(beta)
row.names(high) <- row.names(beta)

low <- low[order(row.names(low)), ]
medium <- medium[order(row.names(medium)), ]
high <- high[order(row.names(high)), ]
snps <- snps[order(snps$IlmnID), ]



snps$SNP_low <- as.character(snps$SNP_low)
snps$SNP_medium <- as.character(snps$SNP_medium)
snps$SNP_high <- as.character(snps$SNP_high)

genotype <- matrix(nrow=nrow(low), ncol=ncol(low))
for (i in 1:ncol(low)){
  data <- data.frame(low[ ,i], medium[ ,i], high[ ,i])
  colnames(data) <- c("low", "medium", "high")
  row.names(data) <- row.names(low)
  for (j in 1:nrow(low)){
    if(data$low[j] != "NaN"){
      if(data$low[j] == max(data[j, ])){
        genotype[j,i] <- snps$SNP_low[j]
      } else if(data$medium[j] == max(data[j, ])){
        genotype[j,i] <- snps$SNP_medium[j]
      } else{
        genotype[j,i] <- snps$SNP_high[j]
      }
    }
  }
}

row.names(genotype) <- row.names(low)
colnames(genotype) <- colnames(low)

#this csv file contains the genotypes for each sample and SNP probe
write.table(genotype, "OUTPUTFILENAME.csv", sep=",", dec=".", col.names=TRUE, row.names=TRUE)

#The following code generates .map and .ped files which are needed for genotype matching with the GSA
mapfile <- matrix(nrow=nrow(low), ncol=4)
mapfile[ ,3] <- 0
mapfile[ ,4] <- snps$Start_hg38
mapfile[ ,1] <- gsub("chr", "", snps$CHR_hg38)
mapfile[ ,2] <- as.character(snps$IlmnID)

mapfile <- as.data.frame(mapfile)

write.table(mapfile, "MAPFILENAME.map", col.names=FALSE, row.names=FALSE, quote=FALSE)

pedfile <- matrix(nrow=ncol(low), ncol=(6+(nrow(low)*2)))

pedfile[ ,1] <- colnames(genotype)
pedfile[ ,2] <- colnames(genotype)
pedfile[ ,3] <- 0
pedfile[ ,4] <- 0
pedfile[ ,5] <- -9
pedfile[ ,6] <- -9
counter=7
for (i in 1:nrow(genotype)){
  pedfile[ ,counter] <- substr(genotype[i, ], 1, 1)
  counter=counter+1
  pedfile[ ,counter] <- substr(genotype[i, ], 3, 3)
  counter=counter+1
}

samplenames <- read.csv("SAMPLESHEETNAME.csv", sep=",", dec=".", header=TRUE)
row.names(samplenames) <- paste(samplenames$SentrixBarcode_A, "_", samplenames$SentrixPosition_A, sep="")
samplenames <- samplenames[row.names(samplenames) %in% pedfile[ ,1], ]
samplenames <- samplenames[order(row.names(samplenames)), ]
pedfile <- pedfile[order(pedfile[ ,1]), ]

pedfile[ ,1] <- as.character(samplenames$Sample_Name)
pedfile[ ,2] <- as.character(samplenames$Sample_Name)

write.table(pedfile, "PEDFILENAME.ped", col.names=FALSE, row.names=FALSE, quote=FALSE)