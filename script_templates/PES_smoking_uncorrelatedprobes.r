#This script can be used to extract uncorrelated CpG-probes for PES-calculations and to calculate PES
#It is similar to the filtering / clumping that is done prior to the PCA
#but instead of choosing a random CpG per 100kb bin, here, the CpG with the most extreme effect size in an EWAS statistic is chosen
#In addition, this process is repeated for several P-value thresholds in the EWAS statistics
#Here, the EWAS statistics from the Joehanes et al. publication for the current vs. never smoking phenotype were used

library("bigmelon")
library("ChAMP")

set.seed(147)

#This function can be used to extract uncorrelated probes from EWAS test-statistics according to a P-value threshold
#Note that depending on the EWAS statistics file, you may need to rename "Effect" to another column name
find_uncorrelated <- function(probes){
  data(probe.features.epic)
  probe.features <- probe.features[row.names(probe.features) %in% row.names(probes), ]
  probe.features <- probe.features[with(probe.features, order(probe.features$CHR, probe.features$MAPINFO)), ]
  probe.features <- probe.features[!is.na(probe.features$MAPINFO), ]
  probe.features$CHR <- as.character(probe.features$CHR)
  probe.features$CHR <- as.factor(probe.features$CHR)
  levels <- levels(probe.features$CHR)

  finalnr = 0
  for (i in 1:length(levels(probe.features$CHR))){
    temp <- probe.features[probe.features$CHR == levels[i], ]
    range <- max(temp$MAPINFO) - min(temp$MAPINFO)
    finalnr = finalnr + ceiling(range/100000)
  }

  finalcpgs <- rep(NA, finalnr)
  counter=1
  for (i in 1:length(levels(probe.features$CHR))){
    tempchr <- probe.features[probe.features$CHR == levels[i], ]
    range <- max(tempchr$MAPINFO) - min(tempchr$MAPINFO)
    for (j in 1:ceiling(range/100000)){
      temprangemin <- min(tempchr$MAPINFO)+((j-1)*100000)
      temprangemax <- temprangemin+100000
      tempcluster <- tempchr[which(tempchr$MAPINFO >= temprangemin & tempchr$MAPINFO < temprangemax), ]
        if(nrow(tempcluster) == 0){
        counter=counter+1
      } else{
        temp <- probes[row.names(probes) %in% row.names(tempcluster), ]
        temp$Effect <- sqrt((temp$Effect)^2)
        temp <- temp[order(-temp$Effect), ] 
        finalcpgs[counter] <- row.names(temp)[1]
        counter=counter+1
      }
    }
  }

  finalcpgs <- finalcpgs[!is.na(finalcpgs)]

  return(finalcpgs)
}

#This function can be used to calculate the PES
#Note that depending on the EWAS statistics file, you may need to rename "Effect" to another column name
multiply <- function(x){ 
    x*probes$Effect
}

#specify the name that should be included in all outputfiles here, e.g. "base2_blood_smokingpes"
outputfilename = "OUTPUTFILENAME"
print(outputfilename)

#specify the P-value thresholds here:
#Note: P<9E-08 is the EWAS epigenome-wide significance threshold for the EPIC array according to Mansell et al.
pvals = c(1E-03, 1E-04, 1E-05, 1E-06, 9E-08)

#read in DNAm data
setwd("PATH TO PROJECT DIRECTORY")
gfile <- openfn.gds("gdsfiles/PROJECTNAME_norm_filtered_celltyperegressedout.gds")
allbeta <- read.gdsn(index.gdsn(gfile, "betas"))
row.names(allbeta) <- read.gdsn(index.gdsn(gfile, "Probe_ID"))
colnames(allbeta) <- read.gdsn(index.gdsn(gfile, "Sample_ID"))
closefn.gds(gfile)

#read in the EWAS statistics you want to use for the PES calculations here. This is an example for the current vs. never smoking EWAS from Joehanes et al.
allprobes <- read.csv("/data/liga_ewas/sommerer_data/EWAS_stats/Joehanetal_stats/currentvsnever_joehanesetal.csv", sep=",", dec=".", header=TRUE, row.names=1)
allprobes$CpG <- row.names(allprobes)
allbeta <- allbeta[row.names(allbeta) %in% row.names(allprobes), ]
allprobes <- allprobes[row.names(allprobes) %in% row.names(allbeta), ]

setwd("PATH TO OUTPUT FILE DIRECTORY")








#AUTOMATED SCRIPT STARTS HERE - EVERYTHING BELOW SHOULD RUN AUTOMATICALLY BASED ON THE INPUT ABOVE
strt <- Sys.time()

###all probes (P<1) = 1###
print("P-value threshold 1: All probes")
print("1:: Number of probes prior to filtering")
print(nrow(allprobes))

finalcpgs <- find_uncorrelated(allprobes)

write.table(finalcpgs, file=paste(outputfilename, "1_cpgs.csv", sep="_"), sep=",", dec=".", col.names=FALSE, row.names=FALSE)

beta <- allbeta[row.names(allbeta) %in% finalcpgs, ]
probes <- allprobes[row.names(allprobes) %in% finalcpgs, ]

probes <- probes[order(row.names(probes)), ]
beta <- beta[order(row.names(beta)), ]

print("1:: Number of probes after filtering")
print(nrow(probes))

scores <- apply(beta, 2, multiply)
pes <- apply(scores, 2, sum)

write.table(pes, file=paste(outputfilename, "1_scores.csv", sep="_"), sep=",", dec=".", quote=FALSE, col.names=FALSE)

#Repeat this while filtering for P-values as specified by the P-value thresholds
for(i in 1:length(pvals)){
  probes <- allprobes[which(allprobes$P < pvals[i]), ]

  print(paste("P-value threshold", i+1, ": P < ", pvals[i], sep=""))
  print(paste(i+1, ":: Number of probes prior to filtering", sep=""))
  print(nrow(probes))

  if(nrow(probes) > 1){
    finalcpgs <- find_uncorrelated(probes)

    write.table(finalcpgs, file=paste(outputfilename, i+1, "cpgs.csv", sep="_"), sep=",", dec=".", col.names=FALSE, row.names=FALSE)

    beta <- allbeta[row.names(allbeta) %in% finalcpgs, ]
    probes <- probes[row.names(probes) %in% finalcpgs, ]

    probes <- probes[order(row.names(probes)), ]
    beta <- beta[order(row.names(beta)), ]

    print(paste(i+1, ":: Number of probes after filtering", sep=""))
    print(nrow(probes))

    scores <- apply(beta, 2, multiply)
    pes <- apply(scores, 2, sum)

    write.table(pes, file=paste(outputfilename, i+1, "scores.csv", sep="_"), sep=",", dec=".", quote=FALSE, col.names=FALSE)
  } else if (nrow(probes) == 0){
    print(paste("No CpGs with a P-value <", pvals[i]))
  } else {
    print(paste("Only one CpG with a P-value <", pvals[i]))
  }
}

print("Time from start to finish of script:")
print(Sys.time()-strt)

warnings()
