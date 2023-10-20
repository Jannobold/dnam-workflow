library("ChAMP")
library("bigmelon")

set.seed(147)

strt <- Sys.time()

setwd("PATH TO PROJECT")

gfile <- openfn.gds("gdsfiles/GDS FILENAME HERE.gds")
probes <- read.gdsn(index.gdsn(gfile, "betas"))
row.names(probes) <- read.gdsn(index.gdsn(gfile, "Probe_ID"))
colnames(probes) <- read.gdsn(index.gdsn(gfile, "Sample_ID"))
closefn.gds(gfile)

print("Number of probes prior to filtering")
print(nrow(probes))

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
    } else if (nrow(tempcluster) == 1){
      finalcpgs[counter] <- row.names(tempcluster)
      counter=counter+1
    } else{
      finalcpgs[counter] <- sample(row.names(tempcluster),1)
      counter=counter+1
    }
  }
}

finalcpgs <- finalcpgs[!is.na(finalcpgs)]

write.table(finalcpgs, file="pca/CPG FILENAME.csv", sep=",", dec=".", col.names=FALSE, row.names=FALSE)

print("Time from start to finish of script:")
print(Sys.time()-strt)

warnings()
