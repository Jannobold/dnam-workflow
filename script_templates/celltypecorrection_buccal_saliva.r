library("bigmelon")
library("limma")

set.seed(147)

setwd("PATH TO PROJECT HERE")

gfile <- openfn.gds("gdsfiles/GDS FILE NAME.gds")
dnam <- read.gdsn(index.gdsn(gfile, "betas"))
row.names(dnam) <- read.gdsn(index.gdsn(gfile, "Probe_ID"))
colnames(dnam) <- read.gdsn(index.gdsn(gfile, "Sample_ID"))

closefn.gds(gfile)

teccovar <- read.csv("results/CELL TYPE ESTIMATION.csv", sep=",", dec=".", header=TRUE, row.names=1)

test <- removeBatchEffect(dnam, covariates=teccovar$Epi)

#create a new gds file with the cell-type corrected data
#ensure that you are not overwriting an already existing gds file
f <- createfn.gds("gdsfiles/GDS FILE NAME.gds")
add.gdsn(f, 'betas', test)
probeids <- row.names(test)
add.gdsn(f, 'Probe_ID', probeids)
sampleids <- colnames(test)
add.gdsn(f, "Sample_ID", sampleids)
closefn.gds(f)
