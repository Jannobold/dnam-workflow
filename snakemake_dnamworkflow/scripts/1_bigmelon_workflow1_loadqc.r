args = commandArgs(trailingOnly=TRUE)
logfile = args[[1]]
idatdir = args[[2]]
workinggds = args[[3]]
rawgds <- args[[4]]
rawfiltered <- args[[5]]
outlyxout <- args[[6]]
bsconout <- args[[7]]
outliernames <- args[[8]]

###this script can be used to load and perform qc on EPIC array data###

set.seed(147)

sink(logfile, append=TRUE, split=TRUE)

print("<<<<<<<<<<<<<<<<<<<<")

print(date())
print("Run Loading & QC")

#The following parameters can be changed depending on the dataset
print("Parameters:")
bsconv_threshold <- 80
outlier_threshold <- 0.15
chunks <- 100
print(paste("Bisulfite conversion efficiency: ", bsconv_threshold, sep=""))
print(paste("Outlier mvP: ", outlier_threshold, sep=""))
print(paste("Chunk size: ", chunks, sep=""))

library("bigmelon")

gfile <- iadd2(idatdir, gds = workinggds, chunksize = chunks)
probeids <- read.gdsn(index.gdsn(gfile, 'fData/Probe_ID'))
print("Number of CpG-probes")
print(length(probeids))
rawmet <- methylated(gfile)[,] 
rawume <- unmethylated(gfile)[,]
newbeta <- rawmet / (rawmet+rawume+100)
add.gdsn(gfile, "betas", newbeta, replace=TRUE)

f <- createfn.gds(rawgds)
backup.gdsn(gds = f, node = index.gdsn(gfile, 'betas'))
probeids <- read.gdsn(index.gdsn(gfile, 'fData/Probe_ID'))
add.gdsn(f, 'Probe_ID', probeids)
sampleids <- read.gdsn(index.gdsn(gfile, 'pData/barcode'))
add.gdsn(f, "Sample_ID", sampleids)
closefn.gds(f) 

sink(stdout(), type = "message")
pfilter(gfile)

f <- createfn.gds(rawfiltered)
backup.gdsn(gds = f, node = index.gdsn(gfile, 'betas'))
probeids <- read.gdsn(index.gdsn(gfile, 'fData/Probe_ID'))
add.gdsn(f, 'Probe_ID', probeids)
sampleids <- read.gdsn(index.gdsn(gfile, 'pData/barcode'))
add.gdsn(f, "Sample_ID", sampleids)
closefn.gds(f) 

beta <- betas(gfile)[,]
outlier <- outlyx(gfile, plot=FALSE, mvP = outlier_threshold)
conversion <- bscon(gfile)
write.table(outlier, outlyxout, col.names=TRUE, row.names=TRUE)
write.table(conversion, bsconout, col.names=TRUE, row.names=TRUE)
closefn.gds(gfile)

conversion <- read.table(bsconout, header=TRUE, row.names=1)
outlier <- read.table(outlyxout, header=TRUE, row.names=1)
row.names(conversion) <- row.names(outlier)
write.table(conversion, bsconout, col.names=TRUE, row.names=TRUE)

print("Samples that have a low bisulfite conversion efficiency:")
print(length(row.names(conversion)[which(conversion$x < bsconv_threshold)]))
print(row.names(conversion)[which(conversion$x < bsconv_threshold)])
print(conversion[which(conversion$x < bsconv_threshold), ])

print("Samples that are outliers:")
print(nrow(outlier[which(outlier$outliers == TRUE), ]))
print(outlier[which(outlier$outliers == TRUE), ])

conversion <- read.table(bsconout, header=TRUE, row.names=1)
outlier <- read.table(outlyxout, header=TRUE, row.names=1)
conversion <- row.names(conversion)[which(conversion$x < bsconv_threshold)]
outlier <- outlier[which(outlier$outliers == TRUE), ]

alloutliers <- c(conversion, row.names(outlier))
alloutliers <- unique(alloutliers)
write.table(alloutliers, outliernames, row.names=FALSE, col.names=FALSE)

print("If any samples are outliers or have a low bisulfite conversion efficiency, re-run this script with problematic samples removed.")
print("Otherwise, proceed with normalization and qc.")

print("<<<<<<<<<<<<<<<<<<<<")

sink()