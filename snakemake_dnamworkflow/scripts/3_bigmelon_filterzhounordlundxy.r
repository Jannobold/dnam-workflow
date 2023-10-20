args = commandArgs(trailingOnly=TRUE)
logfile = args[[1]]
workinggds = args[[2]]
zhouloc = args[[3]]
nordloc = args[[4]]
finalcsv = args[[5]]
normfiltered = args[[6]]

###this script can used to filter CpG probes according to Zhou et al., Nordlund et al., and X & Y chromosomes on EPIC array data###

set.seed(147)

sink(logfile, append=TRUE, split=TRUE)

print("<<<<<<<<<<<<<<<<<<<<")

print(date())
print("Run Filtering according to Zhou et al., Nordlund et al., and X & Y chromosomes")

library("bigmelon")
library("ChAMP")

zhou <- read.table(zhouloc, sep="\t", header=TRUE)
nord <- read.table(nordloc, sep="\t", header=TRUE)

gfile <- openfn.gds(workinggds)
beta <- betas(gfile)[,]
zhou <- zhou[zhou$MASK.general == TRUE, ]
nord <- nord[nord$bwa.multi.hit == 1, ]

print("Number of probes excluded due to SNP alignment (Zhou et al.)")
print(nrow(beta[row.names(beta) %in% zhou$probeID, ]))
beta <- beta[!(row.names(beta) %in% zhou$probeID), ]

print("Number of probes excluded due to aligning to multiple locations (Nordlund et al.)")
print(nrow(beta[row.names(beta) %in% nord$TargetID, ]))
beta <- beta[!(row.names(beta) %in% nord$TargetID), ]

data(probe.features.epic)
probe.features <- probe.features[row.names(probe.features) %in% row.names(beta), ]
probe.features <- probe.features[which(probe.features$CHR != "X"), ]
probe.features <- probe.features[which(probe.features$CHR != "Y"), ]

print("Number of probes excluded from X and Y Chromosome")
print(nrow(beta[!(row.names(beta) %in% row.names(probe.features)), ]))
beta <- beta[row.names(beta) %in% row.names(probe.features), ]

probe.features <- probe.features[which(probe.features$CHR != ""), ]
print("Number of non-cg probes excluded")
print(nrow(beta[!(row.names(beta) %in% row.names(probe.features)), ]))
beta <- beta[row.names(beta) %in% row.names(probe.features), ]

print("final number of CpG probes")
print(nrow(beta))

f <- createfn.gds(normfiltered)
add.gdsn(f, 'betas', beta)
probeids <- row.names(beta)
add.gdsn(f, 'Probe_ID', probeids)
sampleids <- colnames(beta)
add.gdsn(f, "Sample_ID", sampleids)
closefn.gds(f)


write.table(beta, finalcsv, sep=",", dec=".", col.names=TRUE, row.names=TRUE)

closefn.gds(gfile)

print("<<<<<<<<<<<<<<<<<<<<")

sink()