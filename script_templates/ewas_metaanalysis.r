#This script can be used to conduct a fixed-effect inverse-variance meta-analysis of EWAS test statistics

library("meta")
library("ChAMP")

set.seed(147)

setwd("PATH TO OUTPUT DIRECTORY")

#example for 3 ewas test statistics to be combined in the meta-analysis; more files can be added accordingly
#the first column name should contain the CpG ID, and in addition the files should contain the following columns:
#beta (estimate of the effect size), SE (standard error)
data1 <- read.csv("EWAS_STATISTICS1.csv", sep=",", dec=".", header=TRUE, row.names=1)
data2 <- read.csv("EWAS_STATISTICS2.csv", sep=",", dec=".", header=TRUE, row.names=1)
data3 <- read.csv("EWAS_STATISTICS3.csv", sep=",", dec=".", header=TRUE, row.names=1)

#ensure to only include CpGs that are present in all datasets
data1 <- data1[row.names(data1) %in% row.names(data2), ]
data1 <- data1[row.names(data1) %in% row.names(data3), ]
data2 <- data2[row.names(data2) %in% row.names(data1), ]
data3 <- data3[row.names(data3) %in% row.names(data1), ]

#ensure that all CpGs are in the same order in all datasets
data1 <- data1[order(row.names(data1)), ]
data2 <- data2[order(row.names(data2)), ]
data3 <- data3[order(row.names(data3)), ]

#remove non-cg CpG-probes based on EPIC array annotation
data(probe.features.epic)

probe.features <- probe.features[row.names(probe.features) %in% row.names(data1), ]
probe.features <- probe.features[which(probe.features$CHR != ""), ]
data1 <- data1[row.names(data1) %in% row.names(probe.features), ]
data2 <- data2[row.names(data2) %in% row.names(probe.features), ]
data3 <- data3[row.names(data3) %in% row.names(probe.features), ]

#combine effect size estimates (beta) and standard errors (SE) of all datasets
TE <- cbind(data1$beta, data2$beta, data3$beta)
seTE <- cbind(data1$SE, data2$SE, data3$SE)

temp_fixed <- matrix(NA, nrow(TE), 6)
colnames(temp_fixed) <- c("CpG", "beta", "SE", "lowerbCI", "upperbCI", "P")
temp_fixed[ ,1] <- row.names(data1)

#Run the fixed-effect inverse-variance meta-analysis
#Note: in some cases, a random-effect meta-analysis may be needed
#In that case, replace the suffix .fixed by .random
for(i in 1:nrow(TE)){
  meta <- metagen(TE[i, ],  seTE[i, ])
  temp_fixed[i,2] <- meta$TE.fixed
  temp_fixed[i,3] <- meta$seTE.fixed
  temp_fixed[i,4] <- meta$lower.fixed
  temp_fixed[i,5] <- meta$upper.fixed
  temp_fixed[i,6] <- meta$pval.fixed
}

#The output file contains the CpG-ID, effect size estimate (beta), standard error (SE)
#lower boundary of the confidence interval (lowerbCI), upper boundary of the confidence
#interval (upperbCI), and P-value (P) of the meta-analysis
write.table(temp_fixed, "OUTPUTNAME.csv", sep=",", dec=".", col.names=TRUE, row.names=FALSE)
