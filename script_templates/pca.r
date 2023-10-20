library("ggplot2")
library("bigmelon")

set.seed(147)

setwd("PATH TO PROJECT HERE")

gfile <- openfn.gds("gdsfiles/GDS FILE NAME.gds")
beta <- read.gdsn(index.gdsn(gfile, "betas"))
row.names(beta) <- read.gdsn(index.gdsn(gfile, "Probe_ID"))
colnames(beta) <- read.gdsn(index.gdsn(gfile, "Sample_ID"))
closefn.gds(gfile)

probes <- read.csv("pca/CPG FILENAME.csv", sep=",", dec=".", header=FALSE)

beta <- beta[row.names(beta) %in% probes$V1, ]

setwd("PATH TO PHENOTYPE")

pheno <- read.csv("PHENOTYPEFILE.csv", sep=",", dec=".", header=TRUE, row.names=1)

beta <- beta[ ,colnames(beta) %in% pheno$array_id]

pca <- prcomp(t(beta))
pcs <- data.frame(pca$x)

var_explained <- summary(pca)$importance[2,1:20]

plotdata <- data.frame(1:20, var_explained)
colnames(plotdata) <- c("pc", "var")

png("SCREEPLOT FILENAME.png", height=500, width=500, res=100)
ggplot(plotdata, aes(x=pc, y=var)) + 
  geom_point() + 
  geom_line() +
  theme_bw() + 
  xlab("PC") +
  ylab("Explained Variance") +
  scale_x_continuous(breaks=1:19) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

write.table(pcs[ ,1:20], "PCS FILENAME.csv", sep=",", dec=".", col.names=TRUE, row.names=TRUE)
