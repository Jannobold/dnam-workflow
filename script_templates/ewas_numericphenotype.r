library("bigmelon")
library("dplyr")

strt <- Sys.time()

set.seed(147)

setwd("PATH TO PROJECT")

gfile <- openfn.gds("gdsfiles/GDS FILENAME.gds")
dnam <- read.gdsn(index.gdsn(gfile, "betas"))
row.names(dnam) <- read.gdsn(index.gdsn(gfile, "Probe_ID"))
colnames(dnam) <- read.gdsn(index.gdsn(gfile, "Sample_ID"))

closefn.gds(gfile)

pheno <- read.csv("PHENOTYPE.csv", sep=",", dec=".", header=TRUE)
pcs <- read.csv("PCS.csv", sep=",", dec=".", header=TRUE, row.names=1)

pcs$array_id <- row.names(pcs)

pheno <- merge(pheno, pcs, by.x="array_id_blood", by.y="array_id", suffixes=c("",""),
  all.x=FALSE, all.y=FALSE
)
row.names(pheno) <- pheno$array_id_blood

dnam <- dnam[ ,colnames(dnam) %in% row.names(pheno)]
dnam <- dnam[ ,order(colnames(dnam))]
pheno <- pheno[order(row.names(pheno)), ]

#Initializing the matrices where the EWAS results will be saved
temp_male <- matrix(NA, nrow(dnam), 4)
colnames(temp_male) <- c("CpG", "P", "beta", "SE")
temp_male[ ,1] <- row.names(dnam)

temp_female <- matrix(NA, nrow(dnam), 4)
colnames(temp_female) <- c("CpG", "P", "beta", "SE")
temp_female[ ,1] <- row.names(dnam)

temp_both <- matrix(NA, nrow(dnam), 4)
colnames(temp_both) <- c("CpG", "P", "beta", "SE")
temp_both[ ,1] <- row.names(dnam)

for(i in 1:nrow(dnam)){
  probe <- as.vector(t(dnam[i, ]))
  testdata <- cbind(pheno, probe, pcs)
  males <- testdata[which(testdata$geschlecht == "1"), ]
  females <- testdata[which(testdata$geschlecht == "2"), ]
  testdata$geschlecht <- as.factor(testdata$geschlecht)
  #doing as.character first because otherwise, if the variable is coded as a factor,
  #it will take the numeric value of the factor level instead of the numeric value itself <.<
  testdata$NUMERIC_VARIABLE <- as.character(as.numeric(testdata$NUMERIC_VARIABLE))
  #[DNAM PCS HERE] - e.g., PC1 + PC2 + PC3 etc., the number of PCs that should be included can be determined by a scree plot
  #NUMERIC VARIABLE is a numeric, continuous variable
  linreg_m <- lm(formula = probe ~ NUMERIC_VARIABLE + Age_T1 + [DNAM PCS HERE], data=males)
  linreg_f <- lm(formula = probe ~ NUMERIC_VARIABLE + Age_T1 + [DNAM PCS HERE], data=females)
  linreg <- lm(formula = probe ~ NUMERIC_VARIABLE + geschlecht + Age_T1 + [DNAM PCS HERE], data=testdata)
  
  
  coeff <- summary(linreg_m)$coefficients[,4]
  temp_male[i,2] <- unname(coeff[2])
  coeff <- summary(linreg_m)$coefficients[,1]
  temp_male[i,3] <- unname(coeff[2])
  coeff <- summary(linreg_m)$coefficients[,2]
  temp_male[i,4] <- unname(coeff[2])


  coeff <- summary(linreg_f)$coefficients[,4]
  temp_female[i,2] <- unname(coeff[2])
  coeff <- summary(linreg_f)$coefficients[,1]
  temp_female[i,3] <- unname(coeff[2])
  coeff <- summary(linreg_f)$coefficients[,2]
  temp_female[i,4] <- unname(coeff[2])

  coeff <- summary(linreg)$coefficients[,4]
  temp_both[i,2] <- unname(coeff[2])
  coeff <- summary(linreg)$coefficients[,1]
  temp_both[i,3] <- unname(coeff[2])
  coeff <- summary(linreg)$coefficients[,2]
  temp_both[i,4] <- unname(coeff[2])
}

setwd("OUTPUT PATH")

write.table(temp_male, "OUTPUT FILE MALES.csv", sep=",", dec=".", col.names=TRUE, row.names=FALSE)
write.table(temp_female, "OUTPUT FILE FEMALES.csv", sep=",", dec=".", col.names=TRUE, row.names=FALSE)
write.table(temp_both, "OUTPUT FILE ALL.csv", sep=",", dec=".", col.names=TRUE, row.names=FALSE)

print("Time from start to finish of script:")
print(Sys.time()-strt)

