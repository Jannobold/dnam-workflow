args = commandArgs(trailingOnly=TRUE)
logfile = args[[1]]
workinggds = args[[2]]
normunfiltered = args[[3]]
rawfiltered = args[[4]]
qualout = args[[5]]

##this script can used to normalize and perform qc on EPIC array data###

set.seed(147)

sink(logfile, append=TRUE, split=TRUE)

print("<<<<<<<<<<<<<<<<<<<<")

print(date())
print("Run Normalization & QC")

library("bigmelon")

gfile <- openfn.gds(workinggds, readonly=FALSE)

  dasen(gfile)

  f <- createfn.gds(normunfiltered)
    backup.gdsn(gds = f, node = index.gdsn(gfile, 'betas'))
  closefn.gds(f)
  gfile2 <- openfn.gds(rawfiltered)
    raw <- index.gdsn(gfile2, "backup/betas")
    output <- qual(norm=betas(gfile), raw=raw)
    write.table(output, qualout, sep=",", dec=".", col.names=TRUE, row.names=TRUE)
  closefn.gds(gfile2)
closefn.gds(gfile)

print("Samples with a large change of beta values after normalization:")
output <- as.data.frame(output)
print(output[which(output$rmsd > 0.1), ])
output <- output[which(output$rmsd > 0.1), ]
output <- row.names(output)

write.table(output, qualout, sep=",", dec=".", col.names=FALSE, row.names=FALSE)

print("If any samples show a large change of beta values after normalization, re-run everything (loading, qc, normalization, qc) with problematic samples removed.")
print("Otherwise, proceed with filtering.")
print("<<<<<<<<<<<<<<<<<<<<")

sink()