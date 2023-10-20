#The code for the Manhattan plot is based on this: https://r-graph-gallery.com/101_Manhattan_plot.html

set.seed(147)

setwd("PATH TO DIRECTORY WITH TEST STATISTICS") #directory should only contain csvs with test statistics

library("ggplot2")
library("ggrepel")
library("dplyr")
library("ChAMP")
library("qqman")

filelist <- list.files()

lambdas <- rep(NA, length(filelist))

for(i in 1:length(filelist)){
  data(probe.features.epic)
  
  data <- read.csv(filelist[i], sep=",", dec=".", header=TRUE, row.names=1)
  data <- data[which(!(is.na(data$P))), ]
  
  probe.features <- probe.features[row.names(probe.features) %in% row.names(data), ] #data is a matrix of test statistic results
  probe.features <- probe.features[order(row.names(probe.features)), ]

  data <- data[order(data$P), ]
  snpsOfInterest <- row.names(data)[which(data$P < 1E-05)] #is a vector of strings, representing the ids of all probes that should be highlighted in red
  snpsOfGenes <- snpsOfInterest #is a vector of strings, representing the ids of all probes that should be annotated with the gene name

  data <- data[order(row.names(data)), ]

  probe.features$CHR <- factor(probe.features$CHR, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", ""))

  probe.features$MAPINFO <- as.numeric(probe.features$MAPINFO)
  
  spearmanResults <- data.frame(row.names(probe.features), probe.features$CHR, probe.features$MAPINFO, data$P, probe.features$gene)
  colnames(spearmanResults) <- list("SNP", "CHR", "BP", "P", "Gene")

  don <- spearmanResults %>%
    group_by(CHR) %>%
    summarise(chr_len=max(BP)) %>%
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    dplyr::select(-chr_len) %>%
    left_join(spearmanResults, ., by=c("CHR"="CHR")) %>%
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%
    mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
    mutate( is_annotate=ifelse(SNP %in% snpsOfGenes, "yes", "no"))

  axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

  png(filename=gsub(".csv", "_manhattan_forpaper.png", filelist[i]), width=1000, height=500, units="px", res=100)
  print(ggrastr::rasterize(ggplot(don, aes(x=BPcum, y=-log10(P))) +
          geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1) +
          scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
          scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
          scale_y_continuous(limits = c(0, 10) ) +
          geom_point(data=subset(don, is_highlight=="yes"), color="purple", size=1) +
          geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=Gene),
                            size=2.5, min.segment.length = 0,
                            ylim     = c(7.5, NA),
                            direction    = "x",
                            angle        = 90,
                            hjust        = 0) +
          labs(x = "Chromosome") +
          geom_hline(yintercept=-log10(9E-08), color = "red") +
          geom_hline(yintercept=-log10(1E-05), color = "purple") +
          theme_bw() +
          theme(
            legend.position="none",
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()
          ), layers = "point", dpi = 300))
  dev.off()

  if(max(-log10(data$P) >= 10 )){
  png(filename=gsub(".csv", "_manhattan_full.png", filelist[i]), width=1000, height=500, units="px", res=100)
              print(ggplot(don, aes(x=BPcum, y=-log10(P))) +
              geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
              scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
              scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
              scale_y_continuous(expand = c(0, 1) ) +
              geom_point(data=subset(don, is_highlight=="yes"), color="red", size=2) +
              geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=Gene), size=5) +
              labs(x = "Chromosome") +
              geom_hline(yintercept=-log10(1.64E-07), color = "red") +
              theme_bw() +
              theme(
               legend.position="none",
               panel.border = element_blank(),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank()
             ))
     dev.off()
  }

   png(filename=gsub(".csv", "_qq.png", filelist[i]), width=500, height=500)
   print(qq(data$P, main = "")) #pvals is a numerical vector of p-values
   dev.off()
   
   chisq <- qchisq(data$P,1, lower.tail=FALSE)
   lambdas[i] = median(chisq)/qchisq(0.5,1)
}

 names(lambdas) <- filelist
 write.table(lambdas, "lambdas.csv", col.names=FALSE, row.names=TRUE, sep=",", dec=".")
