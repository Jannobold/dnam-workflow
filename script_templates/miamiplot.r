library("hudson")
library("ChAMP")

set.seed(147)

#YS: adapted code from the hudson package function which allows for more in-depth customization of the plot
#YS: In particular, note that you can change the colors for each aspect of the plot in the following variables:
#YS: highlighter, chrcolor1_1, chrcolor1_2, chrcolor2_1, chrcolor2_2
gmirror <- function (top, bottom, tline, bline, log10 = TRUE, yaxis, opacity = 1, 
                     annotate_snp, annotate_p, toptitle = NULL, bottomtitle = NULL, 
                     highlight_snp, highlight_p, highlighter = "red", chrcolor1_1 = "#6f1e71", 
                     chrcolor1_2 = "#b02eb4", chrcolor2_1 = "#212b81", 
                     chrcolor2_2 = "#3d4abd", freey = FALSE, background = "variegated", 
                     chrblocks = FALSE, file, hgt = 7, hgtratio = 0.5, 
                     wi = 12, res = 300) 
{
  topn <- names(top)
  bottomn <- names(bottom)
  top$Location <- "Top"
  bottom$Location <- "Bottom"
  d <- rbind(top, bottom)
  d$POS <- as.numeric(as.character(d$POS))
  d$CHR <- factor(d$CHR, levels = c("1", "2", "3", 
                                    "4", "5", "6", "7", "8", 
                                    "9", "10", "11", "12", "13", 
                                    "14", "15", "16", "17", "18", 
                                    "19", "20", "21", "22", "X", 
                                    "Y"))
  d_order <- d[order(d$CHR, d$POS), ]
  d_order$pos_index <- seq.int(nrow(d_order))
  d_order_sub <- d_order[, c("SNP", "CHR", "POS", 
                             "pvalue", "pos_index")]
  maxRows <- by(d_order_sub, d_order_sub$CHR, function(x) x[which.max(x$pos_index), 
                                                            ])
  minRows <- by(d_order_sub, d_order_sub$CHR, function(x) x[which.min(x$pos_index), 
                                                            ])
  milimits <- do.call(rbind, minRows)
  malimits <- do.call(rbind, maxRows)
  lims <- merge(milimits, malimits, by = "CHR")
  names(lims) <- c("Color", "snpx", "px", 
                   "posx", "posmin", "snpy", "py", 
                   "posy", "posmax")
  lims$av <- (lims$posmin + lims$posmax)/2
  lims <- lims[order(lims$Color), ]
  lims$shademap <- rep(c("shade_ffffff", "shade_ebebeb"), 
                       length.out = nrow(lims), each = 1)
  nchrcolors <- nlevels(factor(lims$Color))
  colnames(d_order)[2] <- "Color"
  newcols_1 <- c(rep(x = c(chrcolor1_1, chrcolor1_2), length.out = nchrcolors, 
                     each = 1), "#FFFFFF", "#EBEBEB")
  names(newcols_1) <- c(levels(factor(lims$Color)), "shade_ffffff", 
                        "shade_ebebeb")
  newcols_2 <- c(rep(x = c(chrcolor2_1, chrcolor2_2), length.out = nchrcolors, 
                     each = 1), "#FFFFFF", "#EBEBEB")
  names(newcols_2) <- c(levels(factor(lims$Color)), "shade_ffffff", 
                        "shade_ebebeb")
  if (log10 == TRUE) {
    d_order$pval <- -log10(d_order$pvalue)
    yaxislab1 <- expression(paste("-log"[10], "(p-value)", 
                                  sep = ""))
    yaxislab2 <- expression(paste("-log"[10], "(p-value)", 
                                  sep = ""))
    if (!missing(tline)) {
      tredline <- -log10(tline)
    }
    if (!missing(bline)) {
      bredline <- -log10(bline)
    }
  }
  else {
    d_order$pval <- d_order$pvalue
    yaxislab1 <- yaxis[1]
    yaxislab2 <- yaxis[2]
    if (!missing(tline)) {
      tredline <- tline
    }
    if (!missing(bline)) {
      bredline <- bline
    }
  }
  yaxismax1 <- ifelse(freey == FALSE, max(d_order$pval[which(d_order$pval < 
                                                               Inf)]), max(d_order$pval[which(d_order$pval < Inf) & 
                                                                                          d_order$Location == "Top"]))
  yaxismax2 <- ifelse(freey == FALSE, max(d_order$pval[which(d_order$pval < 
                                                               Inf)]), max(d_order$pval[which(d_order$pval < Inf) & 
                                                                                          d_order$Location == "Bottom"]))
  yaxismin1 <- ifelse(freey == FALSE, 0, min(d_order$pval[d_order$Location == 
                                                            "Top"]))
  yaxismin2 <- ifelse(freey == FALSE, 0, min(d_order$pval[d_order$Location == 
                                                            "Bottom"]))
  backpanel1 <- ifelse(background == "white", "NULL", 
                       "geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = yaxismin1, ymax = Inf, fill=factor(shademap)), alpha = 0.5)")
  backpanel2 <- ifelse(background == "white", "NULL", 
                       "geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = yaxismin2, ymax = Inf, fill=factor(shademap)), alpha = 0.5)")
  p1 <- ggplot() + eval(parse(text = backpanel1))
  if ("Shape" %in% topn) {
    p1 <- p1 + geom_point(data = d_order[d_order$Location == 
                                           "Top", ], aes(x = pos_index, y = pval, color = factor(Color), 
                                                         shape = factor(Shape)), alpha = opacity)
  }
  else {
    p1 <- p1 + geom_point(data = d_order[d_order$Location == 
                                           "Top", ], aes(x = pos_index, y = pval, color = factor(Color)), 
                          alpha = opacity)
  }
  p1 <- p1 + scale_x_continuous(breaks = lims$av, labels = lims$Color, 
                                expand = c(0, 0))
  if (chrblocks == TRUE) {
    p1 <- p1 + geom_rect(data = lims, aes(xmin = posmin - 
                                            0.5, xmax = posmax + 0.5, ymin = -Inf, ymax = min(d_order$pval), 
                                          fill = as.factor(Color)), alpha = 1)
  }
  p1 <- p1 + scale_colour_manual(name = "Color", values = newcols_1) + 
    scale_fill_manual(name = "Color", values = newcols_1)
  p1 <- p1 + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), 
                   axis.title.x = element_blank(), legend.position = "top", 
                   legend.title = element_blank())
  p2 <- ggplot() + eval(parse(text = backpanel2))
  if ("Shape" %in% bottomn) {
    p2 <- p2 + geom_point(data = d_order[d_order$Location == 
                                           "Bottom", ], aes(x = pos_index, y = pval, color = factor(Color), 
                                                            shape = factor(Shape)), alpha = opacity)
  }
  else {
    p2 <- p2 + geom_point(data = d_order[d_order$Location == 
                                           "Bottom", ], aes(x = pos_index, y = pval, color = factor(Color)), 
                          alpha = opacity)
  }
  p2 <- p2 + scale_x_continuous(breaks = lims$av, labels = lims$Color, 
                                expand = c(0, 0))
  if (chrblocks == TRUE) {
    p2 <- p2 + geom_rect(data = lims, aes(xmin = posmin - 
                                            0.5, xmax = posmax + 0.5, ymin = -Inf, ymax = min(d_order$pval), 
                                          fill = as.factor(Color)), alpha = 1)
  }
  p2 <- p2 + scale_colour_manual(name = "Color", values = newcols_2) + 
    scale_fill_manual(name = "Color", values = newcols_2)
  p2 <- p2 + theme(axis.text.x = element_text(angle = 90), 
                   panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), 
                   axis.title.x = element_blank(), legend.position = "bottom", 
                   legend.title = element_blank())
  if (!missing(highlight_snp)) {
    if ("Shape" %in% topn) {
      p1 <- p1 + geom_point(data = d_order[d_order$SNP %in% 
                                             highlight_snp & d_order$Location == "Top", 
                                           ], aes(x = pos_index, y = pval, shape = Shape), 
                            colour = highlighter)
      p1 <- p1 + guides(shape = guide_legend(override.aes = list(colour = "black")))
    }
    else {
      p1 <- p1 + geom_point(data = d_order[d_order$SNP %in% 
                                             highlight_snp & d_order$Location == "Top", 
                                           ], aes(x = pos_index, y = pval), colour = highlighter)
    }
    if ("Shape" %in% bottomn) {
      p2 <- p2 + geom_point(data = d_order[d_order$SNP %in% 
                                             highlight_snp & d_order$Location == "Bottom", 
                                           ], aes(x = pos_index, y = pval, shape = Shape), 
                            colour = highlighter)
      p2 <- p2 + guides(shape = guide_legend(override.aes = list(colour = "black")))
    }
    else {
      p2 <- p2 + geom_point(data = d_order[d_order$SNP %in% 
                                             highlight_snp & d_order$Location == "Bottom", 
                                           ], aes(x = pos_index, y = pval), colour = highlighter)
    }
  }
  if (!missing(highlight_p)) {
    if ("Shape" %in% topn) {
      p1 <- p1 + geom_point(data = d_order[d_order$pvalue < 
                                             highlight_p[1] & d_order$Location == "Top", 
                                           ], aes(x = pos_index, y = pval, shape = Shape), 
                            colour = highlighter)
      p1 <- p1 + guides(shape = guide_legend(override.aes = list(colour = "black")))
    }
    else {
      p1 <- p1 + geom_point(data = d_order[d_order$pvalue < 
                                             highlight_p[1] & d_order$Location == "Top", 
                                           ], aes(x = pos_index, y = pval), colour = highlighter)
    }
    if ("Shape" %in% bottomn) {
      p2 <- p2 + geom_point(data = d_order[d_order$pvalue < 
                                             highlight_p[2] & d_order$Location == "Bottom", 
                                           ], aes(x = pos_index, y = pval, shape = Shape), 
                            colour = highlighter)
      p2 <- p2 + guides(shape = guide_legend(override.aes = list(colour = "black")))
    }
    else {
      p2 <- p2 + geom_point(data = d_order[d_order$pvalue < 
                                             highlight_p[2] & d_order$Location == "Bottom", 
                                           ], aes(x = pos_index, y = pval), colour = highlighter)
    }
  }
  if (!missing(tline)) {
    p1 <- p1 + geom_hline(yintercept = tredline, colour = "red")
  }
  if (!missing(bline)) {
    p2 <- p2 + geom_hline(yintercept = bredline, colour = "red")
  }
  if (!missing(annotate_p)) {
    if (!requireNamespace(c("ggrepel"), quietly = TRUE) == 
        TRUE) {
      print("Consider installing 'ggrepel' for improved text annotation")
      p1 <- p1 + geom_text(data = d_order[d_order$pvalue < 
                                            annotate_p[1] & d_order$Location == "Top", 
                                          ], aes(pos_index, pval, label = SNP))
      p2 <- p2 + geom_text(data = d_order[d_order$pvalue < 
                                            annotate_p[2] & d_order$Location == "Bottom", 
                                          ], aes(pos_index, pval, label = SNP))
    }
    else {
      p1 <- p1 + ggrepel::geom_text_repel(data = d_order[d_order$pvalue < 
                                                           annotate_p[1] & d_order$Location == "Top", 
                                                         ], aes(pos_index, pval, label = SNP))
      p2 <- p2 + ggrepel::geom_text_repel(data = d_order[d_order$pvalue < 
                                                           annotate_p[2] & d_order$Location == "Bottom", 
                                                         ], aes(pos_index, pval, label = SNP))
    }
  }
  if (!missing(annotate_snp)) {
    if (!requireNamespace(c("ggrepel"), quietly = TRUE) == 
        TRUE) {
      print("Consider installing 'ggrepel' for improved text annotation")
      p1 <- p1 + geom_text(data = d_order[d_order$SNP %in% 
                                            annotate_snp & d_order$Location == "Top", 
                                          ], aes(pos_index, pval, label = SNP))
      p2 <- p2 + geom_text(data = d_order[d_order$SNP %in% 
                                            annotate_snp & d_order$Location == "Bottom", 
                                          ], aes(pos_index, pval, label = SNP))
    }
    else {
      p1 <- p1 + ggrepel::geom_text_repel(data = d_order[d_order$SNP %in% 
                                                           annotate_snp & d_order$Location == "Top", 
                                                         ], aes(pos_index, pval, label = SNP))
      p2 <- p2 + ggrepel::geom_text_repel(data = d_order[d_order$SNP %in% 
                                                           annotate_snp & d_order$Location == "Bottom", 
                                                         ], aes(pos_index, pval, label = SNP))
    }
  }
  p1 <- p1 + ylab(yaxislab1)
  p2 <- p2 + ylab(yaxislab2)
  if (chrblocks == TRUE) {
    if (freey == TRUE) {
      print("Sorry, drawing chrblocks with freey=TRUE is currently unsupported and will be ignored.")
    }
    else {
      p1 <- p1 + theme(axis.text.x = element_text(vjust = 1), 
                       axis.ticks.x = element_blank()) + ylim(c(yaxismin1, 
                                                                yaxismax1))
      p2 <- p2 + scale_y_reverse(limits = c(yaxismax2, 
                                            yaxismin2)) + theme(axis.text.x = element_blank(), 
                                                                axis.ticks.x = element_blank())
    }
  }
  else {
    p1 <- p1 + theme(axis.text.x = element_text(vjust = 1), 
                     axis.ticks.x = element_blank()) + scale_y_continuous(limits = c(yaxismin1, 
                                                                                     yaxismax1), expand = expansion(mult = c(0, 0.1)))
    p2 <- p2 + scale_y_reverse(limits = c(yaxismax2, yaxismin2), 
                               expand = expansion(mult = c(0.1, 0))) + theme(axis.text.x = element_blank(), 
                                                                             axis.ticks.x = element_blank())
  }
  if (background == "white") {
    p1 <- p1 + theme(panel.background = element_rect(fill = "white"))
    p2 <- p2 + theme(panel.background = element_rect(fill = "white"))
  }
  p1 <- p1 + guides(fill = "none", color = "none")
  p2 <- p2 + guides(fill = "none", color = "none")
  print(paste("Saving plot to ", file, ".png", 
              sep = ""))
  p <- grid.arrange(arrangeGrob(p1, top = toptitle), arrangeGrob(p2, 
                                                                 bottom = bottomtitle), padding = 0, heights = c(hgtratio, 
                                                                                                                 1 - hgtratio))
  ggsave(p, filename = paste(file, ".png", sep = ""), 
         dpi = res, units = "in", height = hgt, width = wi)
  return(p)
}

data(probe.features.epic)

female <- read.csv("TEST STATISTICS FEMALES.csv", sep=",", dec=".", header=TRUE, row.names=1)
male <- read.csv("TEST STATISTICS MALES.csv", sep=",", dec=".", header=TRUE, row.names=1)

probe.features <- probe.features[row.names(probe.features) %in% row.names(female), ]
female <- female[order(row.names(female)), ]
male <- male[order(row.names(male)), ]

female_plot <- data.frame(row.names(female), probe.features$CHR, probe.features$MAPINFO, female$P)
colnames(female_plot) <- c("SNP", "CHR", "POS", "pvalue")
male_plot <- data.frame(row.names(male), probe.features$CHR, probe.features$MAPINFO, male$P)
colnames(male_plot) <- c("SNP", "CHR", "POS", "pvalue")

female_plot <- female_plot[with(female_plot, order(female_plot$CHR, female_plot$POS)), ]
male_plot <- male_plot[with(male_plot, order(male_plot$CHR, male_plot$POS)), ]

#If there are epigenome-wide significant hits, annotate them either with the CpG-ID or gene
#You can set here whether you would like to annotate genes or cpgs. For genes, set "cpgannot" to FALSE, and "geneannot" to TRUE; for CpGs vice versa
geneannot = TRUE
cpgannot = FALSE

#Specify the name your generated png files should have here
filename = "FILENAME"

#Specify the P-value threshold for annotations here; here set to Mansell et al. epigenome wide singificance threshold
p_threshold <- 8.62E-9

if(length(which(female$P < p_threshold)>0) | length(which(male$P < p_threshold)>0)){
  if(cpgannot){
    gmirror(file = paste(filename, "cpgannot", sep="_"), top=female_plot, bottom=male_plot, toptitle="Females", bottomtitle = "Males", highlight_p = c(p_threshold, p_threshold), annotate_p = c(p_threshold, p_threshold))
  } else{
    female_plot <- data.frame(probe.features$gene, probe.features$CHR, probe.features$MAPINFO, female$P)
    colnames(female_plot) <- c("SNP", "CHR", "POS", "pvalue")
    male_plot <- data.frame(probe.features$gene, probe.features$CHR, probe.features$MAPINFO, male$P)
    colnames(male_plot) <- c("SNP", "CHR", "POS", "pvalue")

    female_plot <- female_plot[with(female_plot, order(female_plot$CHR, female_plot$POS)), ]
    male_plot <- male_plot[with(male_plot, order(male_plot$CHR, male_plot$POS)), ]

    gmirror(file = paste(filename, "geneannot", sep="_"), top=female_plot, bottom=male_plot, toptitle="Females", bottomtitle = "Males", highlight_p = c(p_threshold, p_threshold), annotate_p = c(p_threshold, p_threshold))
  }
  
} else{
  gmirror(file = filename, top=female_plot, bottom=male_plot, toptitle="Females", bottomtitle = "Males", highlight_p = c(p_threshold, p_threshold))
}
