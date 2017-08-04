### Figure 6e & 6f
### Plotting heatmap of only manually curated contigs in DEA

library(d3heatmap)
library(htmlwidgets)
library(plyr)

setwd("analysis/")
source("helper02_load-exp-des.r")
source("helper03_load-focus-statinf.r")


# abs(logFC) log fold change cut-off.  Anything greater than (-1 x lfc) and less 
# than lfc will be deemed biological insignificant
lfcCutoff <- 2

# p-value cut-off.  Anything > pCutoff will be deemed statistically insignificant.
pCutoff <- 0.05


cpmCounts <- read.table("../results/normalized_cpm.15July.txt")


expDes <- load_expDes()
expDes$grp <-
  with(expDes, factor(grp,
                      levels = paste(levels(gType),
                                     rep(levels(tx), each = 2), sep = ".")))

avgCPM <- ddply(expDes, ~grp, function(x) {
  df <- rowMeans(cpmCounts[, c(expDes[expDes$grp == x$grp,]$sample)])
  return(df)
})

# Shorten contig name for easy display
colnames(avgCPM) <- 
  gsub("WPW_Inoculation_Trinity_BBMerged_31May2016_TRINITY_", "", colnames(avgCPM))
grp <- as.vector(avgCPM$grp)

# Rename group name to be consistent with manuscript
# S = susceptible, Q903
# R = resistance, H898
grp <- gsub("Q903", "S", grp)
grp <- gsub("H898", "R", grp)
grp <- gsub("Control", "C", grp)
grp <- gsub("Wound", "W", grp)
grp <- gsub("Gallery", "G", grp)
grp <- gsub("\\.", "-", grp)

avgCPM <- as.data.frame(t(avgCPM[,-1]))
colnames(avgCPM) <- grp


sidf <- load_focus_statInf()
sidf$cds <- 
  as.character(gsub("WPW_Inoculation_Trinity_BBMerged_31May2016_TRINITY_", "", sidf$cds))


annots <- 
  read.delim("../results/WPW_Inoculation_Trinity_BBMerged_31May2016.transdecoder.cds.nr.diffExp.blastxNR.collapsed.txt", 
  row.names = 1)
row.names(annots) <- 
  gsub("WPW_Inoculation_Trinity_BBMerged_31May2016_TRINITY_", "", row.names(annots))
annots$annot <- as.character(annots$annot)


# Manual curated genes list
manAnn <- scan("../data/manual_annotations/allManuallyAnnotated.contigId", what = "string")
manAnn <- gsub("WPW_Inoculation_Trinity_BBMerged_31May2016_TRINITY_", "", manAnn)


# Constitutive Control (H898C - Q903C)
constDiff <- subset(sidf, sidf$focus_term == "constDiff" & 
                      sidf$adj.P.Val <= pCutoff & abs(sidf$logFC) >= lfcCutoff &
                      sidf$cds %in% manAnn)

# Subset constitutive DE contigs from avgCPM
constDiff_avgCPM <- avgCPM[rownames(avgCPM) %in% constDiff$cds,]


# Create cellnote when mouse hover over
constDiff_cellnote <- constDiff_avgCPM


for (i in 1:6) {
  constDiff_cellnote[, i] <-
    paste("CPM : ", constDiff_cellnote[, i], "\n",
          "Annot : ", annots[rownames(constDiff_cellnote), "annot"])
}

constDiff_row_dist <- as.dist(1 - cor(t(constDiff_avgCPM), method = "spearman"))
constDiff_row_hclust <- hclust(constDiff_row_dist, "complete")
constDiff_row_dendrogram <- as.dendrogram(constDiff_row_hclust)

constDiff_col_dist <- as.dist(1 - cor(constDiff_avgCPM, method = "spearman"))
constDiff_col_hclust <- hclust(constDiff_col_dist, "complete")
constDiff_col_dendrogram <- as.dendrogram(constDiff_col_hclust)


cdHM <- d3heatmap(as.matrix(constDiff_avgCPM),
                  colors = heat.colors(50),
                  scale = "row",
                  Rowv = constDiff_row_dendrogram,
                  Colv = constDiff_col_dendrogram,
                  cexCol = 0.5,
                  cellnote = constDiff_cellnote,
                  anim_duration = 0)
saveWidget(cdHM, "fig6e.3Aug.html", selfcontained = TRUE)


# Weevil Induce Q903 (Q903G - Q903W)
weevilInd_q903 <- subset(sidf, sidf$focus_term == "weevilInd_Q903" & 
                           sidf$adj.P.Val <= pCutoff & abs(sidf$logFC) >= lfcCutoff &
                           sidf$cds %in% manAnn)

# Subset weevilInd_Q903 DE contigs from avgCPM
weevilInd_q903_avgCPM <- avgCPM[rownames(avgCPM) %in% weevilInd_q903$cds,]


# Create cellnote when mouse hover over
weevilInd_q903_cellnote <- weevilInd_q903_avgCPM


for (i in 1:6) {
  weevilInd_q903_cellnote[, i] <-
    paste("CPM : ", weevilInd_q903_cellnote[, i], "\n",
          "Annot : ", annots[rownames(weevilInd_q903_cellnote), "annot"])
}

weevilInd_q903_row_dist <- as.dist(1 - cor(t(weevilInd_q903_avgCPM), method = "spearman"))
weevilInd_q903_row_hclust <- hclust(weevilInd_q903_row_dist, "complete")
weevilInd_q903_row_dendrogram <- as.dendrogram(weevilInd_q903_row_hclust)

weevilInd_q903_col_dist <- as.dist(1 - cor(weevilInd_q903_avgCPM, method = "spearman"))
weevilInd_q903_col_hclust <- hclust(weevilInd_q903_col_dist, "complete")
weevilInd_q903_col_dendrogram <- as.dendrogram(weevilInd_q903_col_hclust)


wiqHM <- d3heatmap(as.matrix(weevilInd_q903_avgCPM),
          colors = heat.colors(50),
          scale = "row",
          Rowv = weevilInd_q903_row_dendrogram,
          Colv = weevilInd_q903_col_dendrogram,
          cexCol = 0.5,
          cellnote = weevilInd_q903_cellnote,
          anim_duration = 0)
saveWidget(wiqHM, "fig6f.3Aug.html", selfcontained = TRUE)

# saveWidget don't handle well with relative path.  Save in "analysis" folder and move
# over to ../results/figure" manually

# Since d3heatmap does not save as image, have to resort to "Export" function in RStudio
# to save as PNG.
