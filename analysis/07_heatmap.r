library(edgeR)
library(vegan)
library(d3heatmap)


setwd("analysis/")
source("helper03_load-focus-statinf.r")


lfc <- 2
pCutoff <- 0.01


rawSailfishCounts <- read.delim("../data/consolidated-Sailfish-results.txt")
str(rawSailfishCounts) # 'data.frame':  483047 obs. of  24 variables:

cpmCounts <- cpm(rawSailfishCounts)
str(cpmCounts)

avgCPM <- data.frame("H898C" = rowMeans(cpmCounts[,1:4]),
                     "H898G" = rowMeans(cpmCounts[,5:8]),
                     "H898W" = rowMeans(cpmCounts[,9:12]),
                     "Q903C" = rowMeans(cpmCounts[,13:16]),
                     "Q903G" = rowMeans(cpmCounts[,17:20]),
                     "Q903W" = rowMeans(cpmCounts[,21:24]))
str(avgCPM) #483047 obs. of  6 variables


sidf <- load_focus_statInf()
str(sidf) # 406728 obs.


annots <- read.delim(
  "../results/WPW_Inoculation_Trinity_C500.diffExp.lfc2.blastxNr.collapsed.txt", row.names = 1)
annots$annot <- as.character(annots$annot)
dim(annots) # [1] 12852     2


### Weevil Induce Q903 (Q903G - Q903W)
weevilInd_q903 <- (subset(sidf, sidf$focus_term == "weevilInd_Q903" & 
                            sidf$adj.P.Val <= pCutoff & abs(sidf$logFC) >= lfc))
dim(weevilInd_q903) # [1] 5802    8


# Subset weevilInd_Q903 DE contigs from avgCPM
weevilInd_q903_avgCPM <- avgCPM[rownames(avgCPM) %in% weevilInd_q903$contig,]
dim(weevilInd_q903_avgCPM) # [1] 5802    6


# Create cellnote when mouse hover over
weevilInd_q903_cellnote <- weevilInd_q903_avgCPM

for (i in 1:6) {
  weevilInd_q903_cellnote[, i] <- 
    paste("CPM : ", weevilInd_q903_cellnote[, i], "\n",
          "Annot : ", annots[rownames(weevilInd_q903_cellnote), "annot"])
}


weevilInd_q903_avgCPM_dist <- vegdist(weevilInd_q903_avgCPM, method = "horn")
weevilInd_q903_avgCPM_hclust <- hclust(weevilInd_q903_avgCPM_dist, "average")
weevilInd_q903_avgCPM_dendrogram <- as.dendrogram(weevilInd_q903_avgCPM_hclust)

weevilInd_q903_col_dist <- vegdist(t(weevilInd_q903_avgCPM), method = "horn")
weevilInd_q903_col_hclust <- hclust(weevilInd_q903_col_dist, "average")
weevilInd_q903_col_dendrogram <- as.dendrogram(weevilInd_q903_col_hclust)


d3heatmap(as.matrix(weevilInd_q903_avgCPM), 
          colors = heat.colors(50),
          scale = "row",
          Rowv = weevilInd_q903_avgCPM_dendrogram,
          Colv = weevilInd_q903_col_dendrogram,
          height = 8000,
          cellnote = weevilInd_q903_cellnote,
          anim_duration = 0)


### Constitutive Difference (H898C - Q903C)
constDiff <- (subset(sidf, sidf$focus_term == "constDiff" & 
                       sidf$adj.P.Val <= pCutoff & abs(sidf$logFC) >= lfc))
dim(constDiff)
# [1] 7012    8


# Subset constDiff DE contigs from avgCPM
constDiff_avgCPM <- avgCPM[rownames(avgCPM) %in% constDiff$contig,]
dim(constDiff_avgCPM)
# [1] 7012    6


# Create cellnote when mouse hover over
constDiff_avgCPM_cellnote <- constDiff_avgCPM

for (i in 1:6) {
  constDiff_avgCPM_cellnote[, i] <- 
    paste("CPM : ", constDiff_avgCPM_cellnote[, i], "\n",
          "Annot : ", annots[rownames(constDiff_avgCPM_cellnote), "annot"])
}


constDiff_avgCPM_dist <- vegdist(constDiff_avgCPM, method = "horn")
constDiff_avgCPM_hclust <- hclust(constDiff_avgCPM_dist, "average")
constDiff_avgCPM_dendrogram <- as.dendrogram(constDiff_avgCPM_hclust)

constDiff_col_dist <- vegdist(t(constDiff_avgCPM), method = "horn")
constDiff_col_hclust<- hclust(constDiff_col_dist, "average")
constDiff_col_dendrogram <- as.dendrogram(constDiff_col_hclust)


d3heatmap(as.matrix(constDiff_avgCPM), 
          colors = heat.colors(50),
          scale = "row", 
          Rowv = constDiff_avgCPM_dendrogram,
          Colv = constDiff_col_dendrogram,
          height = 8000,
          cellnote = constDiff_avgCPM_cellnote,          
          anim_duration = 0)


####


identify(dendrogram)
plot(weevilInd_q903_avgCPM_hclust, hang = -1)


plot(weevilInd_q903_avgCPM_hclust, hang = -1, axes = FALSE)
axis(side = 2, at = seq(0, 0.6, 0.05), col = "#F38630", labels = FALSE, lwd = 2)
mtext(seq(0, 0.6, 0.05), side = 2, at = seq(0, 0.6, 0.05), line = 1, col = "#A38630", las = 2)

branches <- cutree(weevilInd_q903_avgCPM_hclust, h = 0.1)


# Find out the branch from genes identified from heatmap
branches["comp435435_c0_seq4"]

# Extract branch
myBranch <- branches[branches == 11]
length(myBranch)

write.table(names(myBranch), file = "../results/branch.txt", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)
