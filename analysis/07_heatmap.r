library(NMF)
library(vegan)
library(edgeR)

setwd("analysis/")

source("helper03_load-focus-statinf.r")

lfc <- 2
pCutoff <- 0.01

rawSailfishCounts <- read.delim("../data/consolidated-Sailfish-results.txt")
str(rawSailfishCounts) # 'data.frame':  483047 obs. of  24 variables:

cpmCounts <- (cpm(rawSailfishCounts))
str(cpmCounts)

avgCPM <- data.frame("H898C" = rowMeans(cpmCounts[,1:4]),
                     "H898G" = rowMeans(cpmCounts[,5:8]),
                     "H898W" = rowMeans(cpmCounts[,9:12]),
                     "Q903C" = rowMeans(cpmCounts[,13:16]),
                     "Q903G" = rowMeans(cpmCounts[,17:20]),
                     "Q903W" = rowMeans(cpmCounts[,21:24]))
str(avgCPM) #483047 obs. of  6 variables

###

sidf <- load_focus_statInf()
str(sidf)


### Weevil Induce Q903 (Q903G - Q903W)

weevilInd_q903 <- (subset(sidf, sidf$focus_term == "weevilInd_Q903" & 
                            sidf$adj.P.Val <= pCutoff & abs(sidf$logFC) >= lfc))
dim(weevilInd_q903)
# [1] 5802    8

weevilInd_q903_avgCPM <- avgCPM[rownames(avgCPM) %in% weevilInd_q903$contig,]
dim(weevilInd_q903_avgCPM)
# [1] 5802    8

data.dist <- vegdist(weevilInd_q903_avgCPM, method = "jaccard")

row.clus <- hclust(data.dist, "ward.D2")

aheatmap(as.matrix(weevilInd_q903_avgCPM), breaks = 0,
         border_color = NA, cellwidth = 50,
         scale = "row",
         color = "heat", 
         Rowv = row.clus,
         Colv = TRUE,
         width = 8.5, height = 11, 
         filename = "../results/figures/weevilInd_Q903_heatmap.png")


### Constitutive Difference (H898C - Q903C)

constDiff <- (subset(sidf, sidf$focus_term == "constDiff" & 
                       sidf$adj.P.Val <= pCutoff & abs(sidf$logFC) >= lfc))
dim(constDiff)
# [1] 7012    8

constDiff_avgCPM <- avgCPM[rownames(avgCPM) %in% constDiff$contig,]
dim(constDiff_avgCPM)
# [1] 7012    6

data.dist <- vegdist(constDiff_avgCPM, method = "jaccard")

row.clus <- hclust(data.dist, "ward.D2")

aheatmap(as.matrix(constDiff_avgCPM), breaks = 0,
         border_color = NA, cellwidth = 50,
         scale = "row",
         color = "heat", 
         Rowv = row.clus,
         Colv = TRUE,
         width = 8.5, height = 11, 
         filename = "../results/figures/constDiff_heatmap.png")

