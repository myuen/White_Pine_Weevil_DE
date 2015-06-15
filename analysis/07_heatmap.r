library(edgeR)
library(grid)
library(vegan)
library(NMF)
library(RColorBrewer)


# Putative Phenylpropanoid Pathway Enzymes
ppid <- scan("data/putativePhenylpropanoidEnz.id", what = "character")
length(ppid)
# [1] 322

lid <- scan("data/putativeLignanBiosynEnzymes.id", what = "character")
length(lid)
# [1] 68

mevmepid <- scan("data/putativeMevMepEnzymes.id", what = "character")
length(mevmepid)
# [1] 125

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

###

sidf <- load_focus_statInf()
str(sidf)

### Weevil Induce Q903 (Q903G - Q903W)
weevilInd_q903 <- (subset(sidf, sidf$focus_term == "weevilInd_Q903" & 
                            sidf$adj.P.Val <= pCutoff & abs(sidf$logFC) >= lfc))
dim(weevilInd_q903)
# [1] 5802    8


table(weevilInd_q903$contig %in% ppid)
# FALSE  TRUE 
# 5730    72 

table(weevilInd_q903$contig %in% lid)
# FALSE  TRUE 
# 5783    19 

table(weevilInd_q903$contig %in% mevmepid)
# FALSE  TRUE 
# 5797     5 

weevilInd_q903[weevilInd_q903$contig %in% ppid, "Annotation"] <- "Phenylpropanoid Biosyn"
weevilInd_q903[weevilInd_q903$contig %in% lid, "Annotation"] <- "Lignan Biosyn"
weevilInd_q903[weevilInd_q903$contig %in% mevmepid, "Annotation"] <- "MEV/MEP"

table(weevilInd_q903$Annotation)
#          Lignan Biosyn                MEV/MEP                    N/A Phenylpropanoid Biosyn 
#                     19                      5                   5722                     56 

weevilInd_q903[is.na(weevilInd_q903$Annotation), "RowLab"] <- ""

weevilInd_q903[!is.na(weevilInd_q903$Annotation), "RowLab"] <- 
  gsub("WPW_Inoculation_Trinity_C500_", "", 
       weevilInd_q903[!is.na(weevilInd_q903$Annotation), "contig"])


# FOR DEBUG
# head(weevilInd_q903[!is.na(weevilInd_q903$Annotation),])


weevilInd_q903_avgCPM <- avgCPM[rownames(avgCPM) %in% weevilInd_q903$contig,]
dim(weevilInd_q903_avgCPM)
# [1] 5802    6


# Creating dendrogram with just the 2 libraries
# data.dist <- vegdist(weevilInd_q903_avgCPM, method = "bray")
weevilInd_q903_avgCPM_dist <- 
  vegdist(weevilInd_q903_avgCPM[,c("Q903G", "Q903W")], method = "jaccard")

weevilInd_q903_avgCPM_clusters <- hclust(weevilInd_q903_avgCPM_dist, "ward.D2")

weevilInd_q903_avgCPM_dendrogram <- as.dendrogram(weevilInd_q903_avgCPM_clusters)

weevilInd_q903_avgCPM <- as.matrix(weevilInd_q903_avgCPM)

aheatmap(weevilInd_q903_avgCPM, 
         color = "heat", breaks = 0,
         border_color = NA, cellwidth = 40, 
         scale = "row", Rowv = weevilInd_q903_avgCPM_clusters,
         treeheight = c(200, 50), legend = TRUE,
         annRow = weevilInd_q903$Annotation, annColors = "Set1", annLegend = TRUE,
         labRow = weevilInd_q903$RowLab,
         width = 20, height = 20,
         main = "Q903 Weevil Induce (Q903G - Q903W) Comparison", 
         sub = "Libraries", info = TRUE,
         gp = gpar(col = "black", fontsize = 6, cexCol = 0.7, cexRow = 0.5),
         filename = "../results/figures/weevilInd_Q903_annot_heatmap.png")


### Constitutive Difference (H898C - Q903C)
constDiff <- (subset(sidf, sidf$focus_term == "constDiff" & 
                       sidf$adj.P.Val <= pCutoff & abs(sidf$logFC) >= lfc))
dim(constDiff)
# [1] 7012    8

table(constDiff$contig %in% ppid)
# FALSE  TRUE 
#  6985    27 

table(constDiff$contig %in% lid)
# FALSE  TRUE 
#  7005     7 

table(constDiff$contig %in% mevmepid)
# FALSE  TRUE 
#  7008     4 

constDiff[constDiff$contig %in% ppid, "Annotation"] <- "Phenylpropanoid Biosyn"
constDiff[constDiff$contig %in% lid, "Annotation"] <- "Lignan Biosyn"
constDiff[constDiff$contig %in% mevmepid, "Annotation"] <- "MEV/MEP"

table(constDiff$Annotation)
#          Lignan Biosyn                MEV/MEP                    N/A Phenylpropanoid Biosyn 
#                     19                      5                   5722                     56 

constDiff[is.na(constDiff$Annotation), "RowLab"] <- ""

constDiff[!is.na(constDiff$Annotation), "RowLab"] <- 
  gsub("WPW_Inoculation_Trinity_C500_", "", 
       constDiff[!is.na(constDiff$Annotation), "contig"])

# FOR DEBUG
# head(constDiff[!is.na(constDiff$Annotation),])


constDiff_avgCPM <- avgCPM[rownames(avgCPM) %in% constDiff$contig,]
dim(constDiff_avgCPM)
# [1] 7012    6

# constDiff_avgCPM_dist <- vegdist(constDiff_avgCPM, method = "jaccard")
constDiff_avgCPM_dist <- 
  vegdist(constDiff_avgCPM[,c("H898C", "Q903C")], method = "jaccard")

constDiff_avgCPM_clusters <- hclust(data.dist, "ward.D2")

constDiff_avgCPM <- as.matrix(constDiff_avgCPM)

aheatmap(constDiff_avgCPM, 
         color = "heat", breaks = 0,
         border_color = NA, cellwidth = 40, 
         scale = "row", Rowv = constDiff_avgCPM_clusters,
         treeheight = c(200, 50), legend = TRUE,
         annRow = constDiff$Annotation, annColors = "Set1", annLegend = TRUE,
         labRow = constDiff$RowLab,
         width = 20, height = 20,
         main = "Constitutive Difference (H898C - Q903C) Comparison", 
         sub = "Libraries", info = TRUE,
         gp = gpar(col = "black", fontsize = 6, cexCol = 0.7, cexRow = 0.5),
         filename = "../results/figures/constDiff_annot_heatmap.png")