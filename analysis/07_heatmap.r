# library(gplots)
# library(ggplot2)
# library(plyr)
# library(RColorBrewer)
# library(dplyr)


setwd("analysis/")

source("helper03_load-focus-statinf.r")

lfc <- 2
pCutoff <- 0.01


rawSailfishCounts <- read.delim("../data/consolidated-Sailfish-results.txt")
str(rawSailfishCounts) # 'data.frame':  483047 obs. of  24 variables:

avgSailfishCounts <- data.frame("H898C" = rowMeans(rawSailfishCounts[,1:4]),
                                "H898G" = rowMeans(rawSailfishCounts[,5:8]),
                                "H898W" = rowMeans(rawSailfishCounts[,9:12]),
                                "Q903C" = rowMeans(rawSailfishCounts[,13:16]),
                                "Q903G" = rowMeans(rawSailfishCounts[,17:20]),
                                "Q903W" = rowMeans(rawSailfishCounts[,21:24]))
                                
###

sidf <- load_focus_statInf()
str(sidf)

###

library(NMF)

tmp <- (subset(sidf, sidf$focus_term == "weevilInd_Q903" & 
                 sidf$adj.P.Val <= 0.01 & abs(sidf$logFC) >= 2))
# tmp.sorted <- tmp[order(tmp$logFC, decreasing = TRUE),]
tmp.sorted <- tmp[order(tmp$AveExpr, decreasing = TRUE),]
testData <- avgSailfishCounts[tmp.sorted$contig,]

# testData <- avgSailfishCounts[
#   (subset(sidf, sidf$focus_term == "weevilInd_Q903" & 
#            sidf$adj.P.Val <= 0.01 & abs(sidf$logFC) >= 2)[,"contig"]),]
  
testMtrx <- as.matrix(testData)
 
aheatmap(testMtrx, breaks = 0,
         border_color = NA, cellwidth = 50,
         scale = "row",
         # Rowv = NA,
         Colv = NA,
         color = "heat", 
         distfun = "euclidean",
         hclustfun = "ward",
         width = 8.5, height = 11, 
         filename = "../results/test_heatmap.png")

# distfun = "maximum", hclustfun = "ward",         
# , Rowv = 200)


# aheatmap(iris2, color = "-RdBu:50", scale = "col", breaks = 0,
#          annRow = iris["Species"], annColors = "Set2", 
#          distfun = "pearson", treeheight=c(200, 50), 
#          fontsize=13, cexCol=.7, 
#          filename="heatmap.png", width=8, height=16)
###

sd_constDiff <- 
  subset(sidf_wide, (sidf_wide$adj.P.Val.constDiff <= pCutoff & 
                       abs(sidf_wide$logFC.constDiff) >= lfc))

sd_constDiff_mtrx <- as.matrix(sd_constDiff[,c(1, 3, 5, 7, 9, 11, 13)])


aheatmap(sd_constDiff_mtrx, color = "heat", scale = "row", hclustfun = "average", 
         distfun = "manhattan", width = 8.5, height = 11, 
         cellwidth = 25, filename = "../results/test_heatmap.png", Rowv = 200)

###

sd_weevilIndQ903 <- 
  subset(sidf_wide, (sidf_wide$adj.P.Val.weevilInd_Q903 <= pCutoff & 
                       abs(sidf_wide$logFC.weevilInd_Q903) >= lfc))

sd_weevilIndQ903_mtrx <- as.matrix(sd_weevilIndQ903[,c(1, 3, 5, 7, 9, 11, 13)])

aheatmap(sd_weevilIndQ903_mtrx, color = "heat", scale = "row", hclustfun = "average", 
         distfun = "manhattan", width = 8.5, height = 11, 
         cellwidth = 25, filename = "../results/test_heatmap.png", Rowv = 200)

#####

sidf_wide_mtrx <- as.matrix(sidf_wide[,c(1, 3, 5, 7, 9, 11, 13)])


#####

iris2 = iris # prep iris data for plotting

rownames(iris2) = make.names(iris2$Species, unique = T)

iris2 = iris2 %>% select(-Species) %>% as.matrix()

aheatmap(iris2, color = "-RdBu:50", scale = "col", breaks = 0,
         annRow = iris["Species"], annColors = "Set2", 
         distfun = "pearson", treeheight=c(200, 50), 
         fontsize=13, cexCol=.7)


#####

library(Heatplus)
library(vegan)

all.data <- 
  read.delim("~/Documents/UBC_Projects/MPB_Genomics/Transcriptomics/Tria II/2014-analyses/R_analyses/CLC-Assembly-quantification/CompleteDesign/Heatmap.txt",
             header = TRUE, sep = "\t", row.names = "ID")  # load the data

dim(all.data)

data.dist <- vegdist(all.data, method = "jaccard")

row.clus <- hclust(data.dist, "ward.D2")

stat.data <- 
  read.delim("~/Documents/UBC_Projects/MPB_Genomics/Transcriptomics/Tria II/2014-analyses/R_analyses/CLC-Assembly-quantification/CompleteDesign/Stats_data.txt",
             header=TRUE, sep="\t", row.names="ID")  # load the data

anna1 = convAnnData(stat.data)

colnames(anna1) = 
  c("Sex", "Treatment", "Tissue", "Sex x Treatment", "Sex x Tissue",
    "Treatment x Tissue", "Sex x Treatment x Tissue", "MMG", "MFB",
    "FMG", "FFB", "Expression")

# pdf("~/Desktop/Annotated-Heat-Map.pdf", width=8, height=395)

par(pty="m")

plot(
  annHeatmap2(
    as.matrix(all.data), scale="row", 
    col = colorRampPalette(c("#FFFFCC", "#a1dab4", "#41b6c4", "#2c7fb8", "#253494"), 
                           space = "rgb")(96), breaks = 95, 
    dendrogram = list( Row = list(adj = 0,dendro = as.dendrogram(row.clus)),
                       Col=list(status="no")), legend = 3, 
    ann = list(asIs=TRUE, Row = list(data = anna1))), widths=c(1,3,3), heights=c(1,340,1))

# dev.off()

sessionInfo()

