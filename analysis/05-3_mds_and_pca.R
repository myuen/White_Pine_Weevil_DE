library(edgeR)
library(ggplot2)
library(testthat) # facilitate tests that will catch changes on re-analysis

# opts_chunk$set(fig.path = 'figure/03-dea-with-limma-voom-')

### Differential Expression Analysis on Sitka Spruce Weevil 
### Experiment with limma + voom

#' Source purpose-built functions to load and validate the data and the 
#' experimental design. Then call them.
setwd("analysis/")
source("helper01_load-counts.r")
source("helper02_load-exp-des.r")

#' Load counts from Sailfish
x <- load_counts() # takes a few moments
str(x, list.len = 8) # 'data.frame':  65609 obs. of  24 variables:

#' Load experimental design
expDes <- load_expDes()
expDes$grp <- # not really sure that we need this?
  with(expDes, factor(grp,
                      levels = paste(levels(gType),
                                     rep(levels(tx), each = 2), sep = ".")))
str(expDes) # 'data.frame':  24 obs. of  6 variables:

#' Load counts into DGEList object from edgeR package.
y <- DGEList(counts = x, group = expDes$grp)

#' TMM Normalization by Depth
y <- calcNormFactors(y)

#' make model matrix
modMat <- model.matrix(~ gType * tx, expDes)

colnames(modMat)
colnames(modMat) <- gsub(":", "_", colnames(modMat))
colnames(modMat) <- gsub("[()]", "", colnames(modMat))
colnames(modMat)

#' voom transformation
#+ voom-plot
v <- voom(y, modMat, plot = TRUE) # take a couple moments


# MDS analysis
png(filename = "../results/figures/wpw_mds.png", width = 1680, height = 1050)
plotMDS(v, top=Inf)
dev.off()


# PCA analysis
pca <- as.data.frame(prcomp(t(v$E))$x)

wpw_pca <- ggplot(pca, aes(PC1, PC2, color = expDes$grp, shape = expDes$gType)) + 
  geom_point(size = 5) + scale_color_manual(
    limits = c("H898res.Control", "H898res.Gallery", "H898res.Wound",
               "Q903susc.Control", "Q903susc.Gallery", "Q903susc.Wound"),
    values=c("H898res.Control" = "#4d004b",
             "H898res.Gallery" = "#88419d",
             "H898res.Wound" = "#8c96c6",
             "Q903susc.Control" = "#7F0000",
             "Q903susc.Gallery" = "#d7301f",
             "Q903susc.Wound" = "#fc8d59")) +
  labs(title="Voom + Limma Principal Component Analysis") + theme_bw()

ggsave("../results/figures/wpw_pca.png", plot = wpw_pca, width = 16, height = 9)
