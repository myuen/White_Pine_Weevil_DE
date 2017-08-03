library(edgeR)
library(ggplot2)
library(testthat) # facilitate tests that will catch changes on re-analysis

### Differential Expression Analysis on Sitka Spruce Weevil 
### Experiment with limma + voom

#' Source purpose-built functions to load and validate the data and the 
#' experimental design. Then call them.
setwd("analysis/")
source("helper01_load-counts.r")
source("helper02_load-exp-des.r")

#' Load counts from Sailfish
x <- load_counts() # takes a few moments
str(x, list.len = 8) # 'data.frame':	38197 obs. of  24 variables:

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
v <- voom(y, modMat, plot = FALSE) # take a couple moments


# PCA analysis
pca <- as.data.frame(prcomp(t(v$E))$x)

wpw_pca <- 
  ggplot(pca, aes(PC1, PC2, color = expDes$grp, shape = expDes$gType)) + geom_point(size = 2.5) +
  scale_color_manual(
    limits = c("H898.Control", "H898.Gallery", "H898.Wound",
               "Q903.Control", "Q903.Gallery", "Q903.Wound"),
    values=c("H898.Control" = "#4d004b",
             "H898.Gallery" = "#88419d",
             "H898.Wound" = "#8c96c6",
             "Q903.Control" = "#7F0000",
             "Q903.Gallery" = "#d7301f",
             "Q903.Wound" = "#fc8d59")) +
  labs(color = "Groups", shape = NULL) +
  theme_bw() +
  theme(legend.position = "bottom", legend.box.background = element_rect(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(shape = 15)), shape = FALSE)

ggsave("../results/figures/fig6a.3Aug2017.tiff", plot = wpw_pca, width = 9, height = 7)
