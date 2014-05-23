library(edgeR)
library(ggplot2)
library(plyr)
library(reshape2)
library(testthat) # facilitate tests that will catch changes on re-analysis

### Differential Expression Analysis on Sitka Spruce Weevil 
### Experiment with limma + voom

# Load counts from Sailfish
filteredSailfishCounts <- read.delim("consolidated-filtered-Sailfish-results.txt")
str(filteredSailfishCounts) # 'data.frame':  65609 obs. of  24 variables:
test_that("filtered Sailfish data has 65609 rows upon import",
          expect_equal(65609, nrow(filteredSailfishCounts)))
test_that("Sailfish data has data for exactly 24 samples",
          expect_equal(24, ncol(filteredSailfishCounts)))

# Load design matrix
desMat <- read.delim("White_Pine_Weevil_design_matrix.tsv",
                     stringsAsFactors = FALSE)
desMat <-
  mutate(desMat,
         gTypeCode = factor(gTypeCode, levels = c("Q903", "H898")),
         gType = factor(gType, levels = c("Susc", "Res")),
         txCode = factor(txCode, levels = c('C', 'W', 'G')),
         tx = factor(tx, levels = c("Control", "Wound", "Gallery")))
desMat$grp <-
  with(desMat, factor(grp,
                      levels = paste(levels(gTypeCode),
                                     rep(levels(tx), each = 2), sep = ".")))
str(desMat) # 'data.frame':  24 obs. of  7 variables:
test_that("design matrix has 24 rows upon import", expect_equal(24, nrow(desMat)))

# Load counts into DGEList object from edgeR package.
y <- DGEList(counts = filteredSailfishCounts, group = desMat$grp)

# TMM Normalization by Depth
y <- calcNormFactors(y)

# make model matrix
#modMat <- model.matrix(~ gType/tx - 1, desMat)
modMat <- model.matrix(~ tx/gType - 1, desMat)
modMat <- model.matrix(~ gType * tx, desMat)

# voom transformation
v <- voom(y, modMat, plot = TRUE)

# Linear modelling
fit <- lmFit(v, modMat)
fit <- eBayes(fit)

tt <- topTable(fit, coef = grep(":", colnames(modMat)))

data.wide <- cbind(desMat, t(y$counts[rownames(tt)[1:4], ]))
data.tall <- melt(data.wide,
                  id.vars = names(desMat),
                  variable.name = 'contig', value.name = 'Expression')
str(data.tall)

#x <- data.frame(desMat, gExp = y$counts[rownames(tt)[1], ])
p <- ggplot(data.tall, aes(x = gType, y = Expression))
p + geom_point() + facet_grid(tx ~ contig)

# Linear design matrix, simplest approach
design <- model.matrix(~ 0 + grp, desMat)
colnames(design) <- levels(desMat$grp)

# Create contrast matrix
cont.matrix <- makeContrasts(
  Control.H898_vs_Q903 = H898.Control - Q903.Control,
  Induced.H898_vs_Q903 = 
    (H898.Gallery - (H898.Wound - H898.Control)) - 
    (Q903.Gallery - (Q903.Wound - Q903.Control)),
  levels = design)

# Voom transformation
v <- voom(y, modMat, plot = TRUE)



# MDS analysis
plotMDS(v, top=Inf)

# PCA analysis
pca <- as.data.frame(prcomp(t(v$E))$x)
ggplot(pca, aes(PC1, PC2, color = Group)) + geom_point(size=3) + 
  scale_color_manual(name = "", 
                     values=c("H898.Control" = "#4d004b", 
                              "H898.Gallery" = "#88419d", 
                              "H898.Wound" = "#8c96c6", 
                              "Q903.Control" = "#7F0000", 
                              "Q903.Gallery" =  "#d7301f", 
                              "Q903.Wound" = "#fc8d59")) + 
  labs(title="Voom + Limma Principal Component Analysis") + theme_bw()

# Linear modelling
fit <- lmFit(v, design)

fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2)

linear.results <- decideTests(fit2, method="global", 
                              adjust.method="fdr", p.value=0.01)

summary(linear.results)

###

vennDiagram(linear.results, include="both")
vennDiagram(linear.results, include="up")
vennDiagram(linear.results, include="down")

###
