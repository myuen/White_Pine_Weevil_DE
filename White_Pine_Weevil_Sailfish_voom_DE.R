library(edgeR)
library(ggplot2)
library(plyr)
library(reshape2)
library(testthat) # facilitate tests that will catch changes on re-analysis

### Differential Expression Analysis on Sitka Spruce Weevil 
### Experiment with limma + voom

# Load counts from Sailfish
filteredSailfishCounts <- # take a few moments
  read.delim("consolidated-filtered-Sailfish-results.txt")
str(filteredSailfishCounts) # 'data.frame':  65609 obs. of  24 variables:
test_that("filtered Sailfish data has 65609 rows upon import",
          expect_equal(65609, nrow(filteredSailfishCounts)))
test_that("Sailfish data has data for exactly 24 samples",
          expect_equal(24, ncol(filteredSailfishCounts)))

# Load experimental design
expDes <- read.delim("White_Pine_Weevil_exp_design.tsv",
                     stringsAsFactors = FALSE)
expDes <-
  mutate(expDes,
         gType = factor(gType, levels = c("Q903susc", "H898res")),
         txCode = factor(txCode, levels = c('C', 'W', 'G')),
         tx = factor(tx, levels = c("Control", "Wound", "Gallery")))
expDes$grp <-
  with(expDes, factor(grp,
                      levels = paste(levels(gType),
                                     rep(levels(tx), each = 2), sep = ".")))
str(expDes) # 'data.frame':  24 obs. of  6 variables:
test_that("design matrix has 24 rows upon import",
          expect_equal(24, nrow(expDes)))

# Load counts into DGEList object from edgeR package.
y <- DGEList(counts = filteredSailfishCounts, group = expDes$grp)

# TMM Normalization by Depth
y <- calcNormFactors(y)

# make model matrix
#modMat <- model.matrix(~ tx/gType - 1, expDes)
modMat <- model.matrix(~ gType * tx, expDes)

# not sure why this is necessary but I do it to please makeContrasts
colnames(modMat) <- gsub(":", ".", colnames(modMat))

# voom transformation
v <- voom(y, modMat, plot = TRUE) # take a couple moments

# Linear modelling
fit <- lmFit(v, modMat)
fit2 <- eBayes(fit)

# get inferential info on each interaction effect on its own
tt_gTypeH898res.txWound <-
  topTable(fit2, coef = grep("gTypeH898res.txWound", colnames(modMat)))
head(tt_gTypeH898res.txWound)

tt_gTypeH898res.txGallery <-
  topTable(fit2, coef = grep("gTypeH898res.txGallery", colnames(modMat)))
head(tt_gTypeH898res.txGallery)

# create contrast matrix to take difference of interaction terms
cont_matrix <-
  makeContrasts(gTypeH898res.txGallery - gTypeH898res.txWound,
                levels = modMat)
fit3 <- contrasts.fit(fit, cont_matrix)
fit4 <- eBayes(fit3)
tt_diff_interactions <- topTable(fit4, coef = 1)
head(tt)

## I AM JUST ABOVE HERE
## STILL NOT SATISFIED WITH EXTRACTION OF INFERENCE RE: DIFFERENCE OF THE INTERACTION TERMS

## CODE BELOW HERE NOT FRESHED RECENTLY
## SOME WRITTEN BY ME, SOME FRAGMENTS LEFT FROM MACK THAT I MIGHT WANT TO REVIVE

data.wide <- cbind(expDes, t(y$counts[rownames(tt)[1:4], ]))
data.tall <- melt(data.wide,
                  id.vars = names(expDes),
                  variable.name = 'contig', value.name = 'Expression')
str(data.tall)

#x <- data.frame(expDes, gExp = y$counts[rownames(tt)[1], ])
p <- ggplot(data.tall, aes(x = gType, y = Expression))
p + geom_point() + facet_grid(tx ~ contig)


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


linear.results <- decideTests(fit2, method="global", 
                              adjust.method="fdr", p.value=0.01)

summary(linear.results)

###

vennDiagram(linear.results, include="both")
vennDiagram(linear.results, include="up")
vennDiagram(linear.results, include="down")

###
